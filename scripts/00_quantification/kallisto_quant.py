#!/usr/bin/env python3
"""Quantify the raw reads with kallisto and write the gene-level count tables.

This is the step that turns FASTQ files into the per-gene count matrices stored
in data/counts/. It is included so the analysis can be reproduced from the raw
data, but it cannot run without the reads and the reference, which are not part
of this repository (see data/README.md).

Usage
    python scripts/00_quantification/kallisto_quant.py \
        --raw     /path/to/fastq \
        --ref     /path/to/GCA_001048715.1 \
        --outdir  /path/to/kallisto_out \
        [--stranded rf|fr|none] [--threads 8] [--subsample 4000000]

Expected inputs
    --raw    one directory holding <sample>_1.fastq.gz and <sample>_2.fastq.gz
             for each of 3610-1, 3610-2, KO-1, KO-3
    --ref    the NCBI dataset directory for GCA_001048715.1, containing
             cds_from_genomic.fna and genomic.gff

Outputs
    <outdir>/<sample>/abundance.tsv    kallisto output per sample
    <outdir>/tx2gene.tsv               transcript -> locus_tag map, read from
                                       the [locus_tag=...] field of the CDS headers
    <outdir>/counts/htseq_<label>.tab  gene-level counts in the three-column
                                       format used by data/counts/ --
                                       locus_tag, gene length, count
    <outdir>/counts/fungal_read_fraction.csv
                                       the percentage of each library assigned
                                       to fungal transcripts, measured from this
                                       run; copy into data/annotation/ so the QC
                                       step uses it instead of the recorded values
    <outdir>/kallisto_index.idx
    <outdir>/run_params.json           the depth, strandedness, reference md5
                                       and FASTQ manifest this directory was
                                       built from; reused only if they still
                                       match

Strandedness. The libraries are TruSeq Stranded mRNA (dUTP), so read 2 is on
the sense strand and the correct kallisto setting is --rf-stranded; that is the
default here. Pass --stranded none to reproduce an unstranded run.

Depth. The analyses in this repository were run on a 4 million read-pair
subsample of each library because of a limitation of the computing environment
that was available; --subsample reproduces that. Omit it to use the full depth,
which is what a final analysis should do.

Caching. Both the index and the per-sample quantification are skipped when their
output files already exist, which makes an interrupted run cheap to resume but
would otherwise let a rerun silently keep results produced from different
inputs. To prevent that, <outdir>/run_params.json records everything the
directory was built from -- the subsample depth, the strandedness, the md5 of
the reference FASTA, and the size and modification time of every input FASTQ --
and the run stops if any of it has changed. A directory that already holds
quantification output but no run_params.json was built by an older version of
this script and is also refused, because there is no way to know what produced
it.
"""
import argparse, hashlib, json, os, platform, re, shutil, subprocess, sys, tarfile, urllib.request, zipfile
from collections import defaultdict
from pathlib import Path

# sequencing-facility sample name -> the label used in data/counts/
# the four libraries analysed here; see README.md on the two that are excluded
SAMPLES = {"3610-1": "WT1", "3610-2": "WT2", "KO-1": "KO1", "KO-3": "KO3"}
KVER = "0.50.1"
KURL = {
    "Windows": f"https://github.com/pachterlab/kallisto/releases/download/v{KVER}/kallisto_windows-v{KVER}.zip",
    "Linux":   f"https://github.com/pachterlab/kallisto/releases/download/v{KVER}/kallisto_linux-v{KVER}.tar.gz",
    "Darwin":  f"https://github.com/pachterlab/kallisto/releases/download/v{KVER}/kallisto_mac-v{KVER}.tar.gz",
}


def sh(cmd):
    print("\n$ " + " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([str(c) for c in cmd], check=True)


def need(path, what):
    if not Path(path).exists():
        sys.exit(f"[X] {what} not found:\n    {path}")


def find_kallisto(tools: Path):
    exe = "kallisto.exe" if platform.system() == "Windows" else "kallisto"
    found = shutil.which("kallisto") or shutil.which("kallisto.exe")
    if found:
        return found
    if tools.is_dir():
        for root, _, files in os.walk(tools):
            if exe in files:
                return str(Path(root) / exe)
    return None


def install_kallisto(tools: Path):
    k = find_kallisto(tools)
    if k:
        print("kallisto:", k)
        return k
    tools.mkdir(parents=True, exist_ok=True)
    url = KURL[platform.system()]
    pkg = tools / Path(url).name
    print(f"=== Downloading kallisto v{KVER} for {platform.system()} ===\n  {url}")
    urllib.request.urlretrieve(url, pkg)
    if pkg.suffix == ".zip":
        with zipfile.ZipFile(pkg) as z:
            z.extractall(tools)
    else:
        with tarfile.open(pkg) as t:
            t.extractall(tools)
    k = find_kallisto(tools)
    if not k:
        sys.exit(f"[X] Downloaded kallisto but could not find the binary under {tools}.")
    if platform.system() != "Windows":
        os.chmod(k, 0o755)
    print("kallisto:", k)
    return k


def md5(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def fastq_manifest(raw: Path) -> dict:
    """Size and modification time of every input FASTQ. Hashing 18 GB of reads
    on every run would cost minutes; size and mtime together are enough to
    notice that the inputs are not the ones an output directory was built
    from."""
    m = {}
    for s in SAMPLES:
        for r in ("1", "2"):
            p = raw / f"{s}_{r}.fastq.gz"
            st = p.stat()
            m[p.name] = {"bytes": st.st_size, "mtime": int(st.st_mtime)}
    return m


def check_params(outdir: Path, params: dict):
    """An output directory records what it was built from -- depth, strandedness,
    the md5 of the reference FASTA, and the size and mtime of every input FASTQ
    -- and refuses to be reused when any of that has changed. Without this, a
    full-depth rerun would silently keep cached subsampled results, and a
    changed reference would leave stale abundance files in place."""
    f = outdir / "run_params.json"
    if f.is_file():
        old = json.loads(f.read_text())
        differ = {k: (old.get(k), params[k]) for k in params if old.get(k) != params[k]}
        if differ:
            lines = []
            for k, (o, n) in differ.items():
                if isinstance(o, dict) or isinstance(n, dict):
                    o, n = (o or {}), (n or {})
                    changed = sorted(set(o) | set(n))
                    changed = [c for c in changed if o.get(c) != n.get(c)]
                    lines.append(f"    {k}: {len(changed)} file(s) differ, e.g. "
                                 f"{changed[0]}" if changed else f"    {k}: differs")
                else:
                    lines.append(f"    {k}: this directory was built with {o!r}, "
                                 f"you asked for {n!r}")
            sys.exit("[X] --outdir was built from different inputs:\n" + "\n".join(lines) +
                     "\n    Use a different --outdir, or delete this one and start again.")
    else:
        # No manifest. If the directory already holds quantification output it
        # came from an older version of this script, and nothing records what
        # produced it -- writing a manifest now would bless results that may
        # have come from different reads, a different reference or a different
        # depth. Refuse rather than guess.
        stale = [p.name for p in outdir.glob("*/abundance.tsv")]
        if stale or (outdir / "kallisto_index.idx").is_file():
            sys.exit(
                "[X] --outdir already contains quantification output but no "
                "run_params.json,\n    so there is no record of what produced "
                "it. Delete the directory and start\n    again, or point "
                "--outdir somewhere new.")
        outdir.mkdir(parents=True, exist_ok=True)
        f.write_text(json.dumps(params, indent=2))


def build_tx2gene(transcripts: Path, out: Path, tag: str):
    """Map each CDS id to its locus_tag, taken from the [locus_tag=...] field of
    the FASTA header. `tag` records which reference the transcript came from."""
    rows = []
    for line in transcripts.open():
        if not line.startswith(">"):
            continue
        tx = line[1:].split()[0]
        m = re.search(r"\[locus_tag=([^\]]+)\]", line)
        rows.append((tx, m.group(1) if m else tx, tag))
    if not rows:
        sys.exit(f"[X] No FASTA headers found in {transcripts}")
    with out.open("a") as o:
        for r in rows:
            o.write("\t".join(r) + "\n")
    print(f"  {len(rows):,} transcripts from {transcripts.name} ({tag})")
    return {tx for tx, _, _ in rows}


def aggregate(abundance: Path, tx2gene: dict, keep: set):
    """Sum kallisto's estimated counts over the transcripts of each gene.

    kallisto reports est_counts per transcript; a gene with several CDS entries
    therefore needs its transcripts added up. Counts are rounded at the end,
    once, so that the rounding error cannot accumulate across transcripts.

    Returns (counts, lengths). Both dictionaries are built fresh for this one
    abundance file. The gene length is the sum of its transcripts' target
    lengths, which kallisto reports identically for every sample -- so the
    lengths must be recomputed per sample and NOT accumulated across them, or
    the second sample would report twice the true length, the third three
    times, and so on.
    """
    counts = defaultdict(float)
    lengths: dict = {}
    with abundance.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        i_id, i_len, i_cnt = (header.index(c) for c in ("target_id", "length", "est_counts"))
        for line in f:
            p = line.rstrip("\n").split("\t")
            tx = p[i_id]
            if tx not in keep:
                continue
            g = tx2gene[tx]
            counts[g] += float(p[i_cnt])
            lengths[g] = lengths.get(g, 0) + int(p[i_len])
    return {g: int(round(c)) for g, c in counts.items()}, lengths


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", required=True, type=Path, help="directory of FASTQ files")
    ap.add_argument("--ref", required=True, type=Path,
                    help="NCBI dataset directory for GCA_001048715.1")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--stranded", choices=("rf", "fr", "none"), default="rf",
                    help="library strandedness (default rf, i.e. TruSeq dUTP)")
    ap.add_argument("--threads", type=int, default=max(2, os.cpu_count() or 4))
    ap.add_argument("--subsample", type=int, default=0,
                    help="use only the first N read pairs per sample (0 = full depth)")
    a = ap.parse_args()

    transcripts = a.ref / "cds_from_genomic.fna"
    gff = a.ref / "genomic.gff"
    need(transcripts, "CDS FASTA (cds_from_genomic.fna)")
    need(gff, "GFF (genomic.gff)")
    for s in SAMPLES:
        for r in ("1", "2"):
            need(a.raw / f"{s}_{r}.fastq.gz", f"FASTQ {s}_{r}")
    print(f"OK: {len(SAMPLES)} samples and the reference are present.")

    print("Fingerprinting the reference and the reads ...")
    check_params(a.outdir, {
        "subsample":     a.subsample,
        "stranded":      a.stranded,
        "reference":     transcripts.name,
        "reference_md5": md5(transcripts),
        "fastq":         fastq_manifest(a.raw),
    })
    kallisto = install_kallisto(a.outdir / "tools")
    idx = a.outdir / "kallisto_index.idx"

    print("=== 1. transcript -> gene map from the CDS headers ===")
    tx2gene_file = a.outdir / "tx2gene.tsv"
    if tx2gene_file.exists():
        tx2gene_file.unlink()
    tx2gene_file.write_text("transcript\tgene\tsource\n")
    fungal_tx = build_tx2gene(transcripts, tx2gene_file, "fungal")
    tx2gene = {}
    for line in tx2gene_file.open().readlines()[1:]:
        tx, g, _ = line.rstrip("\n").split("\t")
        tx2gene[tx] = g

    print("=== 2. kallisto index ===")
    if idx.is_file():
        print("  (index already present, skipping)")
    else:
        sh([kallisto, "index", "-i", idx, transcripts])

    print("=== 3. kallisto quant ===")
    strand_flag = {"rf": ["--rf-stranded"], "fr": ["--fr-stranded"], "none": []}[a.stranded]
    for s in SAMPLES:
        out = a.outdir / s
        if (out / "abundance.tsv").is_file():
            print(f"  {s}: already quantified, skipping")
            continue
        r1, r2 = a.raw / f"{s}_1.fastq.gz", a.raw / f"{s}_2.fastq.gz"
        if a.subsample:
            # keep the first N read pairs; 4 lines per read in FASTQ format
            import gzip
            sub = a.outdir / "subsampled"
            sub.mkdir(exist_ok=True)
            r1s, r2s = sub / f"{s}_1.fastq.gz", sub / f"{s}_2.fastq.gz"
            for src, dst in ((r1, r1s), (r2, r2s)):
                if dst.is_file():
                    continue
                print(f"  subsampling {src.name} -> {a.subsample:,} reads")
                with gzip.open(src, "rt") as fi, gzip.open(dst, "wt") as fo:
                    for i, line in enumerate(fi):
                        if i >= a.subsample * 4:
                            break
                        fo.write(line)
            r1, r2 = r1s, r2s
        sh([kallisto, "quant", "-i", idx, "-o", out, "-t", a.threads] + strand_flag + [r1, r2])

    print("=== 4. gene-level count tables ===")
    cdir = a.outdir / "counts"
    cdir.mkdir(exist_ok=True)
    reference_lengths = None
    fractions: dict = {}
    for s, label in SAMPLES.items():
        counts, lengths = aggregate(a.outdir / s / "abundance.tsv", tx2gene, fungal_tx)
        # every sample is quantified against the same index, so the gene lengths
        # must agree; if they do not, something is wrong with the inputs
        if reference_lengths is None:
            reference_lengths = lengths
        elif lengths != reference_lengths:
            sys.exit(f"[X] {s}: gene lengths differ from the first sample. "
                     "The abundance files were not produced from the same index.")
        dest = cdir / f"htseq_{label}.tab"
        with dest.open("w") as o:
            for g in sorted(counts):
                o.write(f"{g}\t{lengths[g]}\t{counts[g]}\n")
        print(f"  {dest.name}: {len(counts):,} genes, {sum(counts.values()):,} counts")

        # Fraction of the library that is fungal, measured from this run rather
        # than taken from a hard-coded table, so that the QC step describes the
        # quantification that actually produced these counts.
        info = json.loads((a.outdir / s / "run_info.json").read_text())
        processed = info.get("n_processed") or 0
        fractions[label] = round(100.0 * sum(counts.values()) / processed, 2) if processed else None

    frac_file = cdir / "fungal_read_fraction.csv"
    with frac_file.open("w") as o:
        o.write("sample,fungal_pct\n")
        for label, pct in fractions.items():
            o.write(f"{label},{'' if pct is None else pct}\n")
    print(f"  {frac_file.name}: " +
          ", ".join(f"{k} {v}%" for k, v in fractions.items()))

    print(f"\nDone. Copy {cdir}/*.tab into data/counts/ and "
          f"{frac_file.name} into data/annotation/ (overwriting the distributed "
          "files), then run scripts/01_differential_expression/run_de.R")


if __name__ == "__main__":
    main()

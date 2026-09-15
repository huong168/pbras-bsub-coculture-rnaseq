# Data

## What is in this repository

| Path | Contents |
|---|---|
| `counts/htseq_*.tab` | Per-gene read counts, one file per sample. Three columns: locus tag, gene length, count. These are the input to every analysis in `scripts/01_*` onwards. |
| `annotation/gene_length_Penbrz1.txt` | Gene lengths (sum of CDS lengths per locus tag), used for RPKM. |
| `annotation/gene_annot.tsv` | Locus tag, contig, coordinates, strand and product description for all 11,432 genes, parsed from the reference GFF. |
| `annotation/bgc_genes.csv` | Locus tags assigned to biosynthetic gene clusters (antiSMASH), with the cluster class. |
| `annotation/gene_sets.csv` | Gene sets used for the cluster-level tests: the BGC classes plus two manually curated sets, `Iron_uptake` and `Fe_S_iron_dependent`. |
| `annotation/fungal_read_fraction.csv` | Optional. Percentage of each library assigned to fungal transcripts. If present it replaces the values recorded in `scripts/config.R`, so a re-quantification updates the QC step automatically; `kallisto_quant.py` writes it. |
| `annotation/tx2gene.tsv` | Transcript to locus-tag map derived from the CDS FASTA headers. |

## What is not in this repository

**Raw reads.** All six paired-end libraries (12 FASTQ files, 18.1 GB), too large
to distribute here. That includes 3610-5 and KO-2, which are not analysed in this
repository. Metadata for depositing them in the NCBI SRA is prepared in
`submission/`; file sizes, md5 checksums and per-run statistics are in
`submission/sequencing_run_statistics.tsv`.

> BioProject accession: PRJNA1528400
> it in the Data Availability statement of the manuscript.

**Reference genome.** *Penicillium brasilianum* MG11, assembly
`GCA_001048715.1` (Pbras_Allpaths-LG), 11,432 genes. Download the NCBI dataset
and point `scripts/00_quantification/kallisto_quant.py` at the directory that
contains `cds_from_genomic.fna` and `genomic.gff`.

## Samples

| File | Sample | Group | Co-culture partner | Reads mapping to the fungal genome |
|---|---|---|---|---|
| `htseq_WT1.tab` | 3610-1 | treatment | *B. subtilis* DK1042 wild type | 57.8 % |
| `htseq_WT2.tab` | 3610-2 | treatment | *B. subtilis* DK1042 wild type | 77.8 % |
| `htseq_KO1.tab` | KO-1 | control | *B. subtilis* DK1042 Δ*dhbF* | 19.3 % |
| `htseq_KO3.tab` | KO-3 | control | *B. subtilis* DK1042 Δ*dhbF* | 18.6 % |

The library names `3610-1` and `3610-2` carry the number of the parental strain,
NCIB 3610, because that is how the sequencing facility labelled them; the
bacterium actually used is its naturally competent derivative **DK1042**
(*comI*<sup>Q12L</sup>), and the Δ*dhbF* mutant is in the same DK1042 background.
The FASTQ file names in `submission/` keep the facility labels so that they match
the files as delivered.

Six libraries were sequenced, three per condition. Two — **3610-5** and
**KO-2** — are not analysed here; `README.md` gives the reason and the numbers.
Both are deposited in the SRA along with the four used here, so the full set is
available to anyone who wants to repeat the analysis differently.

The RNA is fungal in both conditions; the only difference between the groups is
the bacterial partner. Δ*dhbF* lacks the NRPS that assembles bacillibactin but
retains the upstream *dhbACEB* genes, so it still produces the precursor
2,3-dihydroxybenzoate. The contrast therefore isolates bacillibactin itself.

The last column matters and is analysed in
`scripts/03_qc/qc_bacterial_load.R`: the libraries are mixed fungal/bacterial,
and the fungal share differs between the groups.

## Library preparation and sequencing

RNA was treated with the RNase-Free DNase I Set (Qiagen) and purified with the
RNA Clean & Concentrator-25 kit (Zymo Research); integrity was checked on an
Agilent 2100 Bioanalyzer, with RIN 9.1–10.0 across all samples.
Poly(A)-selected, strand-specific libraries prepared with the Illumina TruSeq
Stranded mRNA Library Prep Kit (TruSeq Stranded mRNA Reference Guide
#1000000040498 v00, dUTP protocol, so the libraries are reverse-stranded).
Paired-end 101 bp sequencing on an Illumina NovaSeq X at Macrogen Inc., Seoul (order HN00283426,
August 2026). 31.8–38.2 million read pairs per sample, Q30 94.0–94.7 %.

## Sequencing depth

The counts distributed here were generated from a **4 million read-pair
subsample** of each library (about 11 % of the sequenced data), because of a
limitation of the computing environment in which the analysis was carried out.
Re-running `kallisto_quant.py` without `--subsample` reproduces the whole
analysis at full depth; the script refuses to reuse an output directory that was
built at a different depth, so a full-depth rerun cannot silently pick up the
subsampled results.

What full depth will do is not certain in advance. More reads reduce the
sampling noise in each count, which usually lowers the *P* values of genes whose
effect is real — but the fold-change estimates for lowly expressed genes can
move as well, in either direction, and some genes excluded by the expression
filter at this depth will enter the analysis at full depth. The one thing that
should not change is the direction of the *ftm* cluster, which is supported by
eight genes moving together.

One further property of the count tables is worth knowing:

- They were produced **without** `--rf-stranded`, that is, treating the libraries
  as unstranded, even though the TruSeq dUTP protocol makes them
  reverse-stranded. `kallisto_quant.py` now passes `--rf-stranded` by default.
  This matters mainly for overlapping and antisense transcripts.

This was not changed retroactively, because the counts distributed here are the
ones the published analysis used.

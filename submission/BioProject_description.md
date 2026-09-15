# BioProject — text to paste into the NCBI submission portal

Register the BioProject **first**: SRA and BioSample both need its accession.

---

## Project title

> Transcriptional response of *Penicillium brasilianum* to bacillibactin during
> co-culture with *Bacillus subtilis*

## Public description

> *Penicillium brasilianum* SFC20200506-M13 was co-cultured with *Bacillus subtilis*
> DK1042 (the naturally competent comI-Q12L derivative of NCIB 3610)
> or with an isogenic Δ*dhbF* mutant in the same background, which cannot assemble the
> catecholate siderophore bacillibactin but retains the upstream *dhbACEB* genes
> and therefore still produces the precursor 2,3-dihydroxybenzoate. The bacterial
> partner is the only difference between the two conditions, so the comparison
> isolates the effect of bacillibactin itself on the fungus.
>
> Co-cultures were grown on potato dextrose agar. Total RNA was sequenced from
> three biological replicates of each co-culture.
> Reads were mapped to the *P. brasilianum* MG11 reference genome (GCA_001048715.1) to
> measure the fungal transcriptome. The genes of the fumitremorgin/verruculogen
> biosynthetic cluster (*ftmA*–*ftmH*) are the pre-specified focus of the study.

## Relevance

`Environmental` (or `Agricultural` — choose whichever your journal or funder
prefers; this field is not used for retrieval)

## Project data type

`Transcriptome or Gene expression`

## Sample scope

`Multiisolate`

## Organism

`Penicillium brasilianum` (taxid 104259), strain **SFC20200506-M13** — the same
fungal isolate as in the laboratory's earlier BioProject PRJNA1146190. See the
note in `README_SRA_submission.md` about how the co-culture is declared.

## Relationship to the existing BioProject

PRJNA1146190 holds the earlier *B. subtilis* IAM 1145 / 168 Δ*dhbF* co-culture
RNA-seq and is already cited in a published preprint. This submission is a
separate experiment with a different bacterial background (DK1042 WT vs Δ*dhbF*)
and should therefore be registered as its own BioProject; ask
`bioprojecthelp@ncbi.nlm.nih.gov` to group the two under an umbrella BioProject
if you want them linked.

## Grant / funding

TODO — grant number, agency and title, if the work is funded. This is optional
in the portal but journals increasingly ask for it.

## Release date

Choose **hold until publication** and give a date you can extend. The data
become public automatically on that date, or when the accession is cited in a
published paper, whichever comes first.

---

## Linked submissions

| Object | When | Accession goes into |
|---|---|---|
| BioProject | first | `SRA_metadata.tsv` → `bioproject_accession` |
| BioSample × 6 | second | `SRA_metadata.tsv` → `biosample_accession` |
| SRA runs × 6 | last | the paper, and `data/README.md` in this repository |

After the runs are accepted, add the accessions to `data/README.md` and to the
Data Availability statement of the manuscript.

# NCBI SRA submission metadata

The metadata accompanying the deposit of the raw reads. These files describe
**all six sequenced libraries**, including 3610-5 and KO-2, which are not
analysed in this repository — see the root `README.md` for why they were
excluded and `data/README.md` for the library details.

| File | Contents |
|---|---|
| `BioProject_description.md` | The text registered with the BioProject. |
| `BioSample_attributes_Microbe.1.0.tsv` | One row per biological sample, in the NCBI Microbe 1.0 package format. |
| `SRA_metadata.tsv` | One row per sequencing run: library strategy, instrument, design description and FASTQ file names. |
| `sequencing_run_statistics.tsv` | Read pairs, bases, Q20, Q30 and GC content per library, from the sequencing facility report. |
| `md5_checksums.txt` | Checksums of the twelve FASTQ files as delivered. |


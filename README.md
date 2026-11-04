# Gastric MALT lymphoma analysis

R scripts for WES, WTS and H. pylori quantitation in gastric MALT lymphoma.

| File | Analysis |
| --- | --- |
| [WES.R](WES.R) | ANNOVAR/VEP annotations, genotype filters, variant tables and SNV/indel oncoplot |
| [WTS.R](WTS.R) | Raw-count analysis with DESeq2, DEG tables and GO/KEGG GSEA |
| [Microbiome.R](Microbiome.R) | H. pylori read counts, relative abundance and Z-scores from an existing taxonomic table |
| [METHODS.md](METHODS.md) | Input columns, analysis settings, source-code correspondence and remaining study questions |

## Running an analysis

1. Download this repository and open the required R script in RStudio.
2. Set `Project.dir` to the local repository directory. Input paths in `SET` are relative to this directory.
3. Place the input TSV files under `data/`, following the columns described in [METHODS.md](METHODS.md).
4. Edit the settings at the top of the script. Fill required `NULL` settings with confirmed analysis choices; WTS comparisons are in `VAR`, and enrichment settings are in `GO` and `KEGG`.
5. Click Source once. Each script contains the functions needed for that module.

No separate setup script or configuration file is needed. With `Install.packages <- TRUE`, the script installs missing packages from CRAN or Bioconductor without updating installed packages. Package installation requires repository access and a compatible R/Bioconductor environment. Set it to `FALSE` to use an already prepared environment. Routine progress messages are suppressed; errors and calculation warnings remain visible.

## Requirements and inputs

| Module | R packages used directly, in addition to base R | Inputs |
| --- | --- | --- |
| All | jsonlite, digest, dplyr | Tab-delimited files with explicit identifiers and metadata |
| WES | maftools when MAF/oncoplot output is enabled | ANNOVAR and VEP tables, sample-role metadata; verified case map when reporting case aliases |
| WTS | DESeq2; clusterProfiler for GSEA; the selected OrgDb for GO source-API mode | Raw integer counts, sample/patient metadata, local gene-ID mapping and stated reference settings |
| Microbiome | No additional packages beyond the common packages | Species-level taxon table, sample/metric column map, metadata; explicit denominator/reference samples for optional recomputation |

The validation environment used R 4.3.3, DESeq2 1.42.0, maftools 2.26.1 and dplyr 1.2.1. These are measured validation versions, not the unknown original study environment. A historical package lockfile is not available.

## Outputs

Each execution creates `results/<module>/<Run.ID>/`. Existing run directories are not overwritten.

| Module | Main outputs |
| --- | --- |
| WES | Variant selection audit, selected variants, per-sample counts, optional MAF and oncoplot |
| WTS | Complete/selected/upregulated/downregulated DEG tables, model/contrast records, optional GO/KEGG tables and ranking files |
| Microbiome | Imported target metrics, H. pylori quantitative summary and quantitation provenance |

Every run also saves `configuration_used.json`, `input_hashes.tsv`, `run_status.json` and `sessionInfo.txt`.

## Relationship to the reported results

WES methods relate to Supplementary Table S13. WTS methods relate to Table 2 and Supplementary Tables S1-S12. Table 1 contains clinical characteristics. Neither supplied final workbook contains a quantitative microbiome table.

The original analysis objects and processing order are retained where supported by the source. Required metadata and method choices are explicit: the final WES filter profile, case aliases, WTS pairing/model/direction and database snapshots are not established by the supplied code alone. Missing required information stops execution. The S4 raw/adjusted p-value relationship and the source/report KEGG adjustment discrepancy remain unresolved; values and columns have not been exchanged.

## Verification and data availability

The three scripts passed R parsing and synthetic-input checks in the preceding code revision. This packaging revision leaves all three script files unchanged. GO/KEGG statistical-engine execution was not performed because clusterProfiler was unavailable. Actual study execution and numerical reproduction of the final tables have not been performed: study inputs and verified mappings were not supplied.

Original source files, final workbooks, specimen-level data, synthetic fixtures and development checks are kept in separate review material. Local inputs and results are excluded by `.gitignore`; review any local script settings before committing them. The public repository consists of the three scripts, this README, METHODS.md and `.gitignore`.

License, author/citation details and publication identifiers require confirmation by the owner. No repository URL, release or DOI has been assigned.

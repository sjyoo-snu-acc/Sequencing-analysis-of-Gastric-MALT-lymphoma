# Sequencing analysis of gastric MALT lymphoma

R scripts for downstream whole-exome sequencing (WES), whole-transcriptome sequencing (WTS) and microbiome analyses.

| File | Description |
| --- | --- |
| `WES.R` | ANNOVAR/VEP annotation merge, genotype parsing, variant filtering, MAF conversion and oncoplots |
| `WTS.R` | DESeq2 comparisons, volcano plots, heatmaps and GO/KEGG enrichment |
| `Microbiome.R` | Genus composition, alpha diversity and linear regression |
| `METHODS.md` | Input formats and analysis settings |

## Requirements

Use R with a compatible Bioconductor release. Install the required packages before running the scripts.

```r
install.packages(c("dplyr", "ggplot2", "stringr", "reshape2", "RColorBrewer",
                   "ggrepel", "matrixStats", "pheatmap", "writexl"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "DESeq2", "maftools", "clusterProfiler",
                      "AnnotationDbi", "org.Hs.eg.db", "GO.db", "BiocParallel",
                      "enrichplot"), update = FALSE, ask = FALSE)
```

`dplyr >= 1.1.1` is required. The optional live Ensembl mapping also requires `biomaRt`. The manuscript lists R 4.2.2 and DESeq2 1.38.3; the versions used for each execution are recorded in its `sessionInfo.txt`.

## Running the analysis

1. Set the R working directory to the repository root.
2. Prepare the input files listed below and edit `Input.dir` or the named file paths at the beginning of the selected script.
3. Run `source("WES.R")`, `source("WTS.R")` or `source("Microbiome.R")`.

Alternatively, run `Rscript WES.R`, `Rscript WTS.R` or `Rscript Microbiome.R` from the repository root.

Each execution writes to `results/<module>/<Run.ID>/`. `Run.ID` defaults to a timestamp. An existing output directory stops execution; select a new run ID to run again. The scripts do not install or update packages during analysis.

## Input files

All inputs are tab-delimited text files with headers. Column names and sample identifiers are case-sensitive.

| Directory | Required files |
| --- | --- |
| `data/WES/` | `R.GML_H.pylori_n10_Mutect2_paired_merged.AnnoVar.hg19_multianno.txt`, `R.GML_H.pylori_n10_Mutect2_paired_merged.VEP.result.txt`, `sample_info.tsv` |
| `data/WTS/` | `H.pylori-Gastric_MALToma-WTS-COUNT.tsv`, `H.pylori-Gastric_MALToma-WTS-TPM.tsv`, `grouptable_new.tsv`, `gene_symbols.tsv` |
| `data/Microbiome/` | `genus_counts.tsv`, `genus_relative_percent.tsv`, `alpha_diversity.tsv`, `grouptable_new.tsv` |

**WES:** `sample_info.tsv` contains `sample_column`, `case_number`, `include_wes` (`TRUE/FALSE`) and `hp_status` (`positive/negative/unknown`). Optional `include_plot` specifies the plotting cohort. Optional `declared_fusions.tsv` contains `case_number`, `gene_symbol` and provenance for externally determined fusion annotations.

**WTS:** COUNT and TPM use `gene_id`, one empty annotation column, and then the same sample columns in the same order. COUNT contains nonnegative RSEM expected counts. Metadata contains `Name`, `sample` (`Normal/Tumor/Eradicated`) and `HP` (`positive/negative`). For the configured sample-name convention, `01-6001-N-WTS-1` maps to `N1_016001`. `gene_symbols.tsv` contains `gene_id` and `symbol`; record the annotation origin and release in `Annotation.label`.

**Microbiome:** Each matrix begins with a `genus` column followed by sample columns. The two matrices contain the same genus and sample identifiers. Supply untransformed reads in `genus_counts.tsv` and percentages from 0 to 100 in `genus_relative_percent.tsv`. Include the four genus rows `Helicobacter`, `Rothia`, `Veillonella` and `Porphyromonas`, with explicit zeros for absence. `alpha_diversity.tsv` contains `sample`, `Shannon` and `Simpson_index`; metadata follows the WTS schema. Record the upstream classification database, abundance denominator and diversity-index definition in the input documentation and `Input.description`.

## Analysis settings

WES uses DP >= 50, alternate depth >= 5, strand ratio 0.25–4, VAF >= 0.05 and < 0.40, together with the annotation and COSMIC filters in `WES.R`. `AF.column` selects the population-frequency annotation. The default `AF.cutoff=0.1` is an allele fraction, and missing annotations are retained with `Keep.missing.AF=TRUE`.

WTS uses `abs(log2FoldChange) >= 1` and `padj <= 0.05` for up/down lists and volcano colors. `Count.mode="truncate"` applies `as.integer` to expected counts; `"integer"` requires integer input. The configured DEG analysis excludes identifier 10 and follow-up specimens; these selections are recorded in `all_sample_metadata.tsv`.

| Numerator / denominator | Model |
| --- | --- |
| PNormal / NNormal | `~ Type` |
| PTumor / PNormal | `~ ID + Type` |
| PTumor / NTumor | `~ Type` |
| NTumor / NNormal | `~ Type` |

`Run.GSEA=TRUE` runs GO and KEGG enrichment. GO uses positive ranks of LFC with `scoreType="pos"`; KEGG uses signed LFC. The seed is 48575 and the backend is serial. Reference snapshots can be reused with `Use.reference.cache=TRUE`; see METHODS.md. Optional Hallmark analysis is enabled with `Run.Hallmark=TRUE` and a local `hallmark_term2gene.tsv` containing `gs_name` and `human_gene_symbol`.

Microbiome imports alpha-diversity values by default and uses genus percentages as regression predictors. `Alpha.mode="calculate"` instead calculates genus-level Shannon using natural logs and Simpson as `1-sum(p^2)`. `Relative.mode="sum_genus_counts"` derives percentages from the supplied genus-count total. `Predictor.metric="reads"` selects raw read counts instead of percentages. These settings define different calculations and should match the intended input preparation.

## Outputs

WES writes selected variants, per-case counts, MAF data and oncoplots. `BeforeFiltering.tsv` is the selected-record intermediate used for ANNOVAR-to-MAF conversion. Fusion annotations are recorded separately from SNV/indel calls.

WTS writes full/up/down DEG tables, result workbooks, volcano plots and TPM heatmaps. With GSEA enabled, it also saves enrichment tables, ranking vectors, result objects and reference snapshots. `Make.plots=FALSE` skips plotting; `Run.GSEA=FALSE` skips enrichment.

Microbiome writes genus composition and diversity plots, regression inputs, coefficients with 95% confidence intervals, and fitted model objects.

Completed runs also save input file hashes, settings, sample selections and `sessionInfo.txt`. WTS records the gene mapping and comparison design matrices. Repeating a calculation requires the same processed inputs, specimen mapping, settings, annotation references and package versions.

## Analysis scope and data

The scripts begin with processed annotation, expression and microbial tables. Read alignment, variant calling, fusion calling and microbial classification are upstream steps. The sequence accession given in the manuscript is **PRJNA1088908**. Study-specific input tables and specimen mappings are supplied separately from the code.

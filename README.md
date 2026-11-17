# Sequencing analysis of gastric MALT lymphoma

R scripts for downstream WES, WTS and microbiome analyses accompanying the manuscript. `WES.R` follows `Trimmedddd.R`; `WTS.R` retains `ENSG2SYM`, `degsea`, `pheatcustom` and the four original comparison calls from `250326.R`. `Microbiome.R` makes the original genus/diversity/regression analysis independently executable.

| File | Analysis |
| --- | --- |
| `WES.R` | ANNOVAR/VEP merge, genotype parsing, somatic variant filtering, MAF and oncoplots |
| `WTS.R` | Expected-count processing, DESeq2, volcano plots, heatmaps and optional GO/KEGG/Hallmark GSEA |
| `Microbiome.R` | Genus composition, imported or explicitly recalculated alpha diversity, univariate regressions and the original four-genus Shannon model |
| `METHODS.md` | Input schemas, retained analysis choices, changes and limitations |
| `VALIDATION.md` | Checks completed and remaining manuscript-result differences |
| `examples/` | Artificial data demonstrating the input formats; no study specimens |

## Installation

Use R with a compatible Bioconductor release. Install packages before running; the analysis scripts do not install or update packages.

```r
install.packages(c("dplyr", "ggplot2", "stringr", "reshape2", "RColorBrewer",
                   "ggrepel", "matrixStats", "pheatmap", "writexl"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("S4Vectors", "DESeq2", "maftools", "clusterProfiler",
                      "AnnotationDbi", "org.Hs.eg.db", "GO.db", "BiocParallel",
                      "enrichplot"), update = FALSE, ask = FALSE)
```

`dplyr >= 1.1.1` is required for the explicit annotation join. `biomaRt` is needed only if the optional live Ensembl mapping is enabled. A local mapping file is the default.

The manuscript reports R 4.2.2 and DESeq2 1.38.3. The original full environment and database snapshots were not supplied. This revision was checked in R 4.3.3 / DESeq2 1.42.0; those versions describe the new validation environment, not the historical analysis. See `VALIDATION.md` for the checks actually completed.

## Running the scripts

1. Download the repository and set the R working directory to its root.
2. Prepare the corresponding private input directory below. Edit the named paths/settings at the start of the chosen script where necessary.
3. Run `source("WES.R")`, `source("WTS.R")` or `source("Microbiome.R")`.

The equivalent command-line calls are `Rscript WES.R`, `Rscript WTS.R` and `Rscript Microbiome.R`, from the repository directory.

Each script creates `results/<module>/<Run.ID>/`, with a timestamp as the default run ID. An existing run directory stops execution; choose a new `Run.ID` after a failed or completed run. Outputs include the actual input file MD5 hashes, settings, sample selections and `sessionInfo.txt` on successful completion. WTS additionally records gene mappings and comparison metadata/design matrices. There is no separate `SET` configuration object or automatic package installer.

## Inputs

Input column names and sample IDs are case-sensitive. Supply tab-delimited files with headers and no row-number column. Study-specific mappings and patient/specimen data are not included in the repository.

| Directory | Required files |
| --- | --- |
| `data/WES/` | `R.GML_H.pylori_n10_Mutect2_paired_merged.AnnoVar.hg19_multianno.txt`, `R.GML_H.pylori_n10_Mutect2_paired_merged.VEP.result.txt`, `sample_info.tsv` |
| `data/WTS/` | `H.pylori-Gastric_MALToma-WTS-COUNT.tsv`, `H.pylori-Gastric_MALToma-WTS-TPM.tsv`, `grouptable_new.tsv`, `gene_symbols.tsv` |
| `data/Microbiome/` | `genus_counts.tsv`, `genus_relative_percent.tsv`, `alpha_diversity.tsv`, `grouptable_new.tsv` |

**WES.** `sample_info.tsv` contains `sample_column`, `case_number`, `include_wes` (`TRUE/FALSE`), `hp_status` (`positive/negative/unknown`), and optionally `include_plot` (`TRUE/FALSE`). List the genotype columns to be parsed, including controls marked `include_wes=FALSE`. Included cases require an identifier; a fusion-only plotting case can have an empty `sample_column`. Optional `declared_fusions.tsv` contains `case_number`, `gene_symbol` and, preferably, provenance. Declared fusions are external annotations and are never converted into fictitious SNV/indel calls.

**WTS.** COUNT and TPM retain the source layout: `gene_id`, one empty annotation column, then identical sample columns in identical order. COUNT contains nonnegative RSEM expected counts, not TPM. The original `as.integer` truncation is retained explicitly as `Count.mode="truncate"`; `"integer"` instead requires integer input. This is not a tximport/length-offset workflow. Metadata columns are `Name`, `sample` (`Normal/Tumor/Eradicated`) and `HP` (`positive/negative`). The five-part source name `01-6001-N-WTS-1` maps to metadata name `N1_016001`. `gene_symbols.tsv` contains `gene_id` and `symbol`; document its database, release and download date in `Annotation.label`. See the synthetic files for the exact layout.

**Microbiome.** The first column of each matrix is `genus`; subsequent columns are exact metadata sample IDs. Rows are unique genus names, including `Helicobacter`, `Rothia`, `Veillonella` and `Porphyromonas`. The read matrix contains untransformed counts; the relative-abundance matrix contains **percentages from 0 to 100**, with the same genera and samples. An explicit zero is required for a genuinely absent genus. Do not use species-level *H. pylori* rows in place of genus-level *Helicobacter*. `alpha_diversity.tsv` contains `sample`, `Shannon`, `Simpson_index`; `grouptable_new.tsv` uses the WTS metadata schema. Matrix axes and alpha rows can be reordered: alignment is by identifiers.

For an original wide Pavian export, explicitly select its genus rows and sample read/percentage columns to prepare these two matrices. Do not infer sample/metric identity from four-column blocks or truncate sample names. Record the taxonomic database, upstream tool/version, confidence setting, abundance denominator, and alpha-diversity calculation/rank in `Input.description` and your data documentation.

## Retained settings and corrections

- WES retains DP >= 50, alternate depth >= 5, strand ratio 0.25–4, VAF >= 0.05 and < 0.40, and the original annotation/COSMIC/consequence filters. Population AF is compared numerically. `AF.column="AF_eas"` explicitly selects the first duplicated ANNOVAR AF block; the second is `AF_eas_2`. The 0.1 fraction threshold and retention of missing `.` annotations are unchanged. Confirm the database identity/release of the selected block before describing it as a specific gnomAD release.
- WTS uses `abs(log2FoldChange) >= 1` and `padj <= 0.05` consistently for up/down selection and volcano colors. This corrects the source's use of raw p-values for colors/selection while plotting adjusted p-values. DESeq2 raw p-values and adjusted values are both retained in the output.
- WTS retains the source exclusion `Exclude.ID=c("10")`, excludes follow-up samples from DEG analysis, and keeps the original comparison models below. The reason for excluding identifier 10 was not supplied; it is not labelled a confirmed QC decision. Patient IDs in source filenames must be reconciled with the manuscript's case IDs privately.
- Microbiome includes its own supplied samples, including follow-ups. It does not inherit WTS exclusions. Its default regressions use the source's percentage predictor (`Predictor.metric="relative_percent"`), not raw read counts. The manuscript's wording/units require checking against the original microbial data.

| WTS numerator / denominator | Model |
| --- | --- |
| PNormal / NNormal | `~ Type` |
| PTumor / PNormal | `~ ID + Type` |
| PTumor / NTumor | `~ Type` |
| NTumor / NNormal | `~ Type` |

`Run.GSEA=TRUE` runs GO and KEGG. GO retains the original positive LFC-rank score with seeded random ties and `scoreType="pos"`; KEGG retains signed LFC scores. These source choices are documented, not replaced with a new enrichment method. The source seed 48575 is reset before each comparison's GSEA so optional plots do not alter random draws; the backend is serial. Local references can be reused with `Use.reference.cache=TRUE`; see METHODS.md. `Run.Hallmark=FALSE` by default; enabling the retained optional Hallmark calculation requires a named local `hallmark_term2gene.tsv` with `gs_name` and `human_gene_symbol`, and its release provenance.

Microbiome defaults to importing the original alpha values. `Alpha.mode="calculate"` explicitly recalculates genus-level Shannon with natural logs and Simpson as `1-sum(p^2)`. `Relative.mode="sum_genus_counts"` explicitly derives percentages from the supplied genus-count total. These alternatives are not claimed to be the historical methods. The original multivariable formula remains `Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas`; four separate genus predictors are also fitted univariately to each diversity index. No unreported group-comparison test or clinical positivity cutoff is added.

## Outputs and verification limits

WES writes `selected_variants.tsv`, per-case counts, the source-named `BeforeFiltering.tsv` (the **post-filter** ANNOVAR-to-MAF input), a SNV/indel MAF and oncoplot, and a separate declared-event table/plot. Plot cohort size comes from metadata, including specified mutation-negative cases.

WTS writes full/up/down DEG TSVs, full DEG workbooks, volcano plots and TPM heatmaps. With GSEA enabled it also saves enrichment tables, ranking vectors, result objects, reference snapshots and package/database information. `Make.plots=FALSE` skips figures, not statistical fitting. `Run.GSEA=FALSE` skips enrichment and must not be described as a full WTS/GSEA run.

Microbiome writes genus composition and diversity figures, keyed regression inputs, model objects, coefficients with 95% confidence intervals, and Figure 3 analysis plots. These are generated from the supplied data; the manuscript's coefficients are not hardcoded or sign-corrected.

The scripts start from processed annotations, expression matrices and microbial tables. FASTQ processing, Mutect2, fusion callers/intersection, Kraken2 and upstream alpha-diversity tools are not executed here. The manuscript lists sequence accession **PRJNA1088908**; reproducing downstream results additionally requires the processed inputs, specimen mapping and reference versions described above.

The code revision and computational checks do not establish exact numerical reproduction of every published table or figure. Actual microbiome inputs and the historical annotation/environment were not supplied; remaining result discrepancies are recorded in `VALIDATION.md`. No patient data or study result workbooks are bundled.

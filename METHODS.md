# Inputs and analysis methods

## WES

ANNOVAR and VEP tables are joined by `ID`. Multiple VEP transcripts per variant are supported. The ANNOVAR table contains genomic coordinates, refGene/ensGene function and amino-acid annotations, `CLNSIG_220425`, `cosmic96_coding`, VCF fields and genotype columns. VEP provides `CANONICAL`, `BIOTYPE`, `Protein_position` and `Amino_acids`.

`sample_info.tsv` identifies genotype columns and case aliases. FORMAT fields are parsed by tag name; `GT`, `AD`, `DP` and `SB` are required. AD contains reference/alternate counts and SB contains four strand counts. Multiallelic records must be normalized before use.

| Filter | Setting |
| --- | --- |
| Basic annotation | Exclude the intronic/intergenic/upstream/downstream/UTR/noncoding categories listed in WES.R |
| Genotype | Allele 1 present |
| Depth | DP >= 50; alternate AD >= 5 |
| Strand ratio | `(ref_fw + alt_fw)/(ref_rv + alt_rv)` between 0.25 and 4 |
| Transcript | `CANONICAL == "YES"`; `BIOTYPE == "protein_coding"` |
| COSMIC | Haematopoietic/lymphoid tissue or stomach annotation |
| Consequence | Nonsynonymous SNV, frameshift insertion/deletion, stopgain, startloss or stoploss |
| Caller filter | Exclude `normal_artifact` |
| Population AF | Numeric fraction < 0.1; missing values retained when `Keep.missing.AF=TRUE` |
| Tumor VAF | Alternate AD / total AD >= 0.05 and < 0.40 |
| Duplicate records | First qualifying sample/variant record |

`AF.column="AF_eas"` selects the first ANNOVAR frequency block. A second duplicate header is named `AF_eas_2`. Record the database and release corresponding to the selected block. AF is compared numerically, including scientific notation. The cutoff 0.1 is a fraction (10%).

`include_wes` selects cases for variant output. `include_plot` selects the plotting cohort and defaults to `include_wes`. Mutation-negative cases retain zero counts. A separate fusion-only plotting case may have an empty `sample_column`.

The MAF contains SNV/indel calls. `declared_fusions.tsv` provides external fusion annotations for a separate combined gene-event plot. The script performs no fusion calling. `BeforeFiltering.tsv` contains selected records used for MAF conversion.

## WTS

COUNT and TPM share the layout `gene_id`, one empty annotation column, then identically ordered sample columns. `Count.mode="truncate"` converts expected counts with `as.integer`, followed by a per-comparison median-count filter >= 10. `Count.mode="integer"` requires integer input. TPM is used for PCA and heatmaps; expected counts are used for DESeq2.

Metadata columns are `Name`, `sample` and `HP`. The five-part input naming convention is converted to metadata names as described in README.md. `ID` is the two-character identifier in positions 4–5 of `Name`. The configured DEG analysis excludes ID 10 and follow-up specimens. Inclusion and exclusion flags are recorded for every input sample.

| Numerator | Denominator | Design |
| --- | --- | --- |
| PNormal | NNormal | `~ Type` |
| PTumor | PNormal | `~ ID + Type` |
| PTumor | NTumor | `~ Type` |
| NTumor | NNormal | `~ Type` |

Models are fitted separately to the selected samples. Matrices are aligned by exact sample names. The patient-adjusted model requires a full-rank design and unique ID/type observations; incomplete tissue pairs can remain in that model. Other configured comparisons use `~ Type`.

`ENSG2SYM` reads a local `gene_symbols.tsv` snapshot by default. Matching removes optional ENSG version suffixes. The first mapping per ENSG and first input gene per duplicated symbol are retained; empty/unmapped symbols are excluded. The mapping used and a reusable two-column snapshot are exported. Optional live lookup requires `Query.Ensembl=TRUE`; `Ensembl.version` selects a release.

DESeq2 uses `results(..., contrast=c("Type", numerator, denominator), pAdjustMethod="fdr")` and orders results by padj. Up/down classification requires LFC >= 1 or <= -1 with padj <= 0.05. Raw p-values and adjusted values remain separate columns. For plotting, zero padj is capped before the logarithm; result-table values are preserved.

`degsea` performs the comparison and enrichment calculations. `pheatcustom` uses the fitted comparison's metadata, genes in padj order, log2(TPM+1), row scaling and sample annotations. Zero-variance rows with undefined scaling are omitted from the heatmap only.

### GO, KEGG and reference snapshots

With `Run.GSEA=TRUE`, GO BP and GO ALL use positive LFC ranks, seeded random tie ranks, `minGSSize=3`, BH adjustment, `scoreType="pos"` and `eps=0`. This is a positive-enrichment calculation. KEGG uses signed LFC following SYMBOL-to-ENTREZID conversion, the first mapping per symbol, `minGSSize=3`, BH adjustment and `eps=1e-10`.

The seed 48575 is reset before each comparison's GSEA and the backend is serial. Ranking vectors and seed states are saved. Repeating a calculation also requires matching input order, references and package versions.

`Use.reference.cache=FALSE` uses installed GO/OrgDb annotations and the current KEGG service. Result objects and snapshots of the tested gene sets are saved in the run's `GSEA_references/` directory, separately for each comparison.

To reuse those snapshots, set `Reference.dir` to the previous directory and `Use.reference.cache=TRUE`. Restore the recorded package versions and keep the same tested gene universe. With an empty `Explanation`, snapshot names are `GO_BP_<numerator>_<denominator>.rds`, `GO_ALL_<numerator>_<denominator>.rds` and `KEGG_<numerator>_<denominator>.rds`. Whole-reference files named `GO_BP.rds`, `GO_ALL.rds` or `KEGG.rds` are also accepted. Cache objects contain `term2gene`, `term2name`, `package_versions` and `eps`; tested-universe caches also contain `input_genes`. The script checks versions and records reference hashes.

Optional Hallmark analysis uses `Run.Hallmark=TRUE` and a local `hallmark_term2gene.tsv` with `gs_name` and `human_gene_symbol`. Its annotation release should be recorded with the input data.

## Microbiome

Microbiome.R runs independently from WTS.R.

| Input | Columns |
| --- | --- |
| `genus_counts.tsv` | `genus`, followed by sample read-count columns |
| `genus_relative_percent.tsv` | The same genus/sample identifiers, with percentages in [0,100] |
| `alpha_diversity.tsv` | `sample`, `Shannon`, `Simpson_index` |
| `grouptable_new.tsv` | `Name`, `sample` (Normal/Tumor/Eradicated), `HP` (positive/negative) |

Genus rows include Helicobacter, Rothia, Veillonella and Porphyromonas. Reads are nonnegative integers. A missing observation is distinct from an explicit zero. Genus/sample identifiers must be unique, and genus and alpha inputs must contain the same sample set. All matching is by name. Additional metadata rows are recorded as unused. Follow-up specimens are included when supplied; WTS-specific exclusions are not applied.

`Relative.mode="provided"` uses imported percentages without renormalizing a selected subset of genera. `Relative.mode="sum_genus_counts"` derives percentages using the supplied genus read total as denominator. Document the denominator, classification tool/database and taxonomic rank with the inputs.

`Alpha.mode="import"` uses supplied diversity indices. Record their calculation tool, taxonomic rank, log base and Simpson definition. `Alpha.mode="calculate"` computes genus-level Shannon as `-sum(p*log(p))` and Simpson as `1-sum(p^2)`. Diversity for an empty sample is undefined and stops execution.

The default predictor is `relative_percent`. `Predictor.metric="reads"` selects raw counts. These options use different units; input preparation and reported axes/coefficients should use the same metric.

Separate ordinary linear regressions are fitted for each genus against Shannon and Simpson. The multivariable model is:

```r
lm(Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas, data=df2)
```

Outputs include coefficients, standard errors, t-statistics, raw p-values, 95% confidence intervals, sample numbers, design matrices, matched input data and fitted model objects. These models operate at the specimen level. Missing observations, constant predictors, singular designs and insufficient residual degrees of freedom stop execution.

## Scope

The scripts operate on processed annotation, expression and microbial tables. Alignment, variant calling, fusion calling and microbial classification are performed upstream. Repeating the downstream calculations requires the corresponding processed inputs, sample mapping, analysis settings, annotation references and software versions. These are recorded in the run outputs and associated input documentation.

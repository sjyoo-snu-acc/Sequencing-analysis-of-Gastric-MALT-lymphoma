# Inputs and analysis methods

This document describes the three R files in this revision. It replaces documentation for a different configurable implementation. The source functions, expression/count processing order and four main WTS comparisons are retained where stated below; changes affecting selection or presentation are identified explicitly.

## WES

The script reads an ANNOVAR multi-annotation table and a VEP result table joined by `ID`. Multiple VEP transcripts per variant are expected. Annotation headers must otherwise be unique; a second duplicate ANNOVAR header receives `_2`, preserving the original first-block choice. Input columns used for the final table include chromosome/start/end/ref/alt, refGene/ensGene function and amino-acid annotations, `CLNSIG_220425`, `cosmic96_coding`, VEP `CANONICAL`, `BIOTYPE`, `Protein_position`, `Amino_acids`, VCF `FILTER`, `INFO`, `FORMAT`, and metadata-listed genotype columns. The example lists the minimum column layout.

Genotype columns are identified from `sample_info.tsv`, rather than fixed column numbers. FORMAT fields are parsed by their tag names. Optional phasing fields may be absent and required fields may be reordered. Records without allele 1 are excluded before parsing AD; the original workflow also required allele 1. AD must resolve to reference/alternate counts after removal of the source's unused `.` entries. SB has four counts. True multiallelic encodings require prior normalization instead of silent reinterpretation.

| Condition | Rule retained from source |
| --- | --- |
| Basic annotation | Exclude the listed intronic/intergenic/upstream/downstream/UTR/noncoding categories |
| Genotype | Allele 1 present |
| Depth | DP >= 50 and alternate AD >= 5 |
| Strand ratio | `(ref_fw + alt_fw)/(ref_rv + alt_rv)` between 0.25 and 4, inclusive |
| Transcript | `CANONICAL == "YES"`, `BIOTYPE == "protein_coding"` |
| COSMIC | `haematopoietic_and_lymphoid_tissue` or `stomach` in `cosmic96_coding` |
| Consequence | Nonsynonymous SNV, frameshift insertion/deletion, stopgain, startloss or stoploss |
| Caller filter | No `normal_artifact` label; not an added requirement for all records to be PASS |
| Population frequency | Numeric fraction < 0.1; missing annotation retained with `Keep.missing.AF=TRUE` |
| Tumor VAF | Alternate AD / total AD >= 0.05 and < 0.40 |
| Deduplication | First qualifying sample/variant record |

`AF.column="AF_eas"` is the original first ANNOVAR frequency block. The source table contains another `AF_eas` block, exposed as `AF_eas_2`. Their specific database build identities were not established from the column labels alone. The cutoff 0.1 is a fraction, not 0.1%. The numeric correction handles scientific notation and rejects malformed/range-invalid values while preserving the original missing-annotation behavior. It does not tighten the cutoff or substitute the second block.

`include_wes` determines which samples contribute called variants and per-case counts. `include_plot` determines the plotting cohort; it defaults to `include_wes` if absent. Mutation-negative cases receive zero counts. The included cases must be unique. Metadata can include a separate fusion-only plotting case. Its identity and reason for inclusion must come from the author, not from the script.

The MAF contains only called SNVs/indels. The source allowed frameshifts at filtering but omitted them from the non-synonymous plot classes; this revision retains those allowed frameshift classes for plotting. `declared_fusions.tsv` is a separate external annotation table, with no made-up coordinates, AD, DP or VAF. A separate declared-event plot combines gene events for presentation. It does not run fusion detection or the three-caller intersection described in the manuscript. `BeforeFiltering.tsv` is retained as the original intermediate filename but actually holds the selected, post-filter ANNOVAR-to-MAF input.

## WTS

The COUNT and TPM input files have identical gene/sample order: `gene_id`, one empty annotation column, then samples. Expected counts are nonnegative numeric values. With the retained `Count.mode="truncate"`, expected counts are truncated toward zero by `as.integer` before the per-comparison median-count filter (median >= 10). This is the original method; it is not rounded count import or an inferred tximport/length-offset method. `Count.mode="integer"` instead rejects fractional values. TPM is used for PCA and heatmaps, not as the DESeq2 response.

Sample names are transformed using the original five-part convention, with uniqueness and metadata agreement checked. Metadata columns are `Name`, `sample` and `HP`. `ID` is the two-character source identifier in positions 4–5 of `Name`. The source excludes ID 10 and all follow-ups from DEG analysis. Both exclusion flags are recorded for all input samples. The clinical/QC reason for the ID exclusion was not supplied. Source identifiers must not be assumed to equal manuscript case aliases.

| Comparison | Numerator | Denominator | Design |
| --- | --- | --- | --- |
| S1 / S5 / S9 | PNormal | NNormal | `~ Type` |
| S2 / S6 / S10 | PTumor | PNormal | `~ ID + Type` |
| S3 / S7 / S11 | PTumor | NTumor | `~ Type` |
| S4 / S8 / S12 | NTumor | NNormal | `~ Type` |

Each model is fitted separately to its selected samples. Matrices are ordered by exact metadata names. The patient-adjusted design disallows duplicated ID/type observations and must have full rank; the original sample inclusion is retained even if some patients lack one tissue. Other comparisons remain unpaired, including the HP-negative tumor/normal comparison. Those choices are not silently changed to improve agreement with a reported result.

`ENSG2SYM` defaults to the supplied local `gene_symbols.tsv` snapshot. Mapping is matched to input ENSG IDs, with optional version suffixes removed for matching. The original first mapping per ENSG and first input gene per duplicated symbol are retained; empty/unmapped symbols are excluded. The chosen mapping and a two-column reusable snapshot are saved. Live Ensembl lookup requires `Query.Ensembl=TRUE`; record `Ensembl.version` and the actual snapshot origin. A live/current snapshot must not be described as the unknown historical study annotation.

DESeq2 uses `results(..., contrast=c("Type", numerator, denominator), pAdjustMethod="fdr")`, then orders by padj. This revision defines up/down consistently as LFC >= 1 / <= -1 with padj <= 0.05. The source instead used raw p-values for those labels while using padj on the volcano axis and threshold. Consequently, selected lists, counts and point colors can change. The raw p-value column is preserved. No p-value/padj columns are exchanged to match a workbook. Zero padj is capped only for the graphical logarithm, not in result tables.

`degsea` retains its source arguments and result component names. `pheatcustom` retains its three arguments and takes the fitted comparison's metadata rather than a later global variable. The heatmap uses the original top genes in padj order, log2(TPM+1), row scaling and ordered sample annotations. Rows with undefined zero-variance scaling are omitted from the plot only. This does not filter the DESeq2 results.

### Enrichment and reference replay

With `Run.GSEA=TRUE`, GO BP and GO ALL use the source's positive ranks of LFC, random tie ranks, `minGSSize=3`, BH adjustment, `scoreType="pos"` and `eps=0`. This positive-score analysis must not be interpreted as a separate test for negative enrichment. KEGG uses the source's signed LFC ranking after SYMBOL-to-ENTREZID conversion, first mapping per symbol, `minGSSize=3`, BH adjustment and `eps=1e-10`. The source seed 48575 is reset before each comparison's GSEA, and the backend is serial. The reset is a reproducibility correction: the source set the seed only once, which allowed optional plotting or preceding comparisons to affect later Monte Carlo draws. Seed states before GO/KEGG and ranking vectors are saved. Numerical reproducibility also requires the same package versions, references, data/order and execution order.

The default `Use.reference.cache=FALSE` uses the installed GO/OrgDb annotations and current KEGG service. Result objects and snapshots of the tested gene sets are saved under the run's `GSEA_references/`, separately for each comparison. These snapshots describe the new run, not the original study databases. Their tested-universe restriction is recorded.

For a repeat of that calculation, set `Reference.dir` to the previous run's `GSEA_references` and `Use.reference.cache=TRUE`, restore the recorded packages, and keep the tested gene universe unchanged. Per-comparison caches are named `GO_BP_<numerator>_<denominator>.rds`, `GO_ALL_<numerator>_<denominator>.rds` and `KEGG_<numerator>_<denominator>.rds` when `Explanation` is empty. A supplied whole-reference cache named `GO_BP.rds`, `GO_ALL.rds` or `KEGG.rds` is also accepted. Cache objects contain `term2gene`, `term2name`, `package_versions` and `eps`; captured tested-universe caches also contain `input_genes`. Package versions are checked, and input reference files are hashed. Exact statistical replay of this new cache path was not tested in the available environment; see VALIDATION.md.

`Run.Hallmark=FALSE` makes the source's additional Hallmark analysis optional. Enabling it requires a documented local `hallmark_term2gene.tsv` with `gs_name` and `human_gene_symbol`. It reuses the source's positive-rank GSEA settings. No automatically updated MSigDB download is substituted for an unknown historical release.

## Microbiome

The active analysis is extracted from the genus composition, alpha-diversity and regression sections of `250326.R`, which were previously disabled inside WTS.R. It runs independently and does not depend on undefined `Group2`/`Group3` objects or the previous WTS global environment.

| Input | Columns and meaning |
| --- | --- |
| `genus_counts.tsv` | First column `genus`; one nonnegative integer read-count column per sample |
| `genus_relative_percent.tsv` | Same genera/samples; supplied percentages in [0,100], with the actual denominator documented |
| `alpha_diversity.tsv` | `sample`, `Shannon`, `Simpson_index`; original externally calculated indices by default |
| `grouptable_new.tsv` | `Name`, `sample` (Normal/Tumor/Eradicated), `HP` (positive/negative) |

Use genus-level rows, not a species-only H. pylori summary. The four source predictors are Helicobacter, Rothia, Veillonella and Porphyromonas. Sample matching is always by name, including alpha indices; column/row positions are not assumed. Additional metadata samples are recorded as unused. Genus and alpha inputs must describe the same sample set. Follow-ups remain included when supplied, without carrying over the WTS-specific exclusions.

`Relative.mode="provided"` preserves imported percentages; the genus composition plot displays those percentages without renormalizing a selected subset to 100%. Optional `sum_genus_counts` uses the supplied complete genus read total as denominator and is a deliberate alternative. The script rejects negative/nonfinite counts, impossible percentages and duplicate taxon/sample IDs. Missing/unknown observations are not replaced with a zero.

`Alpha.mode="import"` preserves the supplied Shannon and Simpson indices; document their original tool, rank, log base and Simpson definition. No alpha tool version or calculation rank is inferred from a filename. Optional `calculate` computes genus-level Shannon as `-sum(p*log(p))` and Simpson as `1-sum(p^2)` from the count matrix. Empty-sample diversity is undefined and stops rather than being reported as zero. These optional recomputations are not claimed to reproduce the historical analysis.

The default predictor is `relative_percent`, because the source regression used `reads_pct`. The option `Predictor.metric="reads"` is an explicit reanalysis using raw counts, not an automatic substitute for the source. This choice matters for coefficients and must agree with the manuscript axis labels and stated units.

For each of Shannon and Simpson, a separate `lm` is fitted to each of the four genus predictors. The original multivariable model is retained exactly:

```r
lm(Shannon ~ Helicobacter + Rothia + Veillonella + Porphyromonas, data=df2)
```

Coefficients, standard errors, t-statistics, raw p-values, 95% confidence intervals, n, design matrix, keyed data and R model objects are exported. Missing observations, a constant predictor, singular multivariable design or insufficient residual degrees of freedom stop explicitly. No coefficient sign, p-value or manuscript value is hardcoded. These are the source's specimen-level ordinary linear models; repeated-patient adjustment is not added silently. The HP group-comparison test behind the manuscript's p=0.4 was not established from the supplied source and is not invented here.

## Scope and remaining reconciliation

These scripts start with processed tables. They do not run raw-read alignment/calling, the three fusion callers, Kraken2/Bracken, or the historical upstream alpha-diversity workflow. The specimen-to-case map, exclusion reasons, historical package/annotation environment, GO/KEGG/MSigDB references, microbial denominator and alpha-index definition still require author-supplied provenance where absent.

The original source, final workbooks and manuscript may disagree. The revision does not alter workbook values, reverse a contrast, exchange p-value columns, change a regression sign or optimize model/filter settings to force agreement. In particular, the manuscript's positive-normal expression directions, its Shannon coefficient sign/interpretation, source exclusion/cohort counts, and existing WTS table discrepancies require reconciliation against the original analysis outputs. The validation record distinguishes successful computation from exact manuscript-result reproduction.

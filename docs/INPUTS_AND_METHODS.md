# Inputs, methods and outputs

## Common contract

Open WES.R, WTS.R or Microbiome.R, set Project.dir and the Parameters section, then click Source once. The chosen file contains everything required for its analysis; no other R script is loaded. Missing required packages are installed by that same file when Install.packages is TRUE. Missing study choices and incorrect inputs stop with explicit errors. Run.ID and Output.folder define a new output directory. Successful execution means that the configured workflow completed; manuscript_result_reproduction remains NOT_ASSESSED. Normal runs are labelled analysis in the execution record; the test copies declare synthetic_test explicitly.

Every input is tab-delimited with a single header row, unique column names and no quoted embedded tabs/newlines. Identifiers are read as strings. Missing numbers may be `NA`, empty or `.` only where explicitly supported. Counts may never be negative, infinite or silently rounded. Non-count annotations must be explicitly listed. No metadata is inferred from column position, a specimen-name substring or a numerical case suffix.

SET contains the inputs, seed and method choices; WTS comparisons use VAR and enrichment uses GO/KEGG. Execution purpose/status are recorded automatically, with exploratory as the default unconfirmed study status. Production defaults do not read synthetic settings from another file. The resolved R settings are written to configuration_used.json; no YAML file is needed. An actually confirmed author_confirmed declaration must not be substituted for missing evidence.

Run artifacts: `configuration_used.json`, `input_hashes.tsv`, `run_status.json`, `sessionInfo.txt`. Reference snapshots used as files are hashed alongside analysis inputs. Original environment versions are unknown. Package versions in run logs are measurements of the new execution environment. Live annotation labels are declarations of provenance, not a frozen historical database.

## WES

### Input tables

| File | Required columns/content |
| --- | --- |
| ANNOVAR | `ID`, `Chr`, `Start`, `End`, `Ref`, `Alt`, `Func.refGene`, `ExonicFunc.refGene`, `ExonicFunc.ensGene`, `Gene.refGene`, `cosmic96_coding`, `FILTER`, `INFO`, `FORMAT`, and explicitly named sample genotype columns |
| VEP | `ID`, `CANONICAL`, `BIOTYPE`; `Protein_position` and `Amino_acids` for S13 comparison; `AF_eas` for the second profile if not supplied in ANNOVAR |
| MAF conversion annotations | `GeneDetail.refGene`, `AAChange.refGene` in addition to standard ANNOVAR annotation fields |
| Sample metadata | `sample_column`, `sample_id`, `role` (`tumor`/`normal`), `include` (`TRUE`/`FALSE`); every included genotype column must exist |
| Optional verified case map | `sample_id`, `case_number`; unique, nonmissing one-to-one mapping for all included tumor samples, including mutation-negative samples |

Inputs must use the same normalized biallelic IDs/coordinates. ANNOVAR IDs must be unique. Duplicate original headers require an explicit preprocessing decision, not automatic positional renaming. Overlapping ANNOVAR/VEP columns other than ID cause an error; resolve by an authorized, documented rename upstream. VEP one-to-many annotation rows require `vep_duplicate_policy: retain_all` or an upstream confirmed transcript selection; the default choice is unresolved in the study template. Missing VEP annotations fail relevant filter rules. No liftover, re-calling, normal subtraction or genotype imputation is performed.

FORMAT is parsed by tag name for each record, supporting reordered GT/AD/DP/SB tags. Optional AF is preserved separately; VAF uses AD_alt/(AD_ref+AD_alt). SB ratio uses all forward reads / all reverse reads, matching the source. It is not alternate-only strand bias. SB must be the per-sample four-count FORMAT field. INFO DP/AS_SB_TABLE is retained as annotation and never used to replace FORMAT support.

No-call GT is omitted even when it has extra placeholder fields. Called multiallelic GT/AD/ALT causes an error; normalize upstream with a documented method. A trailing placeholder-only AD allele is explicitly removed. Missing metrics/zero VAF denominator fail the corresponding rule. A zero reverse denominator yields infinite ratio when forward support exists, otherwise NA; neither passes a bounded strand-ratio rule.

### Separate observed filter profiles

| Rule | wes_r_first | wes_r_second | wes_rmd_strict |
| --- | --- | --- | --- |
| Source | WES.R 131-162 | WES.R 222-243 | WES.Rmd 111-128 |
| AD_ref | No positive minimum | No positive minimum | >=20 |
| AD_alt | >=3 | >=3 | >=20 |
| Alternate forward/reverse | Each >=1 | No positive minimum | Each >=10 |
| FORMAT DP | >=20 | >=50 | >=40 |
| All-read forward/reverse ratio | [0.125,8] | [0.25,4] | [0.25,4] |
| VAF | [0.1,0.4] OR [0.6,0.9] | [0.1,0.3) | (0.05,0.45) |
| Canonical | YES | YES | YES |
| BIOTYPE | protein_coding | Not required | protein_coding |
| Consequence source | RefGene | ensGene | RefGene |
| Consequences | nonsynonymous SNV, frameshift insertion/deletion, stopgain, startloss | Same plus stoploss | Same as first |
| COSMIC96 | Contains hematopoietic/lymphoid tissue or stomach | Annotation present | Contains hematopoietic/lymphoid tissue or stomach |
| normal_artifact | Excluded | Included exclusively | No FILTER-based exclusion |
| AF_eas | Not used | <0.01 or explicit dot | Not used |

All profiles retain the source's common functional-category exclusions and require alternate allele 1. Exact GT allele tokens replace the original loose string match. Only explicitly labeled tumor samples produce candidate variant results; upstream paired somatic calling provenance remains the user's responsibility. The second historical profile requires a `normal_artifact_rationale`; this is a configuration record, not confirmation that the inclusion is appropriate.

Transcript duplicates are handled by explicit `error` or `first_source_row`. `first_source_row` is a source-order choice after the selected filters and never generates the absent Cartesian combinations produced by the original table()/match() approach. Additional Ensembl consequence exclusions in WES.Rmd 169-171 are a separate `post_ens_exclusions` switch. Original sample exclusions are not imposed automatically. The final profile, exclusion list and case map remain unresolved.

Outputs: `all_genotypes_with_filter_audit.tsv` with per-record selection/reasons; `variants.tsv`; per-tumor `variant_counts.tsv` including zero counts; `filter_profile_used.json`; optional `annovar_for_maf.tsv`, `mutations.maf.tsv`, `oncoplot_snv_indel.pdf`. MAF conversion calls installed `maftools::annovarToMaf`, not the undated supplied copy. Its gene/transcript/classification behavior depends on the measured package version. Cohort size is the number of included tumor samples, not a hard-coded number or the number with variants. This oncoplot is a new SNV/indel visualization, not a reproduction of an unavailable manuscript figure. No WTS fusion calls are manufactured or added to WES MAF.

## WTS

### Inputs and sample/model validation

Counts: `gene_id` followed by sample columns, plus any explicitly declared `count_annotation_columns`. They must be raw nonnegative integers. TPM and fractional expected-count inputs are rejected. The original source coerces to integer; the rewrite stops instead of silently truncating. If expected counts/tximport are needed, confirm that method separately; it is not implemented as an inferred substitute.

Metadata requires `sample_id`, `patient_id`, `type` (`Normal`, `Tumor`, `After`), `hp_status` (`Positive`, `Negative`, `Unknown`), `group`, `include` (`TRUE`/`FALSE`) and all configured design variables. Every count sample, including excluded samples, must have exactly one metadata row. Samples are aligned by ID. Standard `HPpos_`/`HPneg_` group labels must agree with explicit HP/type columns; labels are never generated from specimen IDs. Other explicitly defined group names are supported.

Each comparison requires `id`, contrast `factor`, `numerator`, `denominator`, additive `design`, `model_kind`, `repeated_patients`, and a `reference_levels` entry for every categorical design variable. Numeric covariates must be explicitly listed. Supported formulas are additive fixed effects, such as `~ group` or `~ patient_id + group`; interactions/random effects need a separately reviewed implementation.

`paired` requires the patient term and exactly one sample for each patient in both contrast groups. Incomplete pairs/repeated biopsies are rejected rather than silently discarded or collapsed. `independent` excludes patient blocking and rejects repeated patient IDs unless `repeated_patients: explicit_unpaired` and a recorded `unpaired_rationale` are provided. That option remains an explicitly different analysis. A confounded or rank-deficient design and a model without residual degrees of freedom stop. In particular, between-patient HP comparisons may be confounded with a patient fixed effect; do not repair rank by dropping terms automatically.

`fit_scope: per_comparison` fits each selected contrast population separately. `all_included` fits once to all included samples and requires a consistent model/reference configuration across comparisons. These scopes can change size factors, dispersion estimates, baseMean and adjusted pvalues; neither is selected as the manuscript default.

The original active notebook removes all-zero features and uses `~ ID + Type` with contrasts in the opposite order to some Tumor/Normal table captions. Generic DEG scripts permit externally defined group contrasts and optional additive metadata. Four table captions alone cannot establish which source fit, pairing, sample subset or direction produced the reported numbers.

### Annotation and DEG

`gene_id_policy` explicitly preserves IDs or strips a terminal version suffix. Duplicate features after transformation must error or be explicitly summed at the count level. Different Ensembl features sharing a symbol are not silently collapsed. Optional mapping TSV requires `gene_id`, `symbol`, `entrez_id`; statistics remain indexed by original feature ID. Duplicate mapping rows error or use an explicit first-source-row policy. Missing symbols retain feature IDs.

WTS reads a supplied local mapping snapshot. Keyed mapping and explicit duplicate/label policies retain the needed Toolbox concepts; unused live Get_biomaRt/ENSG2SYM query/export code is not part of the runtime. No new mapping or reference version is invented. The source's source_legacy_integrated issue remains evidence: Toolbox.R:723 replaces nonduplicated integrated labels with Ensembl IDs. Confirm the intended snapshot/label policy before providing the mapping file; do not choose labels to match the final GENE column.

DESeq configuration explicitly records `test` (Wald), `fitType`, `sfType`, `minReplicatesForReplace`, results `alpha`, `pAdjustMethod`, `independentFiltering`, `cooksCutoff`. Parameters omitted by the original notebook are not declared to be historical values. `deg.padj_cutoff` and `deg.abs_log2fc_cutoff` define a separate selected table; full statistics are retained, including NA. No LFC shrinkage or pvalue replacement is added. The observed source selection is |LFC|>=1 and padj<=0.05; this does not prove the final table-selection rule.

Each comparison writes DEG_all.tsv, DEG_selected.tsv, DEG_upregulated.tsv, DEG_downregulated.tsv, metadata_used.tsv, design_matrix.tsv, contrast_used.json and deseq_dataset.rds. The separate UP/DOWN summaries add gene_type; the complete/selected files retain their previous columns and numerical values. Main Table 2 contains selected genes without a deterministic source selection rule; the code does not invent a top-gene list to match it.

### GO/KEGG GSEA

`gsea.enabled` must be explicit. Rankings use configured log2FoldChange (source) or Wald stat (new exploratory choice), finite values, descending order with stable source-order ties. `all_tested` reflects the WTS source. `pvalue_lt` is the separate generic-script prefilter and requires a cutoff. Do not silently treat them as equivalent. KEGG mapping is matched by feature ID; unmapped and duplicate rank IDs are counted, with explicit duplicate policies. A positional mapping assignment in the generic script is not reused.

For each GO/KEGG branch, choose `source_api` or `local_term2gene`, give a reference label, and explicitly configure min/max gene-set size, pvalue cutoff, adjustment method, exponent, eps and backend seed policy. Version-specific settings such as `by`, `method`, `nPerm` must match the installed clusterProfiler API; unknown/incompatible settings fail normally, without a fallback that changes methods. The rewrite calls the real package statistical engine, not an invented GSEA implementation.

`source_api` GO additionally requires OrgDb, ontology and key type. KEGG requires organism, key type and internal/live-data choice. Live KEGG is not a frozen historical snapshot. Local TERM2GENE files require `term_id`, `gene_id` and TERM2NAME requires `term_id`, `description`; duplicates/missing descriptions error. GO ranking IDs must match supplied GO gene-set IDs and KEGG ranking IDs the configured Entrez mapping. Actual local reference files are hashed. Local snapshots can avoid live database changes, but equivalence to the source APIs must be independently established.

WTS.Rmd paired GO uses ALL/ENSEMBL, min 3, max 800, BH. Paired KEGG uses min 3, max 800, nPerm 10000 and `none`. Other HP-read blocks use BH, while generic host scripts use a pvalue<0.05 ranking subset, GO 10-500 and KEGG minimum 8 with adjustment defaults. The final S9-S12 values have adjusted columns different from raw pvalues. No automatic correction from none to BH is made. Database/API versions and complete tested universes remain unavailable.

Outputs: `go_ranking.tsv`, `kegg_ranking.tsv`, `GO_GSEA.tsv`, `KEGG_GSEA.tsv` and reference/ranking summaries. Artificial local GSEA fixture IDs start with `SYN_` and do not represent biological GO/KEGG pathways.

## Microbiome

The module summarizes an existing taxonomic table; it does not run a classifier. Original Microbiome.Rmd imports Pavian metrics and does not define their denominator or upstream Z-score calculation. The unused WTS/TPM, live mapping and MSigDB initialization in that notebook is excluded.

| Input | Required contract |
| --- | --- |
| Wide taxa table | Explicitly configured ID/name/rank column names, species rows and all named metric columns |
| Column map | `input_column`, `sample_id`, `metric`; metric is count, relative_abundance, zscore_count or zscore_relative; every sample/metric pair unique |
| Metadata | `sample_id`, `patient_id`, `group` (Normal/Tumor/After), `include`; all mapped samples represented |
| Optional denominator | `sample_id`, `total_reads`; nonnegative integer denominator for every included sample |
| Optional reference list | `sample_id`; an explicit nonempty subset of included samples, not inferred controls |

Select an exact species name or verified taxon ID at an explicit rank. If both are supplied, they must agree. The target may occur anywhere in the table. Counts, imported fractions and upstream Z-scores remain distinct. Missing taxa error unless an explicit confirmed-zero policy and nondetection basis are supplied; absence is not assumed to be nondetection.

`duplicate_taxa: sum_counts_only` allows an explicit count-only aggregation, propagating NA. Imported fractions/Z-scores cannot be summed. Negative/fractional read counts and out-of-range imported abundances error. NA remains distinct from measured zero. The denominator must be explicitly described; no sum across mixed taxonomic ranks or hierarchical clades is used.

Relative mode preserves imports (percent converted to a separately recorded fraction scale), computes count/explicit total reads, or disables abundance. A zero denominator errors or yields NA by explicit policy. Count-derived abundance never replaces the original imported file. Upstream Z-scores are kept with their stated reference provenance, which may explicitly remain unknown. Imported zscore_relative is not rescaled with percent-to-fraction abundance.

Optional recomputed Z-score requires selected metric, identity/log2(x+1) transformation, explicit reference samples, reference-NA policy and zero-variance policy. It uses reference mean and sample SD (n-1), needs >=2 measured reference values, and writes `zscore_recomputed` as a distinct column. It is a new explicit derivation unless the authors confirm it; it is not inferred to be Pavian's original Z-score.

Optional median grouping uses all included measured abundances or an explicit reference list. The source's `>= median` rule is retained. `median_na_policy: legacy_down` reproduces its missing-to-DOWNmedian behavior; `preserve` keeps missing classification NA. This descriptive split is not a clinical H. pylori positivity cutoff.

Outputs: target_imported_metrics.tsv, H_pylori_metrics.tsv, quantitation_provenance.json. The R objects H.pylori_reads_pct (percent display) and H.pylori_Zscore_pct (upstream Z-score) retain the original summary names; the main saved abundance column remains a fraction. MTD VST/batch-adjusted values are not fed back as counts. No final microbiome quantitation sheet was provided; clinical HP status in Table 1 does not substitute for such a sheet.

## Read-only comparison with supplied final tables

The following is an optional verification tool, separate from routine module execution. Create a private JSON manifest from tools/comparison_manifest.template.json and keep the filled-in copy under data/. Paths refer to already executed outputs. The utility requires explicit matching IDs and numeric tolerances; it never changes a result to make it match. Every comparison directory is new and should be kept outside the public repository or under ignored results/.

```sh
python3 tools/compare_to_excel.py --supplement data/private_supplement.xlsx \
  --manifest data/local_comparison.json --output results/private-comparison-001
```

WES manifest entry:

```json
{"variants":"results/WES/RUN/variants.tsv","columns":{"case_number":"case_number","gene_symbol":"Gene.refGene","variant_function":"ExonicFunc.refGene","protein_position":"Protein_position","amino_acids":"Amino_acids","cosmic96_coding":"cosmic96_coding","INFO":"INFO","Chr":"Chr","Start":"Start","End":"End","Ref":"Ref","Alt":"Alt"}}
```

Put that object under `wes`. Case alias + Chr/Start/End/Ref/Alt is the exact record key. Annotation text differences are reported. Per-sample DP/VAF/strand counts cannot be numerically validated from S13 because it does not report them; INFO fields are compared as reported annotation, not substituted metrics.

WTS manifest entry (place in the `wts` array):

```json
{"table_number":1,"results":"results/WTS/RUN/S1_S5_S9/DEG_selected.tsv","expected_numerator":"AUTHOR_VERIFIED_GROUP_A","expected_denominator":"AUTHOR_VERIFIED_GROUP_B","absolute_tolerance":0,"relative_tolerance":1e-8}
```

These group placeholders must be replaced with verified configured labels. The comparator reads contrast metadata and verifies type/HP semantics against each final caption before numeric comparison. S1-S4 contain symbols rather than Ensembl IDs; only a verified unique symbol mapping is usable, and duplicates error. GO/KEGG comparisons use exact pathway IDs. All full/selected statistics and unmatched records remain visible in the comparison JSON. Choose tolerances appropriate to the reported precision, especially very small pvalues; no rounding or value exchange is performed.

A microbiome manifest can use separately authorized reported/computed TSV paths, an explicit metric list and tolerances. It requires taxon ID + sample ID and taxon name. The supplied package has no such final table, so that comparison remains NOT_RUN. Providing final Excel alone does not constitute raw study inputs or authorize reproducing analysis from the final results.

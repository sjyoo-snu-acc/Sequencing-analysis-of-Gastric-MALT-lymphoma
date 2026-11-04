# Methods and source notes

Input requirements and analysis provenance for the three standalone scripts.

## Source and final-result mapping

All original-source names and line ranges below refer to the untouched handoff originals in 02_original_code, not the newly reconstructed standalone files with similar names.

This audit was performed on the supplied final workbooks and original code. The relationship is an evidence-based implementation link, not proof that a particular script produced the submitted numbers. All source files remain unchanged outside this release candidate.

| Final workbook/sheet | Purpose | Source range | Required input | Decision and evidence |
| --- | --- | --- | --- | --- |
| Table_1-2.xlsx / Table 1 | Clinicopathologic characteristics | None | Clinical records (not provided) | Excluded from omics modules. Age/gender/location/HP columns are descriptive clinical data |
| Table_1-2.xlsx / Table 2 | Selected DEG summaries in four HP comparisons | WTS sources listed for S1-S4 | Same WTS inputs plus author-defined gene selection | Summary relationships verified; computational reproduction not run. All 25 entries agree with S1-S4 within 1e-8 LFC and 1e-15 padj absolute tolerances |
| Additional File 2. supplemental Table.xlsx / Index | Titles of S1-S13 | None | None | Documentation only. All 13 index entries were read and matched to captions |
| Additional File 2. supplemental Table.xlsx / Table S1 | DEG; H. pylori-positive normal group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S2 | DEG; H. pylori-positive cancer group vs H. pylori-positive normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S3 | DEG; H. pylori-positive cancer group vs H. pylori-negative cancer group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S4 | DEG; H. pylori-negative cancer group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S5 | GO GSEA; H. pylori-positive normal group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S6 | GO GSEA; H. pylori-positive cancer group vs H. pylori-positive normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S7 | GO GSEA; H. pylori-positive cancer group vs H. pylori-negative cancer group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S8 | GO GSEA; H. pylori-negative cancer group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S9 | KEGG GSEA; H. pylori-positive normal group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S10 | KEGG GSEA; H. pylori-positive cancer group vs H. pylori-positive normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S11 | KEGG GSEA; H. pylori-positive cancer group vs H. pylori-negative cancer group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S12 | KEGG GSEA; H. pylori-negative cancer group vs H. pylori-negative normal group | WTS.Rmd:305-443; DEG_Anno_Plot.R:64-127,190-214,1356-1433; DEG_Anno_Plot2.R analogous blocks | Raw integer counts, verified metadata/pairing, explicit model/contrast/fit population, annotation and reference DB versions | Implementation selected; generation of these exact four comparisons unverified. Workbook titles specify four HP strata; active WTS notebook instead compares Normal/Tumor/After; generic scripts use an absent sample sheet |
| Additional File 2. supplemental Table.xlsx / Table S13 | WES somatic variant metrics | WES.R:9-90,131-190,222-265; WES.Rmd:11-17,95-150; annovarToMaf.R:31-220 | ANNOVAR, VEP, named genotype columns, verified sample-role and case-alias map, selected filter chain, reference versions | Related source selected; exact filter chain and row provenance unresolved. Field-level overlap; S13 lacks per-sample FORMAT metrics and alias mapping |
| No supplied microbiome result sheet | H. pylori reads, relative quantity, upstream Z-score | Microbiome.Rmd:51-94; supporting WTS.Rmd:698-744 | Pavian/Kraken/Bracken output, taxonomy, explicit metric column map, metadata, denominator and upstream Z-score provenance | Exploratory source selected; final-result numerical comparison unavailable. Neither supplied workbook contains microbiome quantitation; clinical HP status is not a microbiome output |

## Workbook coverage

| Sheet | Data rows | Exact column headers |
| --- | ---: | --- |
| Table 1 | 13 | No; Age; Gender; Location; H. pylori |
| Table 2 | 25 gene entries, four comparison labels | GENE; log2FoldChange; padj (fourth used-range column empty) |
| Index | 13 titles | No header |
| Table S1 | 39 | GENE; baseMean; log2FoldChange; pvalue; padj |
| Table S2 | 391 | GENE; baseMean; log2FoldChange; pvalue; padj |
| Table S3 | 34 | GENE; baseMean; log2FoldChange; pvalue; padj |
| Table S4 | 766 | GENE; baseMean; log2FoldChange; pvalue; padj |
| Table S5 | 609 | ID; Description; NES; pvalue; p.adjust |
| Table S6 | 1298 | ID; Description; NES; pvalue; p.adjust |
| Table S7 | 1102 | ID; Description; NES; pvalue; p.adjust |
| Table S8 | 766 | ID; Description; NES; pvalue; p.adjust |
| Table S9 | 80 | ID; Description2; NES; pvalue; p.adjust |
| Table S10 | 74 | ID; Description2; NES; pvalue; p.adjust |
| Table S11 | 84 | ID; Description2; NES; pvalue; p.adjust |
| Table S12 | 57 | ID; Description2; NES; pvalue; p.adjust |
| Table S13 | 22 | case_number; gene_symbol; variant_function; protein_position; amino_acids; cosmic96_coding; INFO; Chr; Start; End; Ref; Alt |

Each supplementary sheet has a final caption row, excluded from data-row counts. All sheets are visible. Cell values, merged ranges, hidden rows/columns, formulas and comments were inventoried privately. No source workbook was edited.

## Conflicts requiring source verification

- S4: every data row has raw pvalue greater than padj; all raw pvalues exceed 0.05 while all padj values are below 0.05. This is incompatible with interpreting the pair as ordinary BH raw/adjusted pvalues from the same test. No columns or values were swapped. Full test-universe data are absent.
- GO in WTS.Rmd uses BH at lines 384-393. Its paired KEGG block explicitly uses none at lines 420-427. Every reported KEGG row in S9-S12 has p.adjust different from pvalue. The final method cannot be resolved from these tables alone.
- WTS.Rmd has active ~ ID + Type comparisons in the direction Normal/Tumor, Tumor/After, After/Normal. Its four HP-stratified table contrasts are not encoded. Generic DEG scripts permit external group contrasts but no actual sample sheet is supplied.
- The three WES chains disagree on depth, VAF windows, strand support, annotation source and normal_artifact inclusion. None is selected as the final manuscript chain by default.
- S13 lacks verified links between case aliases and original genotype columns. Site INFO AS_SB_TABLE/DP cannot substitute for per-sample FORMAT SB/DP. SNV/indel results are distinct from WTS fusion data; no fusion input or detection pipeline was supplied.

## Source selection and exclusions

| Original file | Decision | Reason |
| --- | --- | --- |
| WES.R | Method reconstruction: first and second filter blocks remain separate explicit profiles | Overlapping S13 annotation fields; inconsistent filters prevent a final default |
| WES.Rmd | Method reconstruction: separate strict profile | Same input/annotation workflow, higher support thresholds; later Ensembl exclusion remains optional |
| annovarToMaf.R | Inspected; not redistributed | Strong similarity to upstream maftools; supplied file has no header establishing revision/license provenance. Call installed maftools through a validated adapter instead |
| WTS.Rmd | DESeq2, rankings, GO/KEGG implementation selected | Active paired workflow retained as a source-supported option, without inventing patient pairing |
| Toolbox.R | Only gene-ID mapping concept retained and rewritten | Mapping functions rely on a global mart and live annotation; unrelated scoring/abstract tools excluded |
| Microbiome.Rmd | Pavian metric extraction selected and rewritten | Explicit taxon/sample/metric matching replaces row and column position assumptions; upstream Z-scores preserved |
| DEG_Anno_Plot.R and DEG_Anno_Plot2.R | Additional method evidence; full scripts excluded | Generic host models and GSEA differ from WTS.Rmd. Version 2 changes CSV name handling, an H. pylori heatmap, plotting and output directory behavior, not a verified final table model |
| FINAL.Rmd | Excluded | Amplicon IDH workflow, cohortSize 236, targeted gene/hotspot selection and manually appended FLT3 records; no support for S13 lineage |
| Early_MALToma.R | Excluded | TPM density/PCA exploration; no final table computation |
| MTD_normalization.R | Excluded | VST/batch-adjusted values rounded into taxonomic read fields. They are not original counts; no final-table contribution established |
| MTD_microbiome_profiling.R | Excluded from executable release | Broad taxonomic/host/differential pipeline with absent inputs; not a verified contributor to a final microbiome table |
| comparison_MTD_SAMSA2_plot.R | Excluded | Tool comparison plotting; no corresponding final table supplied |
| Preliminary draft and previous notes | Reference only | Preliminary unpaired model and first-row taxon assumption are not author-verified methods |

Original line numbers refer to unmodified files in the handoff. No figures/manuscript were supplied; no figure layout, fusion result or prose claim has been independently checked.

## Original analysis structure

The supplied research scripts are the basis of the calculations and thresholds. The reorganized scripts retain the original-style objects and main analysis flow. The source line ranges in this document refer to the untouched handoff, not the current files with the same names. An implementation link does not establish that an unverified branch generated the reported results.

| Module | Retained analysis structure | Later additions or corrections |
| --- | --- | --- |
| WES | ANV/VEP, GML, BEF, FIN-FIN4, tmp, MAF and before; AD-based VAF and all-read strand ratio; distinct source filter chains | Named FORMAT/sample matching, explicit roles/case map, validated joins/deduplication, installed converter and dynamic cohort size |
| WTS | CNT, MET, VAR, pre/des/med/obj/res/wrt/cut/UP/DOWN and gnl/gse/kk2; DESeq2 fit visible in the comparison loop | Count/model/metadata checks, explicit references and fit population, keyed ID mapping and run records |
| Microbiome | pavian_obj, H.pylori, H.pylori_reads_pct and H.pylori_Zscore_pct; imported metrics and optional original median rule | Explicit species/rank/sample/metric selection and separately named optional recomputations |

FIN and FIN2 retain the non-COSMIC intermediate stage; FIN3 contains the complete selected profile; FIN4 applies the stated duplicate/alias policy. UP and DOWN are disjoint directional summaries of the selected DEG rows. The saved microbiome abundance remains a fraction, while H.pylori_reads_pct displays percent and H.pylori_Zscore_pct preserves the imported score.

Input checks, package preparation, run records and safe file handling were added during reorganization. Optional count-derived abundance or recomputed Z-scores are new analyses unless the original authors confirm their equivalence. No original source or final workbook was modified.

## Inputs, methods and outputs

### Common contract

Open WES.R, WTS.R or Microbiome.R, set Project.dir and the Parameters section, then click Source once. The chosen file contains everything required for its analysis; no other R script is loaded. Missing required packages are installed by that same file when Install.packages is TRUE. Missing study choices and incorrect inputs stop with explicit errors. Run.ID and Output.folder define a new output directory. Successful execution means that the configured workflow completed; manuscript_result_reproduction remains NOT_ASSESSED. Normal runs are labelled analysis in the execution record.

Every input is tab-delimited with a single header row, unique column names and no quoted embedded tabs/newlines. Identifiers are read as strings. Missing numbers may be `NA`, empty or `.` only where explicitly supported. Counts may never be negative, infinite or silently rounded. Non-count annotations must be explicitly listed. No metadata is inferred from column position, a specimen-name substring or a numerical case suffix.

SET contains the inputs, seed and method choices; WTS comparisons use VAR and enrichment uses GO/KEGG. Execution purpose/status are recorded automatically, with exploratory as the default unconfirmed study status. Production defaults do not read synthetic settings from another file. The resolved R settings are written to configuration_used.json; no YAML file is needed. An actually confirmed author_confirmed declaration must not be substituted for missing evidence.

Run artifacts: `configuration_used.json`, `input_hashes.tsv`, `run_status.json`, `sessionInfo.txt`. Reference snapshots used as files are hashed alongside analysis inputs. Original environment versions are unknown. Package versions in run logs are measurements of the new execution environment. Live annotation labels are declarations of provenance, not a frozen historical database.

### WES

#### Input tables

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

#### Separate observed filter profiles

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

### WTS

#### Inputs and sample/model validation

Counts: `gene_id` followed by sample columns, plus any explicitly declared `count_annotation_columns`. They must be raw nonnegative integers. TPM and fractional expected-count inputs are rejected. The original source coerces to integer; the rewrite stops instead of silently truncating. If expected counts/tximport are needed, confirm that method separately; it is not implemented as an inferred substitute.

Metadata requires `sample_id`, `patient_id`, `type` (`Normal`, `Tumor`, `After`), `hp_status` (`Positive`, `Negative`, `Unknown`), `group`, `include` (`TRUE`/`FALSE`) and all configured design variables. Every count sample, including excluded samples, must have exactly one metadata row. Samples are aligned by ID. Standard `HPpos_`/`HPneg_` group labels must agree with explicit HP/type columns; labels are never generated from specimen IDs. Other explicitly defined group names are supported.

Each comparison requires `id`, contrast `factor`, `numerator`, `denominator`, additive `design`, `model_kind`, `repeated_patients`, and a `reference_levels` entry for every categorical design variable. Numeric covariates must be explicitly listed. Supported formulas are additive fixed effects, such as `~ group` or `~ patient_id + group`; interactions/random effects need a separately reviewed implementation.

`paired` requires the patient term and exactly one sample for each patient in both contrast groups. Incomplete pairs/repeated biopsies are rejected rather than silently discarded or collapsed. `independent` excludes patient blocking and rejects repeated patient IDs unless `repeated_patients: explicit_unpaired` and a recorded `unpaired_rationale` are provided. That option remains an explicitly different analysis. A confounded or rank-deficient design and a model without residual degrees of freedom stop. In particular, between-patient HP comparisons may be confounded with a patient fixed effect; do not repair rank by dropping terms automatically.

`fit_scope: per_comparison` fits each selected contrast population separately. `all_included` fits once to all included samples and requires a consistent model/reference configuration across comparisons. These scopes can change size factors, dispersion estimates, baseMean and adjusted pvalues; neither is selected as the manuscript default.

The original active notebook removes all-zero features and uses `~ ID + Type` with contrasts in the opposite order to some Tumor/Normal table captions. Generic DEG scripts permit externally defined group contrasts and optional additive metadata. Four table captions alone cannot establish which source fit, pairing, sample subset or direction produced the reported numbers.

#### Annotation and DEG

`gene_id_policy` explicitly preserves IDs or strips a terminal version suffix. Duplicate features after transformation must error or be explicitly summed at the count level. Different Ensembl features sharing a symbol are not silently collapsed. Optional mapping TSV requires `gene_id`, `symbol`, `entrez_id`; statistics remain indexed by original feature ID. Duplicate mapping rows error or use an explicit first-source-row policy. Missing symbols retain feature IDs.

WTS reads a supplied local mapping snapshot. Keyed mapping and explicit duplicate/label policies retain the needed Toolbox concepts; unused live Get_biomaRt/ENSG2SYM query/export code is not part of the runtime. No new mapping or reference version is invented. The source's source_legacy_integrated issue remains evidence: Toolbox.R:723 replaces nonduplicated integrated labels with Ensembl IDs. Confirm the intended snapshot/label policy before providing the mapping file; do not choose labels to match the final GENE column.

DESeq configuration explicitly records `test` (Wald), `fitType`, `sfType`, `minReplicatesForReplace`, results `alpha`, `pAdjustMethod`, `independentFiltering`, `cooksCutoff`. Parameters omitted by the original notebook are not declared to be historical values. `deg.padj_cutoff` and `deg.abs_log2fc_cutoff` define a separate selected table; full statistics are retained, including NA. No LFC shrinkage or pvalue replacement is added. The observed source selection is |LFC|>=1 and padj<=0.05; this does not prove the final table-selection rule.

Each comparison writes DEG_all.tsv, DEG_selected.tsv, DEG_upregulated.tsv, DEG_downregulated.tsv, metadata_used.tsv, design_matrix.tsv, contrast_used.json and deseq_dataset.rds. The separate UP/DOWN summaries add gene_type; the complete/selected files retain their previous columns and numerical values. Main Table 2 contains selected genes without a deterministic source selection rule; the code does not invent a top-gene list to match it.

#### GO/KEGG GSEA

`gsea.enabled` must be explicit. Rankings use configured log2FoldChange (source) or Wald stat (new exploratory choice), finite values, descending order with stable source-order ties. `all_tested` reflects the WTS source. `pvalue_lt` is the separate generic-script prefilter and requires a cutoff. Do not silently treat them as equivalent. KEGG mapping is matched by feature ID; unmapped and duplicate rank IDs are counted, with explicit duplicate policies. A positional mapping assignment in the generic script is not reused.

For each GO/KEGG branch, choose `source_api` or `local_term2gene`, give a reference label, and explicitly configure min/max gene-set size, pvalue cutoff, adjustment method, exponent, eps and backend seed policy. Version-specific settings such as `by`, `method`, `nPerm` must match the installed clusterProfiler API; unknown/incompatible settings fail normally, without a fallback that changes methods. The rewrite calls the real package statistical engine, not an invented GSEA implementation.

`source_api` GO additionally requires OrgDb, ontology and key type. KEGG requires organism, key type and internal/live-data choice. Live KEGG is not a frozen historical snapshot. Local TERM2GENE files require `term_id`, `gene_id` and TERM2NAME requires `term_id`, `description`; duplicates/missing descriptions error. GO ranking IDs must match supplied GO gene-set IDs and KEGG ranking IDs the configured Entrez mapping. Actual local reference files are hashed. Local snapshots can avoid live database changes, but equivalence to the source APIs must be independently established.

WTS.Rmd paired GO uses ALL/ENSEMBL, min 3, max 800, BH. Paired KEGG uses min 3, max 800, nPerm 10000 and `none`. Other HP-read blocks use BH, while generic host scripts use a pvalue<0.05 ranking subset, GO 10-500 and KEGG minimum 8 with adjustment defaults. The final S9-S12 values have adjusted columns different from raw pvalues. No automatic correction from none to BH is made. Database/API versions and complete tested universes remain unavailable.

Outputs: `go_ranking.tsv`, `kegg_ranking.tsv`, `GO_GSEA.tsv`, `KEGG_GSEA.tsv` and reference/ranking summaries.

### Microbiome

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

## Changes with possible analytical effects

| Change | Basis | Potential effect |
| --- | --- | --- |
| Relative configured paths, embedded missing-package setup, unique run directories, errors, input hashes/session/status | Newly written infrastructure | No intended change to a chosen calculation; prevents accidental hidden input/output changes |
| Named FORMAT fields and explicit sample columns | WES.R/Rmd parsing reconstructed | Supports reordered tags; refuses multiallelic/ambiguous inputs instead of positional assumptions |
| Explicit tumor role and inclusion metadata | New input requirement | Can change included records relative to a source that melted every genotype column; study inclusion must be confirmed |
| Three separate WES profiles | Observed contradictory source blocks | No final profile chosen; switching profiles is a method change |
| Keyed transcript deduplication | Replaces source table()/match() Cartesian enumeration | Avoids phantom NA rows; explicit first-source-row behavior may affect which annotation is shown |
| Installed maftools converter; cohort size from included tumors | Source conversion/plotting reconstructed | Package classification/transcript behavior and plotted denominator can differ from supplied undated helper/hard-coded cohort |
| Preserve specimen ID as source_sample_id before MAF conversion; verify genomic and sample identity afterward | Actual artificial conversion test exposed upstream file-origin sample_id collision | Prevents duplicate headers/identifier replacement; variant filtering and statistical values are unchanged |
| Oncoplot explicit selected-variant bar and wider label margins | New presentation settings; no callable territory was supplied | Plots variant counts, never an inferred mutations-per-megabase TMB; mutation-negative samples are retained in cohort denominator/count output but may be absent from the package matrix |
| No automatic specimen exclusion or positional HP annotation | New metadata requirement | Original exclusions/plot labels must be verified rather than inferred |
| Strict integer counts; explicit zero-gene policy | WTS active source plus validation | Rejects fractional inputs instead of source as.integer truncation; remove-all-zero is explicit |
| Explicit model/contrast/reference/fit population | Source paired and generic alternatives | Changes in pairing, design, fit scope or DESeq settings can change all statistics; no unpaired manuscript default |
| Local mapping snapshot and keyed Entrez correspondence | Selected Toolbox concepts/generic mapping reconstructed | Annotation version, duplicates and legacy integrated labels can change reporting/ranking; no online default |
| Configured GO/KEGG method, ranking prefilter and adjustment | Conflicting WTS/generic source blocks | Ranking/adjustment/gene-set/backend changes can change enrichment; no none-to-BH repair |
| Explicit taxon/rank/metric/sample map | Microbiome.Rmd reconstructed | Removes row/position assumptions; wrong original row or metadata correspondence would have changed analysis |
| Preserve imported abundance and upstream Z-scores | Original import-only calculation | No invented denominator/reference population |
| Optional count-derived abundance / separately named recomputed Z-score | Newly written configurable derivation | A new analysis unless author-confirmed; test fixtures only are not evidence of manuscript use |
| Missingness, zero variance, duplicate taxa and median NA policies | New validation plus explicit legacy option | Can change missing classification versus original NA-to-low rule; every choice recorded |

No original source, final Excel cell, pvalue, padj, variant list or method choice was altered to force numerical agreement. Original code/Excel hashes were recorded and checked again during packaging.

## Verification scope

| Verification | Recorded outcome | Interpretation |
| --- | --- | --- |
| R syntax and independent helper resolution | Passed for all three scripts | Performed during the preceding code revision; the present script bytes are identical |
| Synthetic domain assertions | 42 passed | Artificial genotype/count/taxon inputs and invalid-input cases |
| Python verification methods | 22 passed; 1 skipped | The skipped component was the unavailable GO/KEGG engine |
| Numerical regression against the preceding standalone revision | 71 TSV pairs agreed; 9 stored DESeq2 objects agreed in counts/size factors/dispersions/coefficients | Artificial inputs only; four objects are repeated exports of one common fit |
| GO/KEGG statistical-engine execution | Not performed | clusterProfiler unavailable; ranking/ID adapters checked separately |
| Actual study execution | Not performed | Required annotated/count/taxon inputs and verified metadata unavailable |
| Executed study outputs against final Excel | Not performed | No executed study outputs or verified alias/contrast mappings |
| Read-only final-table audit | All 16 sheets inspected; the 25 Table 2 entries agree with S1-S4 at reported precision | Table-to-table consistency, not analytical reproduction |

Validation used R 4.3.3, DESeq2 1.42.0, maftools 2.26.1 and dplyr 1.2.1. Original study package versions and reference snapshots remain unknown. Detailed logs, safe synthetic fixtures, comparison utilities and earlier documentation are preserved in the separate internal validation package; they are not needed for routine execution of these scripts. Original handoff hashes and public payload checks are rechecked during packaging. This revision changes repository presentation and documentation only.

Before a study run, confirm the WES inputs/normalization, profile, exclusions, converter behavior and case map; the WTS counts, full metadata, pairing/repeated samples, contrast direction, fit population, DESeq2 settings and mapping/reference snapshots; the S4 columns/test universe and GO/KEGG settings; and the microbiome export/database, units, denominator and upstream score reference. A final authorized microbiome table would be needed for numerical comparison.

## Code and data provenance

The undated supplied annovarToMaf.R and full Toolbox.R are not redistributed. MAF conversion calls the installed maftools implementation under that package's own terms; its measured version is recorded. No source code from other research repositories is copied into this project.

Original scripts, workbooks, private specimen correspondence and credentials are excluded from the public package. License, authors and exact citation metadata must be confirmed by the owner. This package has not been published to GitHub and has no assigned release or DOI.

API references used for implementation checks, without implying historical-version equivalence:

- [DESeq2 official vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [clusterProfiler GSEA implementation](https://github.com/YuLab-SMU/clusterProfiler/blob/devel/R/gseAnalyzer.R)
- [maftools ANNOVAR converter](https://github.com/PoisonAlien/maftools/blob/master/R/annovarToMaf.R)

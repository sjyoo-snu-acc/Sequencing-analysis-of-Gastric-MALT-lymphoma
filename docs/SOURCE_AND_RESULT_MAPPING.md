# Source and final-result mapping

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
| 03_selected_draft and 91_previous_notes | Reference only | Preliminary unpaired model and first-row taxon assumption are not author-verified methods |

Original line numbers refer to unmodified files in the handoff. No figures/manuscript were supplied; no figure layout, fusion result or prose claim has been independently checked.

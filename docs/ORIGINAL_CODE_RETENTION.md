# Original code and later changes

The analysis concepts and source code supplied by the user remain the basis of these scripts. This revision makes their object names, data flow and calculations visible in the analysis body. It does not assert that any unverified branch was the final manuscript method. Source ranges below refer to the untouched originals in 02_original_code.

| Module / original source | Visible retained structure | Calculations retained | Later handling kept separate |
| --- | --- | --- | --- |
| WES.R:9-90,131-190,222-265; WES.Rmd:111-171 | ANV/VEP → GML → BEF → FIN/FIN2/FIN3/FIN4 → tmp/MAF/before | AD-based VAF; all-read strand ratio; distinct DP/VAF/strand/canonical/consequence/COSMIC/artifact rules | Named-column/FORMAT parsing, explicit sample roles/aliases, error checks, safe deduplication and measured installed converter |
| WTS.Rmd:305-353 | VAR; for-loop; pre/des/med/obj/res/wrt/cut/UP/DOWN | Raw-count subset; supplied design; DESeqDataSetFromMatrix, DESeq and results contrast; directional DEG thresholds | Count/metadata/rank checks, explicit model/reference/fit population and keyed mapping; no guessed patient pairs |
| WTS.Rmd:354-443; selected Toolbox mapping concepts | gnl; gse; kk2 | LFC/stat ranking as explicitly configured, keyed ID correspondence and chosen GO/KEGG parameters | Reference snapshots and adjustment choice remain explicit; source none/BH conflict is not repaired by changing values |
| Microbiome.Rmd:51-94 | pavian_obj; H.pylori_reads_pct; H.pylori_Zscore_pct | Imported reads/abundance/scores and the descriptive median rule when selected | Explicit taxon/rank/sample selection replaces row/column-position assumptions; optional new score derivation has a separate name |

## WES stages

FIN contains the selected profile's non-COSMIC candidates; FIN2 retains that intermediate stage, following the original no-op block. FIN3 applies the same already-computed complete selection, including COSMIC. FIN4 applies the confirmed deduplication policy and receives the verified case aliases. The filter predicates are evaluated once and reused; stages do not maintain competing copies of thresholds. The exact names follow the first WES.R flow and are applied to whichever profile the researcher explicitly selects. The strict notebook and second-script profiles remain distinct.

## WTS objects and summaries

The fitted-model calculation was moved out of the helper and back into the main comparison loop. pre is the checked count matrix for the chosen fit population, des is the checked metadata, med is the unfitted DESeqDataSet and obj is the fitted object. res remains the DESeqResults object. wrt adds reporting IDs without changing statistics. cut/UP/DOWN use the original case_when/positive-negative pattern with the configured thresholds.

DEG_all.tsv and DEG_selected.tsv retain the previous columns/order/values. DEG_upregulated.tsv and DEG_downregulated.tsv restore the source's separate directional summaries; they add gene_type for reporting. UP and DOWN must be disjoint and their union must equal the selected genes. Enrichment ranking is built from the same complete tested rows, with the explicitly selected ranking/filter policy. No gene order or statistic is changed to fit the final tables.

## Microbiome objects and units

The saved main table continues to report relative_abundance as a fraction. H.pylori_reads_pct exposes the same quantities as percent for the original-style object name. H.pylori_Zscore_pct contains the preserved upstream relative-abundance Z-score, not a substituted recomputed score. Group ordering follows Normal/Tumor/After; membership still comes from explicit metadata. Existing quantitative files and reference calculations remain unchanged.

## Deliberate corrections retained

Original fragile behavior is not restored for style: positional specimen inference, treating the first row as H. pylori, truncating fractional counts, silent duplicate joins, the all-zero negative-indexing trap, hard-coded cohort size, automatic exclusions, credential-bearing paths and output overwriting remain corrected or explicit. Those changes are documented in INPUTS_AND_METHODS.md and CHANGELOG_AND_VALIDATION.md.

Input checks, package preparation, run records and file handling are later additions. The three files remain independently executable. Original files were not modified, original authors/license details were not invented, and installed maftools is used under its own provenance rather than republishing the supplied undated converter copy.

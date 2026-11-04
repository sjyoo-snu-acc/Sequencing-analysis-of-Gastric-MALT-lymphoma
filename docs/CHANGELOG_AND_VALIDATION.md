# Changes and validation scope

## Changes by effect

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
| ID-based Excel comparator and output privacy scan | Newly written validation | Reports differences, never rewrites columns/results; does not itself re-run analysis |

No original source, final Excel cell, pvalue, padj, variant list or method choice was altered to force numerical agreement. Original code/Excel hashes were recorded and checked again during packaging.

## Validation levels

The current report is docs/validation/VALIDATION_REPORT.json. Three R files parse and execute independently; 42 domain assertions pass; 22 Python test methods pass and one GO/KEGG engine check is skipped. All 71 TSV comparison pairs agree with the preceding release, including alternative WES/post-filter/zero-result branches and a common WTS fit. Sixty-three helper-expression trees remain unchanged; three differ only in an error label. Nine stored DESeqDataSet objects retain identical count/size-factor/dispersion/coefficient results. These are artificial-input code checks, not manuscript reproduction.

| Level | Meaning | Scope |
| --- | --- | --- |
| 1. R parsing | Actual R parse of release R files | Does not evaluate package statistical functions |
| 2. Static path/dependency/privacy checks | Entry points, definitions/imports, hard-coded paths and public payload inspected | Does not grant data/code rights or prove clinical metadata |
| 3. Unit/artificial CLI tests | Known artificial values, invalid inputs, paired/unpaired distinctions, selected filter branches, sample/taxon/ID joins, no-overwrite behavior | Not study inputs or reported numbers |
| 4. Actual study module execution | Would require source counts/annotations/classifier outputs and metadata | NOT_RUN: source inputs absent |
| 5. Study output vs final Excel | Would require actual executed results and verified alias/contrast/ID mapping | NOT_RUN: no executed study outputs |

Read-only workbook audit: all 16 supplied sheets, titles/captions/headers, actual cell values and workbook structure were examined. Table 2's 25 displayed gene entries agree with S1-S4 within declared display-rounding tolerances. S4's raw/adjusted discrepancy and the KEGG adjustment inconsistency remain unresolved. This is table-to-table checking, not DESeq2/GSEA result reproduction. Private row-level evidence is in the separate internal validation package.

During implementation, tests exposed and corrected no-call placeholder handling, fixture directory setup, and the installed MAF converter's file-origin sample_id collision. These fixes affected code/fixtures, not study statistics. The final reported validation run, optional skipped components and their reasons are recorded without treating skipped tests as passes. The current files remain self-contained. This revision restores original-style analysis objects and direct DESeq2 calculation, and adds the source UP/DOWN summaries without changing the existing quantitative files. Progress output stays suppressed while warnings/errors remain visible. See STYLE_AND_EXECUTION_CHANGES.md for the current three-file execution layout. A prior test-only environment lookup correction did not change any analysis values.

## Remaining required study confirmations

- WES annotated input versions, biallelic normalization, tumor/normal sample list and verified case aliases; author-used filter chain, transcript deduplication and exclusions; supplied converter revision/classification equivalence.
- WTS raw-count provenance and scale, full metadata, patients/time points/repeated biopsies, actual fit population and design/reference levels; correct comparison direction; original DESeq defaults/version; ID mapping snapshots; selected-gene reporting rule.
- Source S4 statistic columns and test-universe export; actual GO/KEGG ranking/filter/size/backend/adjustment choices and reference snapshots.
- Microbiome classifier/export/database and species-rank semantics; metric column mapping, denominator, normalization and upstream Z-score reference; a corresponding authorized final microbiome table if numeric validation is requested.
- Author/institution code rights, project license and citation metadata, then actual publication URL/tag/archive identifier. No GitHub publishing action was taken.

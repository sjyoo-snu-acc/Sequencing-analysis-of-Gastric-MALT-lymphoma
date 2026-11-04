# Code consistency checks

Three R files parsed, their helper functions resolved independently, 42 artificial-input domain assertions passed, and 22 of 23 Python test methods passed. One GO/KEGG statistical-engine test was skipped because clusterProfiler was unavailable. Progress output remained empty in the executed analyses. Original-style objects, directional gene summaries, metadata errors and overwrite protection were checked.

The 38 existing TSV outputs/input-hash records (1362 rows) were checked against the preceding release by sample/variant, gene or taxon/sample identifiers. Five further branch cases were also executed against the preceding release: second and strict WES filters, the additional Ensembl exclusion, zero selected WES variants and an all-included WTS model. Those added 33 byte-identical TSV comparisons. All 71 pairs agreed; no value or analysis setting was changed to force agreement.

The parsed R expressions of 63 retained helper functions are unchanged. Three other helpers differ only in an error-label string. The fitted-model wrapper was moved into the visible WTS loop. Nine stored DESeqDataSet objects match in counts, size factors, dispersions and coefficients. Four of those objects are repeated exports of one all-included fit. The UP/DOWN tables are additional summaries and exactly partition the selected genes.

Tests use artificial input only. Actual study analysis and numerical agreement with the final workbooks remain unverified. The source S4/KEGG issues and unknown study settings remain documented. The installed dplyr version is measured for this validation, not asserted as an original study version. Real small-cohort maftools warnings remain in the execution logs.

# Style and execution changes

The supplied research scripts determine the style: short section labels, <- assignment, %>% processing, ordinary lists/data frames, loops and the original analysis objects. Generic introductory comments and the editable data-purpose/status/schema fields were removed from the front of the scripts. Necessary execution labels remain in run records, with an unconfirmed study run labelled exploratory by default.

WES exposes FIN–FIN4, tmp, MAF and before. WTS exposes the actual pre/des/med/obj/res fit in the comparison loop and restores wrt/cut/UP/DOWN and gse/kk2. Microbiome exposes the original H.pylori_reads_pct and H.pylori_Zscore_pct objects with explicit units and verified metadata. dplyr is used for the original pipe/classification idiom. Source contribution and added checks are listed in ORIGINAL_CODE_RETENTION.md.

Exactly three R files remain. Each includes its own required functions and missing-package preparation. No external R/config/project file or extra execution is needed. Progress/completion messages remain removed; return values, errors, statistical warnings and no-overwrite behavior are retained.

The file layout and presentation changed. Numeric analysis behavior under the same settings was checked against the preceding release, including alternative WES chains, the extra Ensembl exclusion, zero selected variants, paired WTS and all-included WTS fitting. Directional DEG tables are additional summaries, not altered selected-gene results. Current results are in validation/; actual study reproduction remains unverified.

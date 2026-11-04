# Gastric MALT lymphoma

R scripts for WES, WTS and H. pylori quantitation, adapted from the supplied analysis code.

| Script | Analysis | Main objects |
| --- | --- | --- |
| WES.R | Annotation, genotype filters, variant tables and oncoplot | ANV, VEP, GML, BEF, FIN–FIN4, tmp, MAF, before |
| WTS.R | DESeq2, DEG and GO/KEGG GSEA | CNT, MET, VAR, pre, des, med, obj, res, wrt, cut, UP, DOWN, gnl, gse, kk2 |
| Microbiome.R | H. pylori reads, relative abundance and Z-scores | pavian_obj, H.pylori, H.pylori_reads_pct, H.pylori_Zscore_pct |

Open the required file, edit Project.dir and the settings at the top, and click Source. Each file contains its required functions and package preparation. No other R script needs to be run. Missing packages are installed when Install.packages is TRUE; installed packages are not automatically updated. Progress messages are suppressed, while errors and calculation warnings remain visible.

The analysis blocks use the original object names, pipes and processing order. WTS now shows count selection, design, DESeqDataSetFromMatrix, DESeq and up/down DEG classification directly in the main loop. The original-source contribution and later input/output checks are distinguished in [ORIGINAL_CODE_RETENTION.md](docs/ORIGINAL_CODE_RETENTION.md).

## Inputs and results

Supply annotated WES tables, raw WTS counts and sample metadata, or a Pavian-like taxon table as appropriate. Fill the applicable NULL method settings using confirmed information. The final WES filter chain, study case aliases, WTS pairing/model and reference databases remain unconfirmed. The scripts stop when required information is absent. They do not substitute synthetic data automatically.

WTS comparisons are listed in VAR. GO and KEGG hold the enrichment settings. WES keeps the three source filter profiles distinct. Microbiome selects an explicit taxon/rank and preserves upstream scores; any recomputation requires its own stated denominator/reference.

Results are written to results/<module>/<Run.ID>/. The run directory also contains configuration_used.json, input_hashes.tsv, run_status.json and sessionInfo.txt. Existing results are protected. The original UP/DOWN DEG summaries are saved alongside the complete/selected DEG tables. Detailed input columns and outputs are described in [INPUTS_AND_METHODS.md](docs/INPUTS_AND_METHODS.md).

## Source and checks

[SOURCE_AND_RESULT_MAPPING.md](docs/SOURCE_AND_RESULT_MAPPING.md) records the final-sheet/source relationships. Source line numbers refer to the original handoff scripts. The originals and workbooks are preserved separately. Table 1 is clinical characteristics; no final quantitative microbiome sheet was supplied. S4 and KEGG adjustment discrepancies remain unresolved.

Current code, object, branch and numerical-regression checks are recorded in docs/validation/. They use artificial inputs; actual study analysis and manuscript numerical agreement remain unverified because the study inputs and mappings were not supplied. GO/KEGG statistical-engine availability is reported separately.

The archive contains three analysis R files, English documentation, artificial tests and optional audit tools. Original scripts, workbooks, real specimen data and credentials are excluded. Keep local edits containing private paths/identifiers out of public commits. License, authors, citation and publication identifiers remain unconfirmed; see [RIGHTS_AND_PROVENANCE.md](docs/RIGHTS_AND_PROVENANCE.md).

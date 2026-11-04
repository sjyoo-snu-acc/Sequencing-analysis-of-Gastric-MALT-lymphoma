# Rights and code provenance

No LICENSE or CITATION.cff has been assigned. The project authors, exact article metadata and permission to grant a repository license must be confirmed by the owner. The source handoff itself is not a public release authorization.

The supplied original analysis scripts are the source of the analysis methods, thresholds and research workflow. The current files adapt that work for explicit inputs and independent execution. ORIGINAL_CODE_RETENTION.md distinguishes the retained original structure from later validation, data/file handling and optional derivations. No claim is made that an unconfirmed branch was the manuscript-final code. No original script is copied wholesale. In particular, neither the full Toolbox nor the provided annovarToMaf source is redistributed.

The supplied annovarToMaf.R strongly resembles upstream maftools but has no provenance header identifying its revision or applicable redistribution terms. Upstream maftools DESCRIPTION identifies an MIT license plus its LICENSE file; that does not establish the exact history of the supplied local copy. The adapter calls an installed maftools package and records its actual version. Follow that package's own licensing and citation requirements. It is not vendored here.

Original source scripts contain local paths, specimen names and, in one notebook, a plaintext connection credential. Original source copies, private workbooks and specimen-level audit evidence are excluded. The internal finding identifies the file/line range without reproducing the credential. The credential owner should rotate it before reusing or sharing the historical source.

Reference documentation consulted for implementation/API checks, not for reconstructing unknown historical versions:

- [DESeq2 official vignette](https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [clusterProfiler author-maintained GSEA source](https://github.com/YuLab-SMU/clusterProfiler/blob/devel/R/gseAnalyzer.R)
- [maftools author-maintained converter](https://github.com/PoisonAlien/maftools/blob/master/R/annovarToMaf.R)
- [maftools package metadata](https://github.com/PoisonAlien/maftools/blob/master/DESCRIPTION)

Public upload, a GitHub URL, a tagged release, an archive DOI and citation details have not been created or verified by this task.

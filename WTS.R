# WTS

# Parameters ----
Project.dir<-getwd()
Install.packages<-TRUE
Run.ID<-paste0("WTS-", format(Sys.time(), "%Y%m%d-%H%M%S"))
Output.folder<-"results/WTS"

VAR<-list(
   list(
      id = "S1_S5_S9", factor = "group", numerator = NULL, denominator = NULL, design = NULL,
      model_kind = NULL, repeated_patients = "error", reference_levels = NULL
   ),
   list(
      id = "S2_S6_S10", factor = "group", numerator = NULL, denominator = NULL, design = NULL,
      model_kind = NULL, repeated_patients = "error", reference_levels = NULL
   ),
   list(
      id = "S3_S7_S11", factor = "group", numerator = NULL, denominator = NULL, design = NULL,
      model_kind = NULL, repeated_patients = "error", reference_levels = NULL
   ),
   list(
      id = "S4_S8_S12", factor = "group", numerator = NULL, denominator = NULL, design = NULL,
      model_kind = NULL, repeated_patients = "error", reference_levels = NULL
   )
)

GO<-list(
   engine = NULL,
   reference_label = NULL,
   orgdb = "org.Hs.eg.db",
   ontology = "ALL",
   key_type = "ENSEMBL",
   term2gene = NULL,
   term2name = NULL,
   parameters = list(
      minGSSize = NULL, maxGSSize = NULL, pvalueCutoff = NULL, pAdjustMethod = NULL, exponent = NULL,
      eps = NULL, seed = NULL
   )
)

KEGG<-list(
   engine = NULL,
   reference_label = NULL,
   organism = "hsa",
   key_type = "ncbi-geneid",
   use_internal_data = FALSE,
   term2gene = NULL,
   term2name = NULL,
   parameters = list(
      minGSSize = NULL, maxGSSize = NULL, pvalueCutoff = NULL, pAdjustMethod = NULL, exponent = NULL,
      eps = NULL, seed = NULL
   )
)

SET<-list(
   seed = 19970117,
   counts = "data/WTS/counts.tsv",
   metadata = "data/WTS/metadata.tsv",
   input_scale = "raw_integer_counts",
   count_annotation_columns = list(),
   gene_id_policy = NULL,
   duplicate_gene_policy = NULL,
   fit_scope = NULL,
   gene_mapping = "data/WTS/gene_mapping.tsv",
   mapping_duplicate_policy = NULL,
   annotation_label = NULL,
   prefilter = list(mode = NULL, minimum = NULL),
   deseq = list(test = "Wald", fitType = NULL, sfType = NULL, minReplicatesForReplace = NULL),
   results = list(alpha = NULL, pAdjustMethod = NULL, independentFiltering = NULL, cooksCutoff = NULL),
   deg = list(padj_cutoff = 0.05, abs_log2fc_cutoff = 1),
   comparisons = VAR,
   gsea = list(
      enabled = TRUE,
      ranking = list(
         metric = NULL, prefilter = NULL, pvalue_cutoff = NULL, duplicate_rank_ids = NULL,
         tie_policy = "stable_source_order"
      ),
      go = GO,
      kegg = KEGG
   )
)

# Checks and file I/O ----
abort <- function(...) stop(sprintf(...), call. = FALSE)

quiet_value <- function(expr) {
   value <- NULL
   invisible(utils::capture.output(value <- suppressMessages(force(expr)), type = "output"))
   value
}

is_bioc_package <- function(package) {
   package %in% c("DESeq2", "maftools", "clusterProfiler", "AnnotationDbi", "biomaRt") |
      grepl("^org\\..*\\.db$", package)
}

need_packages <- function(packages) {
   loaded <- function(package) suppressPackageStartupMessages(requireNamespace(package, quietly = TRUE))
   absent <- packages[!vapply(packages, loaded, logical(1))]
   if (!length(absent)) return(invisible(NULL))
   if (!Install.packages) abort("Missing R packages: %s. Set Install.packages = TRUE in this file and run it again.", paste(absent, collapse = ", "))
   old_repos <- getOption("repos")
   on.exit(options(repos = old_repos), add = TRUE)
   options(repos = c(CRAN = "https://cloud.r-project.org"))
   bioc <- absent[is_bioc_package(absent)]
   cran <- setdiff(absent, bioc)
   if (length(cran)) quiet_value(utils::install.packages(cran, repos = "https://cloud.r-project.org", quiet = TRUE))
   if (length(bioc)) {
      if (!loaded("BiocManager")) quiet_value(utils::install.packages("BiocManager", repos = "https://cloud.r-project.org", quiet = TRUE))
      if (!loaded("BiocManager")) abort("BiocManager installation failed. Check the R version and package repository connection.")
      quiet_value(BiocManager::install(bioc, ask = FALSE, update = FALSE, quiet = TRUE))
   }
   absent <- packages[!vapply(packages, loaded, logical(1))]
   if (length(absent)) abort("R package installation failed: %s. Check R/Bioconductor compatibility and package build requirements.", paste(absent, collapse = ", "))
   invisible(NULL)
}

required <- function(x, keys, context) {
  absent <- keys[vapply(keys, function(k) is.null(x[[k]]) || !length(x[[k]]) ||
    (is.character(x[[k]]) && any(!nzchar(x[[k]]) | is.na(x[[k]]))), logical(1))]
  if (length(absent)) abort("%s: unresolved/missing fields: %s", context, paste(absent, collapse = ", "))
  invisible(x)
}

columns <- function(x, keys, context) {
  if (anyDuplicated(names(x))) abort("%s has duplicate column names; supply an explicitly renamed input.", context)
  absent <- setdiff(keys, names(x))
  if (length(absent)) abort("%s: missing columns: %s", context, paste(absent, collapse = ", "))
}

unique_key <- function(x, context) {
  if (anyNA(x) || any(!nzchar(as.character(x))) || anyDuplicated(x)) abort("%s must be nonmissing and unique.", context)
}

choose <- function(x, choices, context) {
  if (length(x) != 1L || is.na(x) || !x %in% choices) abort("%s must be explicitly one of: %s", context, paste(choices, collapse = ", "))
  x
}

numeric_strict <- function(x, context, allow_na = FALSE) {
  text <- as.character(x)
  missing <- is.na(x) | text %in% c("", ".", "NA")
  value <- suppressWarnings(as.numeric(text))
  if (any(!missing & (is.na(value) | !is.finite(value)))) abort("%s contains invalid/nonfinite numbers.", context)
  value[missing] <- NA_real_
  if (!allow_na && anyNA(value)) abort("%s contains missing numbers.", context)
  value
}

read_tsv <- function(path, context = basename(path)) {
  if (!file.exists(path)) abort("Missing input: %s", context)
  widths <- count.fields(path, sep = "\t", quote = "", comment.char = "", blank.lines.skip = TRUE)
  if (!length(widths) || anyNA(widths) || any(widths != widths[1])) abort("%s contains inconsistent TSV row widths.", context)
  result <- read.delim(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                       check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = "NA")
  if (anyDuplicated(names(result))) abort("%s contains duplicate headers; positional repair is not supported.", context)
  result
}

repo_path <- function(root, path, context = "path") {
  required(list(path = path), "path", context)
  if (grepl("^(/|~|[A-Za-z]:|\\\\)", path)) abort("%s must be relative to Project.dir.", context)
  pieces <- strsplit(gsub("\\\\", "/", path), "/", fixed = TRUE)[[1]]
  if (".." %in% pieces) abort("%s may not escape Project.dir.", context)
  file.path(root, path)
}

write_tsv <- function(x, path) {
  if (file.exists(path)) abort("Refusing to overwrite file: %s", basename(path))
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, na = "NA")
}

validate_settings <- function(SET) {
   if (!is.logical(Install.packages) || length(Install.packages) != 1L || is.na(Install.packages)) abort("Install.packages must be TRUE or FALSE.")
   required(SET, c("module", "purpose", "method_status", "seed", "run_id", "output_dir"), "Analysis settings")
   choose(SET$purpose, c("analysis", "synthetic_test"), "input purpose")
   choose(SET$method_status, c("exploratory", "author_confirmed", "synthetic_only"), "method_status")
   if (SET$purpose == "synthetic_test" && SET$method_status != "synthetic_only") abort("Test data require synthetic_only status.")
   if (SET$purpose == "analysis" && SET$method_status == "synthetic_only") abort("Study analysis requires exploratory or author_confirmed status.")
   if (length(SET$seed) != 1L || !is.finite(SET$seed) || SET$seed != as.integer(SET$seed)) abort("seed must be an integer.")
   if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", SET$run_id)) abort("Run.ID must be a single safe directory name.")
   invisible(SET)
}

new_run <- function(root, cfg, input_paths) {
  directory <- file.path(repo_path(root, cfg$output_dir, "output_dir"), cfg$run_id)
  if (dir.exists(directory) || file.exists(directory)) abort("Run already exists: %s. Choose a new run_id.", cfg$run_id)
  if (!all(file.exists(input_paths))) abort("One or more configured inputs are missing.")
  if (!dir.create(directory, recursive = TRUE, showWarnings = FALSE)) abort("Cannot create output directory.")
  jsonlite::write_json(cfg, file.path(directory, "configuration_used.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  set.seed(as.integer(cfg$seed))
  hashes <- data.frame(input = names(input_paths), sha256 = vapply(input_paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)), stringsAsFactors = FALSE)
  write_tsv(hashes, file.path(directory, "input_hashes.tsv"))
  state <- list(module = cfg$module, purpose = cfg$purpose, method_status = cfg$method_status,
                run_id = cfg$run_id, seed = cfg$seed, started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                status = "RUNNING", manuscript_result_reproduction = "NOT_ASSESSED")
  jsonlite::write_json(state, file.path(directory, "run_status.json"), pretty = TRUE, auto_unbox = TRUE)
  list(directory = directory, state = state)
}

finish_run <- function(run, status, details = list()) {
  run$state$status <- status
  run$state$completed_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  run$state$details <- details
  jsonlite::write_json(run$state, file.path(run$directory, "run_status.json"), pretty = TRUE, auto_unbox = TRUE, null = "null")
  writeLines(capture.output(sessionInfo()), file.path(run$directory, "sessionInfo.txt"))
}

read_counts <- function(path, metadata, annotation_columns, id_policy, duplicate_policy) {
  raw <- read_tsv(path)
  columns(raw, c("gene_id", annotation_columns), "Raw counts")
  choose(id_policy, c("preserve", "strip_version"), "gene_id_policy")
  choose(duplicate_policy, c("error", "sum_counts"), "duplicate_gene_policy")
  sample_columns <- setdiff(names(raw), c("gene_id", annotation_columns))
  if (!setequal(sample_columns, metadata$sample_id)) abort("Count columns and metadata sample_id must match exactly, including explicitly excluded samples.")
  counts <- matrix(numeric_strict(as.matrix(raw[, metadata$sample_id, drop = FALSE]), "Raw counts"), nrow = nrow(raw),
                   dimnames = list(NULL, metadata$sample_id))
  if (!nrow(counts) || any(counts < 0 | counts != trunc(counts)) || any(counts > .Machine$integer.max))
    abort("DESeq2 input requires nonnegative raw integer counts. TPM/fractional expected counts are not truncated or rounded.")
  id <- raw$gene_id
  if (id_policy == "strip_version") id <- sub("\\.[0-9]+$", "", id)
  if (anyNA(id) || any(!nzchar(id))) abort("Missing gene IDs.")
  if (anyDuplicated(id) && duplicate_policy == "error") abort("Duplicate gene IDs after the configured ID transformation.")
  if (anyDuplicated(id)) counts <- rowsum(counts, group = id, reorder = FALSE) else rownames(counts) <- id
  if (any(counts > .Machine$integer.max)) abort("Summed counts exceed integer storage.")
  storage.mode(counts) <- "integer"
  counts
}

validate_wts_metadata <- function(meta) {
  columns(meta, c("sample_id", "patient_id", "type", "hp_status", "group", "include"), "WTS metadata")
  unique_key(meta$sample_id, "WTS sample_id")
  for (key in c("patient_id", "type", "hp_status", "group", "include")) {
    if (anyNA(meta[[key]]) || any(!nzchar(meta[[key]]))) abort("WTS metadata %s cannot be missing.", key)
  }
  if (any(!meta$type %in% c("Normal", "Tumor", "After")) || any(!meta$hp_status %in% c("Positive", "Negative", "Unknown")))
    abort("Use explicit Normal/Tumor/After and Positive/Negative/Unknown metadata labels.")
  if (any(!meta$include %in% c("TRUE", "FALSE"))) abort("WTS include must be TRUE or FALSE.")
  expected <- paste0(ifelse(meta$hp_status == "Positive", "HPpos_", ifelse(meta$hp_status == "Negative", "HPneg_", "HPunknown_")), meta$type)
  standard <- grepl("^HP(pos|neg|unknown)_", meta$group)
  if (any(standard & meta$group != expected)) abort("HP-stratified group label conflicts with explicit type/hp_status metadata.")
  invisible(meta)
}

prepare_model <- function(meta, contrast) {
  required(contrast, c("id", "factor", "numerator", "denominator", "design", "model_kind", "repeated_patients", "reference_levels"), "WTS comparison")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", contrast$id)) abort("Invalid comparison id.")
  choose(contrast$model_kind, c("paired", "independent"), "model_kind")
  choose(contrast$repeated_patients, c("error", "explicit_unpaired"), "repeated_patients")
  # The supplied sources support additive fixed-effect models only.
  if (!grepl("^~\\s*[A-Za-z][A-Za-z0-9_]*(\\s*\\+\\s*[A-Za-z][A-Za-z0-9_]*)*\\s*$", contrast$design))
    abort("design must be an explicit additive formula such as ~ patient_id + group.")
  formula <- as.formula(contrast$design)
  vars <- all.vars(formula)
  columns(meta, c(vars, contrast$factor), "Model metadata")
  if (!contrast$factor %in% vars || contrast$numerator == contrast$denominator) abort("Contrast factor must be in design with distinct numerator/denominator.")
  if (contrast$model_kind == "paired" && !"patient_id" %in% vars) abort("Paired model must include patient_id.")
  if (contrast$model_kind == "independent" && "patient_id" %in% vars) abort("Independent model cannot silently include patient_id blocking.")
  if (contrast$model_kind == "paired") {
    pair_meta <- meta[meta[[contrast$factor]] %in% c(contrast$numerator, contrast$denominator), , drop = FALSE]
    pairs <- table(pair_meta$patient_id, pair_meta[[contrast$factor]])
    if (ncol(pairs) != 2L || !all(pairs == 1L)) abort("Paired comparison requires exactly one sample per patient in each contrast group. Incomplete pairs/replicates require a separate confirmed model.")
  } else if (anyDuplicated(meta$patient_id)) {
    if (contrast$repeated_patients != "explicit_unpaired") abort("Repeated patients detected in independent fit. Confirm model or explicitly record unpaired treatment.")
    required(contrast, "unpaired_rationale", "Explicit unpaired treatment of repeated patients")
  }
  numeric_vars <- unlist(contrast$numeric_covariates, use.names = FALSE)
  if (length(setdiff(numeric_vars, vars))) abort("numeric_covariates must occur in design.")
  for (v in vars) {
    if (v %in% numeric_vars) meta[[v]] <- numeric_strict(meta[[v]], paste("Covariate", v)) else {
      required(contrast$reference_levels, v, paste("reference_levels for", v))
      values <- meta[[v]]
      if (anyNA(values) || any(!nzchar(values))) abort("Missing design factor: %s", v)
      ref <- contrast$reference_levels[[v]]
      if (!ref %in% values) abort("Reference level for %s is absent from fit population.", v)
      meta[[v]] <- factor(values, levels = c(ref, sort(setdiff(unique(values), ref))))
    }
  }
  if (!all(c(contrast$numerator, contrast$denominator) %in% levels(meta[[contrast$factor]]))) abort("One or both contrast groups are absent.")
  matrix <- model.matrix(formula, data = meta)
  if (qr(matrix)$rank < ncol(matrix)) abort("Design matrix is not full rank; do not resolve confounding by dropping patient/group terms automatically.")
  if (nrow(matrix) <= ncol(matrix)) abort("Model has no residual degrees of freedom.")
  rownames(meta) <- meta$sample_id
  list(metadata = meta, formula = formula, matrix = matrix)
}

read_gene_mapping <- function(path, id_policy, policy) {
  x <- read_tsv(path)
  columns(x, c("gene_id", "symbol", "entrez_id"), "Gene mapping")
  if (id_policy == "strip_version") x$gene_id <- sub("\\.[0-9]+$", "", x$gene_id)
  choose(policy, c("error", "first_source_row"), "mapping_duplicate_policy")
  if (anyNA(x$gene_id) || any(!nzchar(x$gene_id))) abort("Gene mapping requires nonmissing source IDs.")
  if (anyDuplicated(x$gene_id) && policy == "error") abort("One-to-many gene mapping requires an explicit policy.")
  x[!duplicated(x$gene_id), , drop = FALSE]
}

make_ranking <- function(result, settings, mapping = NULL, target = "gene_id") {
  required(settings, c("metric", "prefilter", "duplicate_rank_ids", "tie_policy"), "GSEA ranking")
  choose(settings$metric, c("log2FoldChange", "stat"), "ranking metric")
  choose(settings$prefilter, c("all_tested", "pvalue_lt"), "GSEA prefilter")
  choose(settings$duplicate_rank_ids, c("error", "first_source_row"), "duplicate_rank_ids")
  choose(settings$tie_policy, "stable_source_order", "tie_policy")
  keep <- is.finite(result[[settings$metric]])
  if (settings$prefilter == "pvalue_lt") {
    required(settings, "pvalue_cutoff", "pvalue_lt ranking")
    if (settings$pvalue_cutoff <= 0 || settings$pvalue_cutoff > 1) abort("Invalid ranking pvalue cutoff.")
    keep <- keep & !is.na(result$pvalue) & result$pvalue < settings$pvalue_cutoff
  }
  r <- result[keep, , drop = FALSE]
  ids <- r$gene_id
  if (target == "entrez_id") {
    if (is.null(mapping)) abort("KEGG ranking requires an explicit gene mapping snapshot.")
    ids <- mapping$entrez_id[match(r$gene_id, mapping$gene_id)]
  }
  valid <- !is.na(ids) & !ids %in% c("", ".", "NA")
  dropped_mapping <- sum(!valid)
  ids <- ids[valid]
  values <- r[[settings$metric]][valid]
  if (anyDuplicated(ids) && settings$duplicate_rank_ids == "error") abort("Duplicate GSEA ranking identifiers after mapping.")
  keep_unique <- !duplicated(ids)
  ranking <- values[keep_unique]
  names(ranking) <- ids[keep_unique]
  ranking <- ranking[order(-ranking, seq_along(ranking))]
  if (length(ranking) < 2L) abort("Too few finite, mapped genes for GSEA.")
  attr(ranking, "mapping_dropped") <- dropped_mapping
  attr(ranking, "duplicate_ids_dropped") <- sum(!keep_unique)
  ranking
}

validate_gsea <- function(root, g) {
  required(g, "enabled", "GSEA")
  if (!is.logical(g$enabled) || length(g$enabled) != 1L) abort("gsea.enabled must be logical.")
  if (!g$enabled) return(character())
  required(g, c("ranking", "go", "kegg"), "GSEA")
  inputs <- character()
  for (kind in c("go", "kegg")) {
    part <- g[[kind]]
    required(part, c("engine", "reference_label", "parameters"), paste("GSEA", kind))
    choose(part$engine, c("source_api", "local_term2gene"), paste(kind, "engine"))
    required(part$parameters, c("minGSSize", "maxGSSize", "pvalueCutoff", "pAdjustMethod", "exponent", "eps", "seed"), paste(kind, "parameters"))
    choose(part$parameters$pAdjustMethod, stats::p.adjust.methods, paste(kind, "pAdjustMethod"))
    if (part$parameters$minGSSize < 1 || part$parameters$maxGSSize < part$parameters$minGSSize ||
        part$parameters$pvalueCutoff <= 0 || part$parameters$pvalueCutoff > 1 || part$parameters$eps < 0 || part$parameters$exponent < 0) abort("Invalid GSEA parameters.")
    if (part$engine == "local_term2gene") {
      required(part, c("term2gene", "term2name"), paste(kind, "local reference"))
      inputs[paste0(kind, "_term2gene")] <- repo_path(root, part$term2gene)
      inputs[paste0(kind, "_term2name")] <- repo_path(root, part$term2name)
    } else if (kind == "go") required(part, c("orgdb", "ontology", "key_type"), "GO source API") else
      required(part, c("organism", "key_type", "use_internal_data"), "KEGG source API")
  }
  need_packages("clusterProfiler")
  inputs
}

gsea_one <- function(kind, ranking, part, inputs) {
  args <- part$parameters
  args$geneList <- ranking
  args$verbose <- FALSE
  if (part$engine == "local_term2gene") {
    t2g <- read_tsv(inputs[paste0(kind, "_term2gene")])
    t2n <- read_tsv(inputs[paste0(kind, "_term2name")])
    columns(t2g, c("term_id", "gene_id"), "TERM2GENE")
    columns(t2n, c("term_id", "description"), "TERM2NAME")
    if (anyNA(t2g$term_id) || anyNA(t2g$gene_id) || any(!nzchar(t2g$term_id)) || any(!nzchar(t2g$gene_id))) abort("TERM2GENE contains missing identifiers.")
    unique_key(t2n$term_id, "TERM2NAME term_id")
    if (anyDuplicated(paste(t2g$term_id, t2g$gene_id))) abort("TERM2GENE contains duplicate term/gene pairs.")
    if (!all(t2g$term_id %in% t2n$term_id)) abort("TERM2NAME must describe every supplied gene set.")
    args$TERM2GENE <- t2g[, c("term_id", "gene_id")]
    args$TERM2NAME <- t2n[, c("term_id", "description")]
    answer <- do.call(clusterProfiler::GSEA, args)
  } else if (kind == "go") {
    need_packages(part$orgdb)
    args$OrgDb <- get(part$orgdb, asNamespace(part$orgdb))
    args$ont <- part$ontology
    args$keyType <- part$key_type
    answer <- do.call(clusterProfiler::gseGO, args)
  } else {
    args$organism <- part$organism
    args$keyType <- part$key_type
    args$use_internal_data <- part$use_internal_data
    answer <- do.call(clusterProfiler::gseKEGG, args)
  }
  if (is.null(answer)) abort("GSEA returned NULL; no result object was produced.")
  as.data.frame(answer)
}

validate_wts_config <- function(root, cfg) {
  required(cfg, c("counts", "metadata", "input_scale", "gene_id_policy", "duplicate_gene_policy",
    "fit_scope", "comparisons", "prefilter", "deseq", "results", "deg", "gsea", "mapping_duplicate_policy", "annotation_label"), "WTS configuration")
  choose(cfg$input_scale, "raw_integer_counts", "input_scale")
  if (!"count_annotation_columns" %in% names(cfg)) abort("Specify count_annotation_columns, using [] when there are none.")
  choose(cfg$fit_scope, c("per_comparison", "all_included"), "fit_scope")
  if (!is.list(cfg$comparisons) || !length(cfg$comparisons)) abort("At least one fully specified comparison is required.")
  ids <- vapply(cfg$comparisons, function(c) { required(c, c("id", "design", "model_kind", "factor", "numerator", "denominator", "reference_levels", "repeated_patients"), "comparison")
  c$id }, character(1))
  unique_key(ids, "Comparison id")
  required(cfg$prefilter, "mode", "Count prefilter")
  choose(cfg$prefilter$mode, c("remove_all_zero", "sum_at_least"), "prefilter.mode")
  if (cfg$prefilter$mode == "sum_at_least") required(cfg$prefilter, "minimum", "Count prefilter")
  required(cfg$deseq, c("test", "fitType", "sfType", "minReplicatesForReplace"), "DESeq settings")
  choose(cfg$deseq$test, "Wald", "test")
  choose(cfg$deseq$fitType, c("parametric", "local", "mean"), "fitType")
  choose(cfg$deseq$sfType, c("ratio", "poscounts", "iterate"), "sfType")
  required(cfg$results, c("alpha", "pAdjustMethod", "independentFiltering", "cooksCutoff"), "results settings")
  required(cfg$deg, c("padj_cutoff", "abs_log2fc_cutoff"), "DEG selection")
  if (cfg$results$alpha <= 0 || cfg$results$alpha >= 1 || cfg$deg$padj_cutoff <= 0 || cfg$deg$padj_cutoff >= 1 || cfg$deg$abs_log2fc_cutoff < 0) abort("Invalid DEG/alpha settings.")
  choose(cfg$results$pAdjustMethod, stats::p.adjust.methods, "DESeq pAdjustMethod")
  paths <- c(counts = repo_path(root, cfg$counts), metadata = repo_path(root, cfg$metadata))
  if (!is.null(cfg$gene_mapping)) paths <- c(paths, gene_mapping = repo_path(root, cfg$gene_mapping))
  paths <- c(paths, validate_gsea(root, cfg$gsea))
  if (cfg$gsea$enabled && !"gene_mapping" %in% names(paths)) abort("Enabled KEGG GSEA requires a mapping snapshot.")
  need_packages("DESeq2")
  paths
}

# Packages ----
SET$module<-"WTS"
if (is.null(SET$purpose)) SET$purpose<-"analysis"
if (is.null(SET$method_status)) SET$method_status<-"exploratory"
SET$run_id<-Run.ID
SET$output_dir<-Output.folder
validate_settings(SET)
Input.files<-validate_wts_config(Project.dir, SET)
need_packages(c("jsonlite", "digest", "dplyr"))
suppressPackageStartupMessages(library(dplyr))
Run<-new_run(Project.dir, SET, Input.files)
OUT<-Run$directory

# Analysis ----
tryCatch({
   # Count and sample information
   MET<-read_tsv(Input.files["metadata"])
   validate_wts_metadata(MET)
   CNT<-read_counts(Input.files["counts"], MET, unlist(SET$count_annotation_columns), SET$gene_id_policy, SET$duplicate_gene_policy)
   MET<-MET[MET$include == "TRUE", , drop = FALSE]
   CNT<-CNT[, MET$sample_id, drop = FALSE]
   if (!nrow(MET)) abort("No included WTS samples.")
   TRANSFORM<-if ("gene_mapping" %in% names(Input.files)) read_gene_mapping(Input.files["gene_mapping"], SET$gene_id_policy, SET$mapping_duplicate_policy) else NULL
   VAR<-SET$comparisons

   if (SET$fit_scope == "all_included") {
      signature<-function(x) x[c("design", "model_kind", "factor", "reference_levels", "numeric_covariates", "repeated_patients", "unpaired_rationale")]
      if (!all(vapply(VAR, function(x) identical(signature(x), signature(VAR[[1]])), logical(1)))) {
         abort("all_included requires the same model and factor/reference configuration for every contrast.")
      }
      for (var in VAR) prepare_model(MET, var)
   }

   # DESeq2
   obj<-NULL
   totals<-list()
   for (i in seq_along(VAR)) {
      var<-VAR[[i]]
      subdir<-file.path(OUT, var$id)
      dir.create(subdir)

      if (SET$fit_scope == "per_comparison" || is.null(obj)) {
         des<-if (SET$fit_scope == "per_comparison") MET[MET[[var$factor]] %in% c(var$numerator, var$denominator), , drop = FALSE] else MET
         MODEL<-prepare_model(des, var)
         des<-MODEL$metadata
         pre<-CNT[, des$sample_id, drop = FALSE]
         keep<-if (SET$prefilter$mode == "remove_all_zero") rowSums(pre) > 0 else rowSums(pre) >= SET$prefilter$minimum
         if (!any(keep)) abort("No genes pass the configured count filter.")
         pre<-pre[keep, , drop = FALSE]
         if (any(colSums(pre) == 0)) abort("An included sample has zero library size after filtering.")

         med<-DESeq2::DESeqDataSetFromMatrix(countData = pre, colData = des, design = MODEL$formula)
         obj<-quiet_value(DESeq2::DESeq(med, test = SET$deseq$test, fitType = SET$deseq$fitType,
            sfType = SET$deseq$sfType, minReplicatesForReplace = SET$deseq$minReplicatesForReplace,
            quiet = TRUE, parallel = FALSE))
      }

      args<-SET$results
      args$object<-obj
      args$contrast<-c(var$factor, var$numerator, var$denominator)
      res<-do.call(DESeq2::results, args)
      wrt<-res %>% as.data.frame()
      wrt$gene_id<-rownames(wrt)
      wrt$GENE<-wrt$gene_id
      if (!is.null(TRANSFORM)) {
         match_id<-match(wrt$gene_id, TRANSFORM$gene_id)
         symbol<-TRANSFORM$symbol[match_id]
         valid<-!is.na(symbol) & !symbol %in% c("", ".", "NA")
         wrt$GENE[valid]<-symbol[valid]
      }
      wrt<-wrt[, c("gene_id", "GENE", setdiff(names(wrt), c("gene_id", "GENE"))), drop = FALSE]

      # Upregulated / downregulated genes
      cut<-wrt %>%
         mutate(gene_type = case_when(
            log2FoldChange >= SET$deg$abs_log2fc_cutoff & padj <= SET$deg$padj_cutoff ~ "up",
            log2FoldChange <= -SET$deg$abs_log2fc_cutoff & padj <= SET$deg$padj_cutoff ~ "down",
            TRUE ~ "ns"))
      UP<-cut %>% filter(gene_type == "up")
      DOWN<-cut %>% filter(gene_type == "down")
      significant<-cut$gene_type != "ns"

      write_tsv(wrt, file.path(subdir, "DEG_all.tsv"))
      write_tsv(wrt[significant, , drop = FALSE], file.path(subdir, "DEG_selected.tsv"))
      write_tsv(UP, file.path(subdir, "DEG_upregulated.tsv"))
      write_tsv(DOWN, file.path(subdir, "DEG_downregulated.tsv"))
      write_tsv(des, file.path(subdir, "metadata_used.tsv"))
      write_tsv(data.frame(sample_id = rownames(MODEL$matrix), MODEL$matrix, check.names = FALSE), file.path(subdir, "design_matrix.tsv"))
      saveRDS(obj, file.path(subdir, "deseq_dataset.rds"))
      jsonlite::write_json(list(factor = var$factor, numerator = var$numerator, denominator = var$denominator,
         log2FoldChange = "log2(numerator / denominator)", design = var$design, model_kind = var$model_kind,
         fit_scope = SET$fit_scope, reference_levels = var$reference_levels),
         file.path(subdir, "contrast_used.json"), pretty = TRUE, auto_unbox = TRUE)

      # GSEA: GO / KEGG
      gsea_status<-"DISABLED_BY_EXPLICIT_CONFIGURATION"
      if (SET$gsea$enabled) {
         for (kind in c("go", "kegg")) {
            gnl<-make_ranking(cut, SET$gsea$ranking, TRANSFORM, if (kind == "kegg") "entrez_id" else "gene_id")
            write_tsv(data.frame(gene_id = names(gnl), ranking = as.numeric(gnl)), file.path(subdir, paste0(kind, "_ranking.tsv")))
            if (kind == "go") {
               gse<-quiet_value(gsea_one(kind, gnl, SET$gsea[[kind]], Input.files))
               write_tsv(gse, file.path(subdir, "GO_GSEA.tsv"))
            } else {
               kk2<-quiet_value(gsea_one(kind, gnl, SET$gsea[[kind]], Input.files))
               write_tsv(kk2, file.path(subdir, "KEGG_GSEA.tsv"))
            }
            jsonlite::write_json(list(reference_label = SET$gsea[[kind]]$reference_label,
               engine = SET$gsea[[kind]]$engine, mapped_genes = length(gnl), mapping_dropped = attr(gnl, "mapping_dropped"),
               duplicate_ids_dropped = attr(gnl, "duplicate_ids_dropped"), tied_statistics = sum(duplicated(as.numeric(gnl)))),
               file.path(subdir, paste0(kind, "_reference_and_ranking.json")), pretty = TRUE, auto_unbox = TRUE)
         }
         gsea_status<-"COMPLETED_WITH_CONFIGURED_REFERENCES"
      }
      totals[[var$id]]<-list(tested_rows = nrow(wrt), selected_rows = sum(significant), fit_samples = nrow(des), gsea = gsea_status)
   }
   Details<-list(comparisons = totals, annotation_label = SET$annotation_label, original_environment = "UNKNOWN")
   finish_run(Run, "COMPLETED", Details)
}, error = function(e) {
   finish_run(Run, "FAILED", list(error = conditionMessage(e)))
   stop(e)
})

# Microbiome

# Parameters ----
Project.dir<-getwd()
Install.packages<-TRUE
Run.ID<-paste0("Microbiome-", format(Sys.time(), "%Y%m%d-%H%M%S"))
Output.folder<-"results/Microbiome"

SET<-list(
   seed = 19970117, taxa_table = "data/Microbiome/pavian.tsv",
   column_map = "data/Microbiome/metric_columns.tsv", metadata = "data/Microbiome/metadata.tsv",
   input_tool = NULL, input_database_label = NULL,
   taxon_columns = list(id = NULL, name = NULL, rank = NULL), taxonomic_rank = "species",
   input_rank_value = NULL, target = list(id = NULL, name = "Helicobacter pylori"),
   duplicate_taxa = "error", missing_taxon = "error", relative_mode = NULL, relative_units = NULL,
   relative_denominator_description = NULL, denominator = NULL, zero_denominator = "error",
   zscore_mode = NULL, upstream_zscore_reference = NULL, zscore_reference_samples = NULL,
   zscore_metric = NULL, zscore_transform = NULL, zscore_na_policy = NULL, zero_variance = "error",
   median_mode = NULL, median_reference_samples = NULL, median_na_policy = "preserve"
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

validate_microbiome_config <- function(root, cfg) {
  required(cfg, c("taxa_table", "column_map", "metadata", "input_tool", "input_database_label", "taxon_columns",
    "taxonomic_rank", "input_rank_value", "target", "duplicate_taxa", "missing_taxon", "relative_mode",
    "relative_units", "relative_denominator_description", "zero_denominator", "zscore_mode",
    "median_mode", "median_na_policy"), "Microbiome configuration")
  required(cfg$taxon_columns, c("id", "name", "rank"), "Taxon columns")
  if (is.null(cfg$target$id) && is.null(cfg$target$name)) abort("Target requires an exact taxon ID or species name.")
  choose(cfg$duplicate_taxa, c("error", "sum_counts_only"), "duplicate_taxa")
  choose(cfg$missing_taxon, c("error", "confirmed_zero"), "missing_taxon")
  if (cfg$missing_taxon == "confirmed_zero") required(cfg, "nondetection_basis", "Confirmed zero policy")
  choose(cfg$relative_mode, c("imported", "count_over_denominator", "none"), "relative_mode")
  choose(cfg$relative_units, c("percent", "fraction"), "relative_units")
  choose(cfg$zero_denominator, c("error", "NA"), "zero_denominator")
  choose(cfg$zscore_mode, c("imported", "recompute", "none"), "zscore_mode")
  choose(cfg$median_mode, c("none", "all_included", "explicit_reference"), "median_mode")
  choose(cfg$median_na_policy, c("preserve", "legacy_down"), "median_na_policy")
  inputs <- c(taxa = repo_path(root, cfg$taxa_table), column_map = repo_path(root, cfg$column_map), metadata = repo_path(root, cfg$metadata))
  if (cfg$relative_mode == "count_over_denominator") {
    required(cfg, "denominator", "Count-based relative abundance")
    inputs <- c(inputs, denominator = repo_path(root, cfg$denominator))
  }
  if (cfg$zscore_mode == "imported") required(cfg, "upstream_zscore_reference", "Imported Z-score provenance (explicit unknown is allowed)")
  if (cfg$zscore_mode == "recompute") {
    required(cfg, c("zscore_reference_samples", "zscore_metric", "zscore_transform", "zscore_na_policy", "zero_variance"), "Recomputed Z-score")
    choose(cfg$zscore_metric, c("count", "relative_abundance"), "zscore_metric")
    choose(cfg$zscore_transform, c("identity", "log2p1"), "zscore_transform")
    choose(cfg$zscore_na_policy, c("error", "exclude_reference_NA"), "zscore_na_policy")
    choose(cfg$zero_variance, c("error", "NA"), "zero_variance")
    inputs <- c(inputs, zscore_reference = repo_path(root, cfg$zscore_reference_samples))
  }
  if (cfg$median_mode == "explicit_reference") {
    required(cfg, "median_reference_samples", "Median reference")
    inputs <- c(inputs, median_reference = repo_path(root, cfg$median_reference_samples))
  }
  inputs
}

select_taxon <- function(taxa, cfg) {
  col <- cfg$taxon_columns
  columns(taxa, unlist(col, use.names = FALSE), "Taxa table")
  rank <- taxa[[col$rank]] == cfg$input_rank_value & !is.na(taxa[[col$rank]])
  id_match <- if (!is.null(cfg$target$id)) !is.na(taxa[[col$id]]) & taxa[[col$id]] == as.character(cfg$target$id) else rep(TRUE, nrow(taxa))
  name_match <- if (!is.null(cfg$target$name)) !is.na(taxa[[col$name]]) & trimws(taxa[[col$name]]) == cfg$target$name else rep(TRUE, nrow(taxa))
  if (!is.null(cfg$target$id) && !is.null(cfg$target$name) && any(rank & id_match & !name_match)) abort("Configured taxon ID and exact species name disagree in input.")
  hit <- which(rank & id_match & name_match)
  if (!length(hit) && cfg$missing_taxon == "error") abort("Target taxon absent at the configured rank; absence is not assumed to be zero.")
  if (length(hit) > 1L && cfg$duplicate_taxa == "error") abort("Duplicate target taxon rows; select a documented aggregation policy.")
  hit
}

extract_taxon_metrics <- function(taxa, map, meta, cfg) {
  columns(map, c("input_column", "sample_id", "metric"), "Metric column map")
  unique_key(map$input_column, "Metric input_column")
  unique_key(paste(map$sample_id, map$metric), "Sample/metric combination")
  if (any(!map$metric %in% c("count", "relative_abundance", "zscore_count", "zscore_relative"))) abort("Unsupported metric name in column map.")
  if (!setequal(map$sample_id, meta$sample_id)) abort("Every metadata sample needs mapped metric columns; no unlisted samples are allowed.")
  columns(taxa, map$input_column, "Taxa metric columns")
  hit <- select_taxon(taxa, cfg)
  if (length(hit) > 1L && (any(map$metric != "count") || cfg$relative_mode == "imported" || cfg$zscore_mode == "imported"))
    abort("Duplicate taxa may be summed only for counts. Imported fractions/Z-scores cannot be aggregated as counts.")
  result <- meta
  result$taxon_id <- if (!length(hit)) as.character(cfg$target$id %||% NA_character_) else taxa[[cfg$taxon_columns$id]][hit[1]]
  result$taxon_name <- if (!length(hit)) as.character(cfg$target$name %||% NA_character_) else taxa[[cfg$taxon_columns$name]][hit[1]]
  result$taxonomic_rank <- cfg$taxonomic_rank
  for (metric in c("count", "relative_abundance", "zscore_count", "zscore_relative")) result[[metric]] <- NA_real_
  for (j in seq_len(nrow(map))) {
    value <- if (!length(hit)) { if (map$metric[j] == "count") 0 else NA_real_ } else {
      v <- numeric_strict(taxa[[map$input_column[j]]][hit], paste("Metric", map$metric[j]), allow_na = TRUE)
      if (anyNA(v)) NA_real_ else sum(v)
    }
    result[[map$metric[j]]][match(map$sample_id[j], result$sample_id)] <- value
  }
  if (any(result$count < 0 | result$count != trunc(result$count), na.rm = TRUE)) abort("Mapped count must be a nonnegative integer read count, not transformed abundance.")
  cap <- if (cfg$relative_units == "percent") 100 else 1
  if (any(result$relative_abundance < 0 | result$relative_abundance > cap, na.rm = TRUE)) abort("Imported relative abundance is outside its declared units.")
  result
}

`%||%` <- function(x, y) if (is.null(x)) y else x

reference_ids <- function(path, known_samples, context) {
  ref <- read_tsv(path)
  columns(ref, "sample_id", context)
  unique_key(ref$sample_id, context)
  if (!nrow(ref) || !all(ref$sample_id %in% known_samples)) abort("%s must be a nonempty subset of included sample_id values.", context)
  ref$sample_id
}

compute_zscore <- function(values, reference_positions, transform, na_policy, variance_policy) {
  x <- if (transform == "log2p1") log2(values + 1) else values
  ref <- x[reference_positions]
  if (anyNA(ref) && na_policy == "error") abort("Missing Z-score reference values.")
  ref <- ref[!is.na(ref)]
  if (length(ref) < 2L) abort("Z-score needs at least two nonmissing reference values.")
  center <- mean(ref)
  spread <- stats::sd(ref)
  if (spread == 0) {
    if (variance_policy == "error") abort("Z-score reference has zero variance.")
    value <- rep(NA_real_, length(x))
  } else value <- (x - center) / spread
  list(value = value, reference_n = length(ref), mean = center, sd = spread)
}

# Packages ----
SET$module<-"Microbiome"
if (is.null(SET$purpose)) SET$purpose<-"analysis"
if (is.null(SET$method_status)) SET$method_status<-"exploratory"
SET$run_id<-Run.ID
SET$output_dir<-Output.folder
validate_settings(SET)
Input.files<-validate_microbiome_config(Project.dir, SET)
need_packages(c("jsonlite", "digest", "dplyr"))
suppressPackageStartupMessages(library(dplyr))
Run<-new_run(Project.dir, SET, Input.files)
OUT<-Run$directory

# Analysis ----
tryCatch({
   pavian_obj <- read_tsv(Input.files["taxa"])
   Metric.map <- read_tsv(Input.files["column_map"])
   MET <- read_tsv(Input.files["metadata"])
   columns(MET, c("sample_id", "patient_id", "group", "include"), "Microbiome metadata")
   unique_key(MET$sample_id, "Microbiome sample_id")
   if (anyNA(MET$group) || anyNA(MET$patient_id) || any(!MET$group %in% c("Normal", "Tumor", "After")) || any(!MET$include %in% c("TRUE", "FALSE"))) abort("Explicit sample groups, patient IDs and TRUE/FALSE inclusion are required.")

   # Select the configured taxon and map quantitative metrics
   H.pylori <- extract_taxon_metrics(pavian_obj, Metric.map, MET, SET)
   H.pylori <- H.pylori[H.pylori$include == "TRUE", , drop = FALSE]
   if (!nrow(H.pylori)) abort("No included microbiome samples.")
   # Preserve supplied metrics in a separate file before any optional derivation.
   write_tsv(H.pylori, file.path(OUT, "target_imported_metrics.tsv"))
   Quantitation <- list(input_tool = SET$input_tool, database_label = SET$input_database_label,
   relative_denominator = SET$relative_denominator_description, relative_mode = SET$relative_mode, zscore_mode = SET$zscore_mode)

   # Preserve imported abundance or use an explicit count denominator
   if (SET$relative_mode == "imported") {
      if (!"relative_abundance" %in% Metric.map$metric) abort("relative_mode imported requires explicitly mapped abundance columns.")
      if (SET$relative_units == "percent") H.pylori$relative_abundance <- H.pylori$relative_abundance / 100
   } else if (SET$relative_mode == "count_over_denominator") {
      if (!all(H.pylori$sample_id %in% Metric.map$sample_id[Metric.map$metric == "count"])) abort("Count-based abundance requires counts for every included sample.")
      d <- read_tsv(Input.files["denominator"])
      columns(d, c("sample_id", "total_reads"), "Denominator")
      unique_key(d$sample_id, "Denominator sample_id")
      if (!all(H.pylori$sample_id %in% d$sample_id)) abort("Denominator does not cover all included samples.")
      values <- numeric_strict(d$total_reads[match(H.pylori$sample_id, d$sample_id)], "Denominator total_reads")
      if (any(values < 0 | values != trunc(values))) abort("Denominator must be nonnegative integer total reads.")
      if (any(H.pylori$count > values, na.rm = TRUE)) abort("Target count exceeds its explicitly supplied denominator.")
      if (any(values == 0) && SET$zero_denominator == "error") abort("Zero denominator.")
      H.pylori$total_reads <- values
      H.pylori$relative_abundance <- H.pylori$count / values
      H.pylori$relative_abundance[values == 0] <- NA_real_
   } else H.pylori$relative_abundance <- NA_real_

   # Preserve upstream Z-scores or derive a separately named Z-score
   if (SET$zscore_mode == "imported") {
      if (!any(Metric.map$metric %in% c("zscore_count", "zscore_relative"))) abort("Imported Z-score mode requires mapped upstream Z-score columns.")
      Quantitation$upstream_zscore_reference <- SET$upstream_zscore_reference
   } else if (SET$zscore_mode == "recompute") {
      Ref.ids <- reference_ids(Input.files["zscore_reference"], H.pylori$sample_id, "Z-score reference")
      Zscore <- compute_zscore(H.pylori[[SET$zscore_metric]], match(Ref.ids, H.pylori$sample_id), SET$zscore_transform, SET$zscore_na_policy, SET$zero_variance)
      # New names keep recomputed values separate from upstream Z-scores.
      H.pylori$zscore_recomputed <- Zscore$value
      Quantitation$zscore_recomputed <- list(metric = SET$zscore_metric, transform = SET$zscore_transform,
      reference_n = Zscore$reference_n, mean = Zscore$mean, sd = Zscore$sd, sd_definition = "sample standard deviation, n-1")
   } else { H.pylori$zscore_count <- NA_real_
      H.pylori$zscore_relative <- NA_real_ }

   # Optional descriptive median split
   if (SET$median_mode != "none") {
      positions <- if (SET$median_mode == "all_included") seq_len(nrow(H.pylori)) else
      match(reference_ids(Input.files["median_reference"], H.pylori$sample_id, "Median reference"), H.pylori$sample_id)
      if (all(is.na(H.pylori$relative_abundance[positions]))) abort("Median reference has no measured relative abundance.")
      Median.value <- stats::median(H.pylori$relative_abundance[positions], na.rm = TRUE)
      H.pylori$median_group <- ifelse(is.na(H.pylori$relative_abundance), if (SET$median_na_policy == "legacy_down") "DOWNmedian" else NA_character_,
      ifelse(H.pylori$relative_abundance >= Median.value, "UPmedian", "DOWNmedian"))
      Quantitation$median <- list(value_fraction = Median.value, reference_n = sum(!is.na(H.pylori$relative_abundance[positions])), na_policy = SET$median_na_policy,
      interpretation = "Descriptive source-inspired split, not a clinical HP positivity cutoff")
   }

   # Save quantitative results and provenance
   # Objects used by the original Pavian summaries
   H.pylori_reads_pct<-H.pylori %>%
      transmute(sample_id, group = factor(group, levels = c("Normal", "Tumor", "After")),
         reads_percent = relative_abundance * 100)
   if ("median_group" %in% names(H.pylori)) H.pylori_reads_pct$group2<-H.pylori$median_group
   H.pylori_Zscore_pct<-H.pylori %>%
      transmute(sample_id, group = factor(group, levels = c("Normal", "Tumor", "After")),
         Zscore_pct = zscore_relative)

   write_tsv(H.pylori, file.path(OUT, "H_pylori_metrics.tsv"))
   jsonlite::write_json(Quantitation, file.path(OUT, "quantitation_provenance.json"), pretty = TRUE, auto_unbox = TRUE)
   Details <- list(included_samples = nrow(H.pylori), final_microbiome_table_available = FALSE,
   quantitation_provenance = Quantitation, clinical_positive_negative_call = "NOT_DERIVED")
   finish_run(Run, "COMPLETED", Details)
}, error = function(e) {
   finish_run(Run, "FAILED", list(error = conditionMessage(e)))
   stop(e)
})

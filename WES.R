# WES

# Parameters ----
Project.dir<-getwd()
Install.packages<-TRUE
Run.ID<-paste0("WES-", format(Sys.time(), "%Y%m%d-%H%M%S"))
Output.folder<-"results/WES"

SET<-list(
   seed = 19970117, annovar = "data/WES/annovar.tsv", vep = "data/WES/vep.tsv",
   sample_metadata = "data/WES/samples.tsv", case_alias_map = NULL, profile = NULL,
   vep_duplicate_policy = NULL, deduplicate_policy = NULL, post_ens_exclusions = NULL,
   reference_build = "hg19", annotation_versions = NULL, make_maf = TRUE, make_oncoplot = TRUE,
   normal_artifact_rationale = NULL
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

wes_profiles <- function() {
  common <- list(function_exclude = c("intronic", "intergenic", "downstream", "upstream", "upstream;downstream",
      "UTR3", "UTR5", "ncRNA_exonic", "ncRNA_exonic;splicing", "ncRNA_intronic", "ncRNA_splicing"),
      canonical = TRUE, require_alt_gt = TRUE)
  ref_types <- c("nonsynonymous SNV", "frameshift insertion", "frameshift deletion", "stopgain", "startloss")
  list(
    wes_r_first = c(common, list(ad_ref_min = 0, ad_alt_min = 3, alt_strand_min = 1, dp_min = 20,
      sb_min = 0.125, sb_max = 8, vaf_ranges = list(c(0.1, 0.4), c(0.6, 0.9)), vaf_open = FALSE,
      coding = TRUE, consequence_column = "ExonicFunc.refGene", consequence = ref_types,
      cosmic_mode = "tissue", normal_artifact = "exclude", population_filter = FALSE)),
    wes_r_second = c(common, list(ad_ref_min = 0, ad_alt_min = 3, alt_strand_min = 0, dp_min = 50,
      sb_min = 0.25, sb_max = 4, vaf_ranges = list(c(0.10, 0.30)), vaf_open = c(FALSE, TRUE),
      coding = FALSE, consequence_column = "ExonicFunc.ensGene", consequence = c(ref_types, "stoploss"),
      cosmic_mode = "present", normal_artifact = "include_only", population_filter = TRUE)),
    wes_rmd_strict = c(common, list(ad_ref_min = 20, ad_alt_min = 20, alt_strand_min = 10, dp_min = 40,
      sb_min = 0.25, sb_max = 4, vaf_ranges = list(c(0.05, 0.45)), vaf_open = TRUE,
      coding = TRUE, consequence_column = "ExonicFunc.refGene", consequence = ref_types,
      cosmic_mode = "tissue", normal_artifact = "ignore", population_filter = FALSE)))
}

parse_genotype <- function(format, value) {
  if (is.na(value) || value %in% c(".", "./.", ".|.")) return(NULL)
  keys <- strsplit(format, ":", fixed = TRUE)[[1]]
  vals <- strsplit(value, ":", fixed = TRUE)[[1]]
  if (anyDuplicated(keys)) abort("Malformed FORMAT: duplicate tags.")
  gt_position <- match("GT", keys)
  if (!is.na(gt_position) && length(vals) >= gt_position && vals[gt_position] %in% c("./.", ".|.", ".")) return(NULL)
  if (length(keys) != length(vals)) abort("Malformed FORMAT/genotype field count.")
  names(vals) <- keys
  required(as.list(vals), c("GT", "AD", "DP", "SB"), "FORMAT")
  if (vals["GT"] %in% c("./.", ".|.", ".")) return(NULL)
  gt <- strsplit(vals["GT"], "[/|]")[[1]]
  if (any(!gt %in% c("0", "1", "."))) abort("Multiallelic GT is unsupported. Normalize biallelic records upstream.")
  ad <- strsplit(vals["AD"], ",", fixed = TRUE)[[1]]
  if (length(ad) > 2L && all(ad[-c(1, 2)] == ".")) ad <- ad[1:2]
  if (length(ad) != 2L) abort("AD must contain reference and one alternate allele only.")
  sb <- strsplit(vals["SB"], ",", fixed = TRUE)[[1]]
  if (length(sb) != 4L) abort("FORMAT SB must be ref_fw,ref_rv,alt_fw,alt_rv; INFO is not a substitute.")
  n <- numeric_strict(c(ad, vals["DP"], sb), "Genotype depths", allow_na = TRUE)
  if (any(n < 0 | n != trunc(n), na.rm = TRUE)) abort("Genotype depths must be nonnegative integers.")
  vaf <- if (anyNA(n[1:2]) || sum(n[1:2]) == 0) NA_real_ else n[2] / sum(n[1:2])
  f <- n[4] + n[6]
  r <- n[5] + n[7]
  ratio <- if (is.na(f) || is.na(r)) NA_real_ else if (r == 0) { if (f == 0) NA_real_ else Inf } else f / r
  data.frame(GT = unname(vals["GT"]), AD_ref = n[1], AD_alt = n[2], DP = n[3],
    ref_fw = n[4], ref_rv = n[5], alt_fw = n[6], alt_rv = n[7], VAF = vaf,
    SBratio = ratio, AF_reported = if ("AF" %in% keys) unname(vals["AF"]) else NA_character_,
    has_alt = "1" %in% gt, stringsAsFactors = FALSE)
}

wes_join <- function(ann, vep, sample_columns, policy) {
  columns(ann, c("ID", "Chr", "Start", "End", "Ref", "Alt", "FORMAT", "FILTER", "INFO", sample_columns), "ANNOVAR")
  columns(vep, "ID", "VEP")
  unique_key(ann$ID, "ANNOVAR ID")
  choose(policy, c("error", "retain_all"), "vep_duplicate_policy")
  if (anyNA(vep$ID) || any(!nzchar(vep$ID))) abort("VEP IDs must be present.")
  if (anyDuplicated(vep$ID) && policy == "error") abort("VEP has multiple annotations per ID; choose an explicit transcript policy.")
  overlap <- intersect(setdiff(names(ann), "ID"), setdiff(names(vep), "ID"))
  if (length(overlap)) abort("ANNOVAR/VEP ambiguous shared columns: %s. Rename explicitly upstream.", paste(overlap, collapse = ", "))
  if (any(grepl(",", ann$Alt, fixed = TRUE))) abort("Multiallelic Alt input is unsupported.")
  for (p in c("Start", "End")) {
    value <- numeric_strict(ann[[p]], paste("ANNOVAR", p))
    if (any(value <= 0 | value != trunc(value))) abort("ANNOVAR coordinates must be positive integers.")
  }
  if (any(as.numeric(ann$End) < as.numeric(ann$Start))) abort("End cannot precede Start.")
  ann$annotation_row <- seq_len(nrow(ann))
  vep$vep_row <- seq_len(nrow(vep))
  joined <- merge(ann, vep, by = "ID", all.x = TRUE, sort = FALSE)
  joined[order(joined$annotation_row, joined$vep_row, na.last = TRUE), , drop = FALSE]
}

wes_long <- function(joined, meta) {
  out <- list()
  n <- 0L
  for (j in seq_len(nrow(meta))) {
    for (i in seq_len(nrow(joined))) {
      p <- parse_genotype(joined$FORMAT[i], joined[[meta$sample_column[j]]][i])
      if (is.null(p)) next
      n <- n + 1L
      annotation <- joined[i, setdiff(names(joined), meta$sample_column), drop = FALSE]
      if (length(intersect(names(annotation), names(p)))) abort("Annotation columns collide with parsed FORMAT fields; rename them.")
      out[[n]] <- cbind(annotation, sample_id = meta$sample_id[j], role = meta$role[j], p)
    }
  }
  if (!length(out)) abort("No callable genotype records found.")
  do.call(rbind, out)
}

wes_filter <- function(x, profile) {
  p <- wes_profiles()[[profile]]
  if (is.null(p)) abort("Unknown WES filter profile.")
  columns(x, c("Func.refGene", "CANONICAL", "BIOTYPE", p$consequence_column, "cosmic96_coding", "FILTER"), "WES annotation")
  rules <- list(tumor_role = x$role == "tumor", eligible_function = !is.na(x$Func.refGene) & !x$Func.refGene %in% p$function_exclude,
    alternate_genotype = x$has_alt, AD_ref = x$AD_ref >= p$ad_ref_min,
    AD_alt = x$AD_alt >= p$ad_alt_min, DP = x$DP >= p$dp_min,
    strand_ratio = x$SBratio >= p$sb_min & x$SBratio <= p$sb_max,
    alternate_strands = x$alt_fw >= p$alt_strand_min & x$alt_rv >= p$alt_strand_min,
    canonical = x$CANONICAL == "YES", consequence = x[[p$consequence_column]] %in% p$consequence)
  ranges <- p$vaf_ranges
  vaf_ok <- rep(FALSE, nrow(x))
  open <- rep(p$vaf_open, length.out = 2)
  for (r in ranges) vaf_ok <- vaf_ok | (if (open[1]) x$VAF > r[1] else x$VAF >= r[1]) &
    (if (open[2]) x$VAF < r[2] else x$VAF <= r[2])
  rules$VAF <- vaf_ok
  if (p$coding) rules$protein_coding <- x$BIOTYPE == "protein_coding"
  rules$COSMIC <- if (p$cosmic_mode == "tissue") grepl("haematopoietic_and_lymphoid_tissue|stomach", x$cosmic96_coding) else
    !is.na(x$cosmic96_coding) & !x$cosmic96_coding %in% c(".", "", "NA")
  artifact <- grepl("normal_artifact", x$FILTER, fixed = TRUE)
  if (p$normal_artifact == "exclude") rules$normal_artifact <- !artifact & !is.na(x$FILTER)
  if (p$normal_artifact == "include_only") rules$normal_artifact <- artifact
  if (p$population_filter) {
    columns(x, "AF_eas", "Second WES chain")
    af <- numeric_strict(x$AF_eas, "AF_eas", allow_na = TRUE)
    if (any(af < 0 | af > 1, na.rm = TRUE)) abort("AF_eas must be a fraction from 0 to 1.")
    rules$population_frequency <- (!is.na(af) & af < 0.01) | (!is.na(x$AF_eas) & x$AF_eas == ".")
  }
  r <- as.data.frame(rules, check.names = FALSE)
  r[is.na(r)] <- FALSE
  reason <- apply(r, 1, function(v) paste(names(r)[!v], collapse = ";"))
  list(pass = rowSums(!r) == 0, reason = reason, rules = r)
}

wes_deduplicate <- function(x, policy) {
  choose(policy, c("error", "first_source_row"), "deduplicate_policy")
  key <- paste(x$sample_id, x$ID, sep = "::")
  if (anyDuplicated(key) && policy == "error") abort("Multiple selected transcripts per sample/variant. Confirm deduplication policy.")
  x[!duplicated(key), , drop = FALSE]
}

validate_wes_config <- function(root, cfg) {
  required(cfg, c("annovar", "vep", "sample_metadata", "profile", "vep_duplicate_policy", "deduplicate_policy",
    "post_ens_exclusions", "reference_build", "annotation_versions", "make_maf", "make_oncoplot"), "WES configuration")
  choose(cfg$profile, names(wes_profiles()), "profile")
  choose(cfg$vep_duplicate_policy, c("error", "retain_all"), "vep_duplicate_policy")
  choose(cfg$deduplicate_policy, c("error", "first_source_row"), "deduplicate_policy")
  if (!is.logical(cfg$post_ens_exclusions) || !is.logical(cfg$make_maf) || !is.logical(cfg$make_oncoplot)) abort("WES switches must be logical.")
  if (cfg$make_oncoplot && !cfg$make_maf) abort("Oncoplot requires MAF conversion.")
  if (cfg$profile == "wes_r_second") required(cfg, "normal_artifact_rationale", "This historical profile explicitly includes normal_artifact")
  paths <- c(annovar = repo_path(root, cfg$annovar), vep = repo_path(root, cfg$vep), metadata = repo_path(root, cfg$sample_metadata))
  if (!is.null(cfg$case_alias_map)) paths <- c(paths, alias_map = repo_path(root, cfg$case_alias_map))
  if (cfg$make_maf) need_packages("maftools")
  paths
}

# Packages ----
SET$module<-"WES"
if (is.null(SET$purpose)) SET$purpose<-"analysis"
if (is.null(SET$method_status)) SET$method_status<-"exploratory"
SET$run_id<-Run.ID
SET$output_dir<-Output.folder
validate_settings(SET)
Input.files<-validate_wes_config(Project.dir, SET)
need_packages(c("jsonlite", "digest", "dplyr"))
suppressPackageStartupMessages(library(dplyr))
Run<-new_run(Project.dir, SET, Input.files)
OUT<-Run$directory

# Analysis ----
tryCatch({
   ANV <- read_tsv(Input.files["annovar"])
   VEP <- read_tsv(Input.files["vep"])
   MET <- read_tsv(Input.files["metadata"])
   columns(MET, c("sample_column", "sample_id", "role", "include"), "WES sample metadata")
   unique_key(MET$sample_column, "sample_column")
   unique_key(MET$sample_id, "sample_id")
   if (any(!MET$role %in% c("tumor", "normal")) || any(!MET$include %in% c("TRUE", "FALSE"))) abort("Invalid WES role/include values.")
   MET <- MET[MET$include == "TRUE", , drop = FALSE]
   if (!nrow(MET) || !any(MET$role == "tumor")) abort("No included tumor samples.")

   # ANNOVAR / VEP
   GML <- wes_join(ANV, VEP, MET$sample_column, SET$vep_duplicate_policy)
   BEF <- wes_long(GML, MET)

   # Filtering
   Filter.result <- wes_filter(BEF, SET$profile)
   BEF$selected <- Filter.result$pass
   BEF$exclusion_reason <- Filter.result$reason
   if (SET$post_ens_exclusions) {
      columns(BEF, "ExonicFunc.ensGene", "post_ens_exclusions")
      reject <- BEF$ExonicFunc.ensGene %in% c("synonymous SNV", "nonframeshift substitution", "unknown")
      BEF$selected[reject] <- FALSE
      BEF$exclusion_reason[reject] <- paste0(BEF$exclusion_reason[reject], ";post_ens_exclusion")
   }

   # FIN / FIN2 / FIN3 / FIN4
   write_tsv(BEF, file.path(OUT, "all_genotypes_with_filter_audit.tsv"))
   keep<-rowSums(!Filter.result$rules[, names(Filter.result$rules) != "COSMIC", drop = FALSE]) == 0
   if (SET$post_ens_exclusions) keep<-keep & !reject
   FIN<-BEF[keep, , drop = FALSE]
   FIN2<-FIN
   FIN3<-FIN2 %>% filter(selected)
   FIN4<-wes_deduplicate(FIN3, SET$deduplicate_policy)
   if ("alias_map" %in% names(Input.files)) {
      aliases <- read_tsv(Input.files["alias_map"])
      columns(aliases, c("sample_id", "case_number"), "Case alias map")
      unique_key(aliases$sample_id, "Alias sample_id")
      unique_key(aliases$case_number, "case_number")
      if (!all(MET$sample_id[MET$role == "tumor"] %in% aliases$sample_id)) abort("Alias map must cover all included tumor samples, including mutation-negative samples.")
      FIN4$case_number <- aliases$case_number[match(FIN4$sample_id, aliases$sample_id)]
   }
   write_tsv(FIN4, file.path(OUT, "variants.tsv"))
   Variant.count <- data.frame(sample_id = MET$sample_id[MET$role == "tumor"], stringsAsFactors = FALSE)
   Variant.count$variant_count <- vapply(Variant.count$sample_id, function(s) sum(FIN4$sample_id == s), integer(1))
   write_tsv(Variant.count, file.path(OUT, "variant_counts.tsv"))
   jsonlite::write_json(wes_profiles()[[SET$profile]], file.path(OUT, "filter_profile_used.json"), pretty = TRUE, auto_unbox = TRUE)

   # MAF / oncoplot
   if (SET$make_maf && nrow(FIN4)) {
      columns(FIN4, c("Gene.refGene", "GeneDetail.refGene", "AAChange.refGene", "ExonicFunc.refGene"), "MAF conversion")
      tmp <- FIN4
      tmp$Tumor_Sample_Barcode <- tmp$sample_id
      # annovarToMaf creates a file-origin sample_id column. Keep the actual specimen
      # identity under a distinct name so upstream cannot duplicate/replace it.
      if ("source_sample_id" %in% names(tmp)) abort("Reserved MAF input column source_sample_id already exists.")
      tmp$source_sample_id <- tmp$sample_id
      tmp$sample_id <- NULL
      maf_input <- file.path(OUT, "annovar_for_maf.tsv")
      write_tsv(tmp, maf_input)
      # Installed upstream implementation, not the unlicensed/undated supplied copy.
      MAF <- quiet_value(maftools::annovarToMaf(annovar = maf_input, tsbCol = "Tumor_Sample_Barcode", table = "refGene",
      refBuild = SET$reference_build, MAFobj = FALSE))
      MAF <- as.data.frame(MAF)
      columns(MAF, c("ID", "Tumor_Sample_Barcode", "source_sample_id", "Chromosome", "Start_Position",
      "End_Position", "Reference_Allele", "Tumor_Seq_Allele2"), "Installed MAF converter output")
      expected_key <- paste(FIN4$sample_id, FIN4$ID, FIN4$Chr, FIN4$Start, FIN4$End, FIN4$Ref, FIN4$Alt, sep = "::")
      converted_key <- paste(MAF$Tumor_Sample_Barcode, MAF$ID, MAF$Chromosome, MAF$Start_Position, MAF$End_Position,
      MAF$Reference_Allele, MAF$Tumor_Seq_Allele2, sep = "::")
      if (nrow(MAF) != nrow(FIN4) || anyDuplicated(converted_key) || !setequal(expected_key, converted_key) ||
      !identical(as.character(MAF$source_sample_id), as.character(MAF$Tumor_Sample_Barcode)))
      abort("Installed MAF converter changed variant/sample identity. Inspect inputs and package behavior; no automatic repair applied.")
      write_tsv(MAF, file.path(OUT, "mutations.maf.tsv"))
      if (SET$make_oncoplot) {
         before <- maftools::read.maf(maf = MAF, verbose = FALSE)
         grDevices::pdf(file.path(OUT, "oncoplot_snv_indel.pdf"), width = 11, height = 8)
         tryCatch(quiet_value(maftools::oncoplot(maf = before, cohortSize = sum(MET$role == "tumor"),
         top = min(30L, length(unique(MAF$Hugo_Symbol))), removeNonMutated = FALSE,
         topBarData = Variant.count[, c("sample_id", "variant_count")], topBarLabel = "Selected variants",
         gene_mar = 8, fontSize = 0.7, titleFontSize = 1, legendFontSize = 0.85)), finally = grDevices::dev.off())
      }
   }
   Details <- list(records = nrow(FIN4), included_tumor_samples = sum(MET$role == "tumor"),
   profile = SET$profile, annotation_versions = SET$annotation_versions,
   maf = if (!SET$make_maf) "DISABLED" else if (!nrow(FIN4)) "NO_SELECTED_VARIANTS" else "COMPLETED",
   fusion = "NOT_RUN_NO_FUSION_PIPELINE_OR_INPUT")
   finish_run(Run, "COMPLETED", Details)
}, error = function(e) {
   finish_run(Run, "FAILED", list(error = conditionMessage(e)))
   stop(e)
})

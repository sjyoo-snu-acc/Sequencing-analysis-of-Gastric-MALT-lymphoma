"""Optional checks; the three analysis scripts do not depend on this file.

Only artificial input tables are used. Each execution directory contains one
configured analysis script and the required TSV inputs, without helper R files.
"""
import copy
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get("RSCRIPT_BIN") or shutil.which("Rscript")
CASES = json.loads((ROOT / "tests/synthetic/cases.json").read_text())
MODULES = ("WES", "WTS", "Microbiome")


def r_literal(value):
    if value is None:
        return "NULL"
    if value is True:
        return "TRUE"
    if value is False:
        return "FALSE"
    if isinstance(value, dict):
        if set(value) == {"__r_number__"}:
            if value["__r_number__"] not in ("Inf", "-Inf"):
                raise ValueError("Unsupported tagged R number")
            return value["__r_number__"]
        return "list(" + ", ".join(json.dumps(k) + " = " + r_literal(v) for k, v in value.items()) + ")"
    if isinstance(value, list):
        return "list(" + ", ".join(r_literal(v) for v in value) + ")"
    if isinstance(value, str):
        return json.dumps(value)
    return str(value)


def run_r_text(text, directory):
    driver = directory / "validation_driver.R"
    driver.write_text(text)
    env = os.environ.copy()
    env["MALTOMA_UNIT_REPORT"] = str(directory / "domain_unit_report.json")
    return subprocess.run([RSCRIPT, str(driver)], cwd=ROOT, capture_output=True, text=True, env=env)


def read_tsv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


@unittest.skipUnless(RSCRIPT, "Rscript is unavailable; analysis execution NOT_RUN")
class StandaloneTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from datetime import datetime, timezone
        cls.output = ROOT / "test-output" / ("standalone-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S") + "-" + str(os.getpid()))
        cls.output.mkdir(parents=True)
        result = subprocess.run([RSCRIPT, "-e", 'cat(jsonlite::toJSON(as.list(setNames(vapply(c("maftools", "DESeq2", "clusterProfiler"), function(p) suppressPackageStartupMessages(requireNamespace(p, quietly=TRUE)), logical(1)), c("maftools", "DESeq2", "clusterProfiler"))), auto_unbox=TRUE))'],
                                cwd=ROOT, capture_output=True, text=True)
        if result.returncode:
            raise unittest.SkipTest("jsonlite or the native R runtime is unavailable")
        cls.available = json.loads(result.stdout)

    @staticmethod
    def object_checks(module):
        common = '\nstopifnot(SET$purpose == "synthetic_test")\n'
        if module == "WES":
            return common + '''
stopifnot(all(vapply(c("ANV", "VEP", "GML", "BEF", "FIN", "FIN2", "FIN3", "FIN4"), exists, logical(1), envir=.GlobalEnv)))
stopifnot(nrow(FIN) >= nrow(FIN3), nrow(FIN3) >= nrow(FIN4))
stopifnot(all(FIN3$selected))
stopifnot(identical(as.character(FIN4$ID), read_tsv(file.path(OUT, "variants.tsv"))$ID))
'''
        if module == "WTS":
            return common + '''
stopifnot(all(vapply(c("CNT", "MET", "VAR", "pre", "des", "med", "obj", "res", "wrt", "cut", "UP", "DOWN"), exists, logical(1), envir=.GlobalEnv)))
stopifnot(setequal(c(UP$gene_id, DOWN$gene_id), wrt$gene_id[cut$gene_type != "ns"]))
stopifnot(length(intersect(UP$gene_id, DOWN$gene_id)) == 0L)
stopifnot(all(UP$gene_type == "up"), all(DOWN$gene_type == "down"))
'''
        return common + '''
stopifnot(all(vapply(c("pavian_obj", "H.pylori", "H.pylori_reads_pct", "H.pylori_Zscore_pct"), exists, logical(1), envir=.GlobalEnv)))
stopifnot(isTRUE(all.equal(H.pylori_reads_pct$reads_percent / 100, H.pylori$relative_abundance, check.attributes=FALSE)))
stopifnot(isTRUE(all.equal(H.pylori_Zscore_pct$Zscore_pct, H.pylori$zscore_relative, check.attributes=FALSE)))
'''

    def compare_legacy(self, module, directory, settings, out):
        # This is a private validation reference to the previous release, not a runtime dependency.
        previous = os.environ.get("MALTOMA_PREVIOUS_RELEASE")
        if not previous:
            return
        text = (Path(previous) / (module + ".R")).read_text()
        first, rest = text.split("# BEGIN USER SETTINGS", 1)
        _, last = rest.split("# END USER SETTINGS", 1)
        edit = "\nProject.dir <- getwd()\nInstall.packages <- FALSE\nData.type <- \"synthetic_test\"\n"
        edit += "Run.ID <- " + r_literal(settings["run_id"]) + "\nOutput.folder <- \"previous_results\"\nSET <- " + r_literal(settings) + "\n"
        reference = directory / "previous_reference.R"
        reference.write_text(first + "# BEGIN USER SETTINGS" + edit + "# END USER SETTINGS" + last)
        process = subprocess.run([RSCRIPT, str(reference)], cwd=directory, capture_output=True, text=True)
        self.assertEqual(process.returncode, 0, process.stderr)
        baseline = directory / "previous_results" / settings["run_id"]
        records = []
        for old_file in baseline.rglob("*.tsv"):
            new_file = out / old_file.relative_to(baseline)
            self.assertTrue(new_file.is_file())
            self.assertEqual(old_file.read_bytes(), new_file.read_bytes(), str(old_file.relative_to(baseline)))
            records.append(str(old_file.relative_to(baseline)))
        (directory / "additional_regression.json").write_text(json.dumps({"data":"SYNTHETIC_ONLY", "status":"PASS", "files":records}, indent=2))

    def configured_case(self, module, name, extra=None, label=None, edit_inputs=None, check_legacy=False):
        settings = copy.deepcopy(CASES[name])
        if extra:
            settings.update(extra)
        directory = self.output / (label or name)
        directory.mkdir()
        inputs = directory / "tests/synthetic"
        inputs.mkdir(parents=True)
        for path in (ROOT / "tests/synthetic").glob("*.tsv"):
            shutil.copyfile(path, inputs / path.name)
        if edit_inputs:
            edit_inputs(inputs)
        settings["output_dir"] = "results"
        text = (ROOT / (module + ".R")).read_text()
        begin, tail = text.split("# Parameters ----", 1)
        _, end = tail.split("# Checks and file I/O ----", 1)
        edit = "\nProject.dir <- getwd()\nInstall.packages <- FALSE\n"
        edit += "Run.ID <- " + r_literal(settings["run_id"]) + "\nOutput.folder <- \"results\"\n"
        edit += "SET <- " + r_literal(settings) + "\n"
        script = directory / (module + ".R")
        script.write_text(begin + "# Parameters ----" + edit + "# Checks and file I/O ----" + end + self.object_checks(module))
        self.assertEqual([p.name for p in directory.glob("*.R")], [module + ".R"])
        result = subprocess.run([RSCRIPT, str(script)], cwd=directory, capture_output=True, text=True)
        (directory / "stdout.txt").write_text(result.stdout)
        (directory / "stderr.txt").write_text(result.stderr)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "", "Normal analysis must not print progress/completion text")
        self.assertNotIn("completed:", result.stderr)
        self.assertNotIn("NOT_ASSESSED", result.stderr)
        out = directory / "results" / settings["run_id"]
        state = json.loads((out / "run_status.json").read_text())
        self.assertEqual(state["status"], "COMPLETED")
        self.assertEqual(state["purpose"], "synthetic_test")
        self.assertTrue((out / "configuration_used.json").is_file())
        self.assertTrue((out / "sessionInfo.txt").is_file())
        if check_legacy:
            self.compare_legacy(module, directory, settings, out)
        return directory, script, out

    def test_01_three_scripts_parse_and_resolve_their_own_functions(self):
        text = '''
for (module in c("WES", "WTS", "Microbiome")) {
   file <- paste0(module, ".R")
   parse(file)
   lines <- readLines(file)
   first <- which(lines == "# Checks and file I/O ----") + 1L
   last <- which(lines == "# Packages ----") - 1L
   env <- new.env(parent = globalenv())
   eval(parse(text = lines[first:last]), envir = env)
   for (name in ls(env, all.names = TRUE)) {
      fun <- get(name, envir = env)
      if (!is.function(fun)) next
      for (global in codetools::findGlobals(fun, merge = FALSE)$functions) {
         if (!exists(global, envir = env, mode = "function", inherits = TRUE)) stop(paste(module, name, global))
      }
   }
}
'''
        result = run_r_text(text, self.output)
        self.assertEqual(result.returncode, 0, result.stderr)
        files = [str(p.relative_to(ROOT)) for p in ROOT.rglob("*.R") if "test-output" not in p.parts and "results" not in p.parts]
        self.assertEqual(set(files), {"WES.R", "WTS.R", "Microbiome.R"})
        for module in MODULES:
            text = (ROOT / (module + ".R")).read_text()
            self.assertIsNone(re.search(r"\b(?:source|sys\.source|message|cat|print)\s*\(", text))
            self.assertNotIn("yaml::", text)
            self.assertNotIn("maltoma.validation", text)

    def test_02_domain_units(self):
        text = DOMAIN_UNIT_R
        result = run_r_text(text, self.output)
        (self.output / "domain_units.stdout.txt").write_text(result.stdout)
        (self.output / "domain_units.stderr.txt").write_text(result.stderr)
        self.assertEqual(result.returncode, 0, result.stderr)
        report = json.loads((self.output / "domain_unit_report.json").read_text())
        self.assertEqual(report["fail"], 0)
        self.assertGreaterEqual(report["pass"], 42)

    def test_03_wes(self):
        directory, script, out = self.configured_case("WES", "synthetic_wes")
        variants = read_tsv(out / "variants.tsv")
        self.assertEqual([x["ID"] for x in variants], ["SYN_VAR_1", "SYN_VAR_3"])
        hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob("*.tsv")}
        repeat = subprocess.run([RSCRIPT, str(script)], cwd=directory, capture_output=True, text=True)
        self.assertNotEqual(repeat.returncode, 0)
        self.assertIn("already exists", repeat.stderr)
        self.assertEqual(hashes, {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.glob("*.tsv")})

    def test_04_wes_maf_oncoplot(self):
        if not self.available["maftools"]:
            self.skipTest("maftools unavailable; MAF/oncoplot NOT_RUN")
        _, _, out = self.configured_case("WES", "synthetic_wes_maf")
        records = read_tsv(out / "mutations.maf.tsv")
        self.assertEqual([x["ID"] for x in records], ["SYN_VAR_1", "SYN_VAR_3"])
        self.assertTrue(all(x["Tumor_Sample_Barcode"] == x["source_sample_id"] == "SYN_WES_A" for x in records))
        self.assertGreater((out / "oncoplot_snv_indel.pdf").stat().st_size, 0)

    def test_05_wts_four_contrasts(self):
        if not self.available["DESeq2"]:
            self.skipTest("DESeq2 unavailable; WTS NOT_RUN")
        _, _, out = self.configured_case("WTS", "synthetic_wts")
        records = read_tsv(out / "S1_S5_S9/DEG_all.tsv")
        signal = [float(x["log2FoldChange"]) for x in records if x["gene_id"] in {"SYN_FEATURE_" + str(i) for i in range(1, 31)}]
        self.assertGreater(sorted(signal)[len(signal) // 2], 1)
        self.assertEqual(len(list(out.glob("*/DEG_all.tsv"))), 4)

    def test_06_wts_paired(self):
        if not self.available["DESeq2"]:
            self.skipTest("DESeq2 unavailable; paired WTS NOT_RUN")
        _, _, out = self.configured_case("WTS", "synthetic_wts_paired")
        state = json.loads((out / "S1_S5_S9/contrast_used.json").read_text())
        self.assertEqual(state["design"], "~ patient_id + group")
        self.assertEqual(state["model_kind"], "paired")

    def test_07_microbiome_imports(self):
        _, _, out = self.configured_case("Microbiome", "synthetic_microbiome")
        rows = read_tsv(out / "H_pylori_metrics.tsv")
        self.assertEqual([float(x["relative_abundance"]) for x in rows[:3]], [0.01, 0, 0.03])
        self.assertEqual(rows[3]["count"], "NA")
        self.assertEqual(rows[3]["median_group"], "NA")

    def test_08_microbiome_derivation(self):
        _, _, out = self.configured_case("Microbiome", "synthetic_microbiome_recomputed")
        rows = read_tsv(out / "H_pylori_metrics.tsv")
        values = [0.01, 0, 0.03]
        mean = sum(values) / len(values)
        sd = math.sqrt(sum((x - mean) ** 2 for x in values) / 2)
        for row, value in zip(rows, values):
            self.assertAlmostEqual(float(row["zscore_recomputed"]), (value - mean) / sd)
        self.assertEqual([float(x["zscore_relative"]) for x in rows[:3]], [-0.2, -0.8, 1.4])

    def test_09_gsea_statistical_engine(self):
        if not self.available["clusterProfiler"]:
            self.skipTest("clusterProfiler unavailable; GO/KEGG engine NOT_RUN; ranking adapters have separate unit checks")
        _, _, out = self.configured_case("WTS", "synthetic_wts_local_gsea")
        for kind in ("GO", "KEGG"):
            rows = read_tsv(out / "S1_S5_S9" / (kind + "_GSEA.tsv"))
            self.assertTrue(rows)
            self.assertTrue(all(x["ID"].startswith("SYN_" + kind + "_") for x in rows))

    def test_10_unresolved_study_inputs_stop(self):
        for module in MODULES:
            with self.subTest(module=module), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                text = re.sub(r"Install\.packages\s*<-\s*TRUE", "Install.packages <- FALSE", (ROOT / (module + ".R")).read_text(), count=1)
                script = directory / (module + ".R")
                script.write_text(text)
                result = subprocess.run([RSCRIPT, str(script)], cwd=directory, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn("unresolved/missing fields", result.stderr)
                self.assertFalse((directory / "results").exists())

    def test_11_quiet_value_keeps_values_errors_and_warnings(self):
        text = '''
lines <- readLines("Microbiome.R")
eval(parse(text = lines[(which(lines == "# Checks and file I/O ----")+1):(which(lines == "# Packages ----")-1)]))
Install.packages <- FALSE
value <- quiet_value({message("artificial progress"); cat("artificial output"); 42})
stopifnot(identical(value, 42))
warning_seen <- FALSE
withCallingHandlers(quiet_value({warning("artificial warning"); 1}), warning = function(w) {warning_seen <<- TRUE; invokeRestart("muffleWarning")})
stopifnot(warning_seen)
error <- tryCatch(quiet_value(stop("artificial error")), error = conditionMessage)
stopifnot(error == "artificial error")
stopifnot(identical(is_bioc_package(c("org.Hs.eg.db", "DESeq2", "jsonlite")), c(TRUE, TRUE, FALSE)))
need_packages(c("base", "stats"))
missing <- tryCatch(need_packages("MALTOMA_SYNTHETIC_MISSING_PACKAGE"), error = conditionMessage)
stopifnot(grepl("Missing R packages", missing))
'''
        result = run_r_text(text, self.output)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotIn("artificial progress", result.stderr)
        self.assertNotIn("artificial output", result.stdout)

    def test_12_wes_other_profiles_and_post_filter(self):
        for profile, expected in (("wes_r_second", ["SYN_VAR_2"]), ("wes_rmd_strict", ["SYN_VAR_1", "SYN_VAR_2", "SYN_VAR_3"])):
            with self.subTest(profile=profile):
                _, _, out = self.configured_case("WES", "synthetic_wes", extra={"profile":profile, "normal_artifact_rationale":"Artificial branch test only"},
                    label=profile, check_legacy=True)
                self.assertEqual([x["ID"] for x in read_tsv(out / "variants.tsv")], expected)

        def synonymous(inputs):
            file=inputs/"wes_annovar.tsv"
            with file.open(newline="") as stream:
                reader=csv.DictReader(stream, delimiter="\t"); fields=reader.fieldnames; rows=list(reader)
            rows[0]["ExonicFunc.ensGene"]="synonymous SNV"
            with file.open("w", newline="") as stream:
                writer=csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n");writer.writeheader();writer.writerows(rows)

        _, _, out=self.configured_case("WES", "synthetic_wes", extra={"post_ens_exclusions":True},
            label="post_ens_exclusion", edit_inputs=synonymous, check_legacy=True)
        self.assertEqual([x["ID"] for x in read_tsv(out/"variants.tsv")], ["SYN_VAR_3"])

    def test_13_wes_zero_selected_variants(self):
        def absent_cosmic(inputs):
            file=inputs/"wes_annovar.tsv"
            with file.open(newline="") as stream:
                reader=csv.DictReader(stream,delimiter="\t");fields=reader.fieldnames;rows=list(reader)
            for row in rows:row["cosmic96_coding"]="."
            with file.open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=fields,delimiter="\t",lineterminator="\n");writer.writeheader();writer.writerows(rows)
        _, _, out=self.configured_case("WES", "synthetic_wes", label="zero_variants", edit_inputs=absent_cosmic, check_legacy=True)
        self.assertEqual(read_tsv(out/"variants.tsv"), [])
        self.assertTrue(all(row["variant_count"]=="0" for row in read_tsv(out/"variant_counts.tsv")))

    def test_14_wts_all_included_fit(self):
        if not self.available["DESeq2"]:
            self.skipTest("DESeq2 unavailable; all_included branch NOT_RUN")
        comparisons=copy.deepcopy(CASES["synthetic_wts"]["comparisons"])
        for comparison in comparisons:comparison["reference_levels"]={"group":"HPneg_Normal"}
        _, _, out=self.configured_case("WTS", "synthetic_wts", extra={"fit_scope":"all_included", "comparisons":comparisons},
            label="all_included", check_legacy=True)
        self.assertEqual(len(list(out.glob("*/DEG_all.tsv"))), 4)
        for metadata in out.glob("*/metadata_used.tsv"):
            self.assertEqual(len(read_tsv(metadata)), 16)


# Generated below from the established artificial-input domain assertions.
DOMAIN_UNIT_R = 'root <- getwd()\nInstall.packages <- FALSE\nfor (module in c("WES", "WTS", "Microbiome")) {\n   lines <- readLines(paste0(module, ".R"))\n   first <- which(lines == "# Checks and file I/O ----") + 1L\n   last <- which(lines == "# Packages ----") - 1L\n   eval(parse(text = lines[first:last]))\n}\ndecode_numbers <- function(x) {\n   if (is.list(x) && identical(names(x), "__r_number__")) return(as.numeric(x[[1]]))\n   if (is.list(x)) return(lapply(x, decode_numbers))\n   x\n}\ncases <- decode_numbers(jsonlite::read_json("tests/synthetic/cases.json", simplifyVector = FALSE))\nread_case <- function(name) cases[[name]]\nneed_packages(c("digest", "jsonlite"))\nreport <- list(); report_n <- 0L\ntest <- function(name, expression) {\n  report_n <<- report_n + 1L\n  tryCatch({ force(expression); report[[report_n]] <<- data.frame(test = name, status = "PASS", detail = "") },\n    error = function(e) report[[report_n]] <<- data.frame(test = name, status = "FAIL", detail = conditionMessage(e)))\n  message(report[[report_n]]$status, " ", name)\n}\nexpect_error <- function(expression, pattern) {\n  msg <- tryCatch({ force(expression); NULL }, error = function(e) conditionMessage(e))\n  if (is.null(msg) || !grepl(pattern, msg, ignore.case = TRUE)) abort("Expected error matching \'%s\'; received \'%s\'.", pattern, msg %||% "no error")\n}\nequal <- function(x, y, tolerance = 1e-10) if (!isTRUE(all.equal(x, y, tolerance = tolerance, check.attributes = FALSE))) abort("Unexpected value in independent assertion.")\nbase <- file.path(root, "tests", "synthetic")\nif (!dir.exists(base)) abort("Run Restore the supplied tests/synthetic fixtures.")\nann <- read_tsv(file.path(base, "wes_annovar.tsv")); vep <- read_tsv(file.path(base, "wes_vep.tsv"))\nmeta_wes <- read_tsv(file.path(base, "wes_metadata.tsv"))\ntest("FORMAT parsing by named keys", equal(parse_genotype("DP:SB:AD:GT", "100:40,40,10,10:80,20:0/1")$VAF, 0.2))\ntest("Trailing placeholder AD is explicit", equal(parse_genotype("GT:AD:DP:SB", "0/1:80,20,.:100:40,40,10,10")$AD_alt, 20))\ntest("No-call genotype is not a variant", stopifnot(is.null(parse_genotype("GT:AD:DP:SB", "./.:.:.:."))))\ntest("Malformed FORMAT stops", expect_error(parse_genotype("GT:AD:DP:SB", "0/1:80,20:100"), "FORMAT"))\ntest("Multiallelic GT stops", expect_error(parse_genotype("GT:AD:DP:SB", "0/2:80,20:100:40,40,10,10"), "Multiallelic"))\ntest("Negative read depth stops", expect_error(parse_genotype("GT:AD:DP:SB", "0/1:80,-1:100:40,40,10,10"), "nonnegative"))\ntest("Zero VAF denominator remains NA", stopifnot(is.na(parse_genotype("GT:AD:DP:SB", "0/1:0,0:0:0,0,0,0")$VAF)))\ntest("INFO cannot replace missing FORMAT SB", expect_error(parse_genotype("GT:AD:DP", "0/1:80,20:100"), "missing fields"))\njoined <- wes_join(ann, vep, meta_wes$sample_column, "error"); long <- wes_long(joined, meta_wes)\ntest("First WES chain keeps two expected artificial variants", {\n  keep <- wes_filter(long, "wes_r_first")$pass\n  equal(long$ID[keep], c("SYN_VAR_1", "SYN_VAR_3")); stopifnot(all(long$role[keep] == "tumor"))\n})\ntest("Second WES chain explicitly includes normal artifact", equal(long$ID[wes_filter(long, "wes_r_second")$pass], "SYN_VAR_2"))\ntest("Strict WES chain remains a distinct method", equal(sum(wes_filter(long, "wes_rmd_strict")$pass), 3))\ntest("Missing/zero metrics do not become accepted variants", { x <- long; x$VAF <- NA_real_; stopifnot(!any(wes_filter(x, "wes_r_first")$pass)) })\ntest("Duplicate VEP default stops", expect_error(wes_join(ann, rbind(vep, vep[1, ]), meta_wes$sample_column, "error"), "multiple annotations"))\ntest("Named join collision stops", { v <- vep; v$Chr <- "chr1"; expect_error(wes_join(ann, v, meta_wes$sample_column, "error"), "shared columns") })\ntest("Variant dedup never creates Cartesian NA rows", {\n  x <- rbind(long[1, ], long[1, ]); stopifnot(nrow(wes_deduplicate(x, "first_source_row")) == 1L)\n  expect_error(wes_deduplicate(x, "error"), "Multiple selected")\n})\nmeta_wts <- read_tsv(file.path(base, "wts_metadata.tsv"))\ncounts <- read_counts(file.path(base, "wts_counts.tsv"), meta_wts, "transcript_id", "preserve", "error")\nscratch <- tempfile(tmpdir = tempdir()); dir.create(scratch)\non.exit_cleanup <- function() unlink(scratch, recursive = TRUE)\ntest("WTS counts preserve raw integer scale and sample order", { stopifnot(is.integer(counts)); equal(colnames(counts), meta_wts$sample_id) })\ntest("Fractional counts rejected without rounding", {\n  raw <- read_tsv(file.path(base, "wts_counts.tsv")); raw[1, 3] <- "1.5"; path <- file.path(scratch, "fractional.tsv"); write_tsv(raw, path)\n  expect_error(read_counts(path, meta_wts, "transcript_id", "preserve", "error"), "raw integer")\n})\ntest("Missing count values rejected", {\n  raw <- read_tsv(file.path(base, "wts_counts.tsv")); raw[1, 3] <- "NA"; path <- file.path(scratch, "missing.tsv"); write_tsv(raw, path)\n  expect_error(read_counts(path, meta_wts, "transcript_id", "preserve", "error"), "missing numbers")\n})\ntest("Sample mismatch rejected", { m <- meta_wts[-1, ]; expect_error(read_counts(file.path(base, "wts_counts.tsv"), m, "transcript_id", "preserve", "error"), "match exactly") })\ntest("HP metadata conflict rejected", { m <- meta_wts; m$type[1] <- "After"; expect_error(validate_wts_metadata(m), "conflicts") })\npaired <- read_tsv(file.path(base, "wts_paired_metadata.tsv"))\npc <- read_case("synthetic_wts_paired")$comparisons[[1]]\ntest("Complete artificial pairs yield full-rank model", { x <- prepare_model(paired, pc); stopifnot(qr(x$matrix)$rank == ncol(x$matrix)) })\ntest("Incomplete pairing stops", expect_error(prepare_model(paired[-1, ], pc), "exactly one"))\ntest("Independent model cannot silently treat repeated patients", {\n  c <- pc; c$design <- "~ group"; c$model_kind <- "independent"; expect_error(prepare_model(paired, c), "Repeated patients")\n})\ntest("Confounded patient and between-patient group design stops", {\n  m <- paired; m$patient_id <- paste0("SYN_SINGLE_", seq_len(nrow(m))); expect_error(prepare_model(m, pc), "exactly one")\n})\ntest("Unresolved design is an explicit error", { c <- pc; c$design <- NULL; expect_error(prepare_model(paired, c), "unresolved") })\nranking_cfg <- list(metric = "log2FoldChange", prefilter = "all_tested", duplicate_rank_ids = "error", tie_policy = "stable_source_order")\nresult <- data.frame(gene_id = c("SYN_A", "SYN_B", "SYN_C"), log2FoldChange = c(2, -3, 1), stat = c(1, -1, 0.5), pvalue = c(0.01, 0.02, 0.2))\ntest("GSEA ranking order is metric-based with no significance prefilter", equal(names(make_ranking(result, ranking_cfg)), c("SYN_A", "SYN_C", "SYN_B")))\ntest("Source-alternative pvalue prefilter is separate", {\n  c <- ranking_cfg; c$prefilter <- "pvalue_lt"; c$pvalue_cutoff <- 0.05; equal(names(make_ranking(result, c)), c("SYN_A", "SYN_B"))\n})\ntest("KEGG mapping is keyed by gene ID, never row position", {\n  mapping <- data.frame(gene_id = c("SYN_C", "SYN_A", "SYN_B"), entrez_id = c("E3", "E1", "E2"))\n  equal(names(make_ranking(result, ranking_cfg, mapping, "entrez_id")), c("E1", "E3", "E2"))\n})\ntest("Duplicate KEGG rank IDs stop by default", {\n  mapping <- data.frame(gene_id = result$gene_id, entrez_id = c("E1", "E1", "E2"))\n  expect_error(make_ranking(result, ranking_cfg, mapping, "entrez_id"), "Duplicate")\n})\nmicro_cfg <- read_case("synthetic_microbiome")\ntaxa <- read_tsv(file.path(base, "microbiome_taxa.tsv")); map <- read_tsv(file.path(base, "microbiome_column_map.tsv")); mm <- read_tsv(file.path(base, "microbiome_metadata.tsv"))\ntest("Target identity does not depend on first row", equal(select_taxon(taxa, micro_cfg), 2L))\ntest("Count zero and missing remain distinct", { x <- extract_taxon_metrics(taxa, map, mm, micro_cfg); equal(x$count[2], 0); stopifnot(is.na(x$count[4])) })\ntest("Taxon ID/name disagreement stops", { t <- taxa; t$taxon_name[2] <- "Synthetic wrong name"; expect_error(select_taxon(t, micro_cfg), "disagree") })\ntest("Duplicate target taxon stops", expect_error(select_taxon(rbind(taxa, taxa[2, ]), micro_cfg), "Duplicate"))\ntest("Absent target not assumed to be zero", expect_error(select_taxon(taxa[1, , drop = FALSE], micro_cfg), "absent"))\ntest("Imported fractions cannot be summed across duplicate taxa", { c <- micro_cfg; c$duplicate_taxa <- "sum_counts_only"; expect_error(extract_taxon_metrics(rbind(taxa, taxa[2, ]), map, mm, c), "only for counts") })\ntest("Incorrect abundance units rejected", { t <- taxa; t$SYN_A_pct[2] <- "101"; expect_error(extract_taxon_metrics(t, map, mm, micro_cfg), "outside") })\ntest("Z-score mean/sd checked independently", { z <- compute_zscore(c(1, 2, 3), 1:3, "identity", "error", "error"); equal(z$value, c(-1, 0, 1)); equal(z$sd, 1) })\ntest("Zero reference variance error is explicit", expect_error(compute_zscore(c(1, 1, 1), 1:3, "identity", "error", "error"), "zero variance"))\ntest("Zero variance NA policy does not invent zeros", stopifnot(all(is.na(compute_zscore(c(1, 1, 1), 1:3, "identity", "error", "NA")$value))))\ntest("Missing reference values are not imputed", expect_error(compute_zscore(c(1, NA, 3), 1:3, "identity", "error", "error"), "Missing"))\ntest("Absolute paths are rejected", expect_error(repo_path(root, "/absolute/input.tsv"), "relative"))\ntest("Directory escape rejected", expect_error(repo_path(root, "../input.tsv"), "escape"))\n\n\nresult <- do.call(rbind, report)\njsonlite::write_json(list(data = "SYNTHETIC_ONLY", pass = sum(result$status == "PASS"), fail = sum(result$status == "FAIL"), results = result), Sys.getenv("MALTOMA_UNIT_REPORT"), pretty = TRUE, auto_unbox = TRUE)\non.exit_cleanup()\nif (any(result$status == "FAIL")) quit(status = 1L)\n'


if __name__ == "__main__":
    unittest.main(verbosity=2)

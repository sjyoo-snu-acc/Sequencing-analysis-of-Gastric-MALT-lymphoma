"""Artificial records and temporary test workbooks only; no study inputs used."""
import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

path = Path(__file__).parents[1] / "tools/compare_to_excel.py"
spec = importlib.util.spec_from_file_location("comparison", path)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def write_tsv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def write_contrast(directory, reverse=False):
    state = {"factor": "group", "numerator": "SYN_POS", "denominator": "SYN_NEG"}
    if reverse:
        state["numerator"], state["denominator"] = state["denominator"], state["numerator"]
    (directory / "contrast_used.json").write_text(json.dumps(state))
    write_tsv(directory / "metadata_used.tsv", [
        {"sample_id": "SYN_A", "group": "SYN_POS", "type": "Normal", "hp_status": "Positive"},
        {"sample_id": "SYN_B", "group": "SYN_NEG", "type": "Normal", "hp_status": "Negative"},
    ])


class ComparisonTests(unittest.TestCase):
    def test_match_by_id_not_row_order(self):
        a = [{"ID": "SYN_A", "NES": 1.2}, {"ID": "SYN_B", "NES": -2.0}]
        self.assertEqual(module.compare_records(a, list(reversed(a)), ["ID"], ["NES"])["differences"], [])

    def test_differences_are_not_repaired(self):
        a = [{"ID": "SYN_A", "pvalue": 0.8, "padj": 0.01}]
        b = [{"ID": "SYN_A", "pvalue": 0.01, "padj": 0.8}]
        report = module.compare_records(a, b, ["ID"], ["pvalue", "padj"], atol=0, rtol=0)
        self.assertEqual(len(report["differences"]), 2)
        self.assertEqual(a[0]["pvalue"], 0.8)

    def test_duplicate_identity_is_ambiguous(self):
        a = [{"ID": "SYN_A", "NES": 1}]*2
        with self.assertRaisesRegex(ValueError, "Ambiguous duplicate"):
            module.compare_records(a, a, ["ID"], ["NES"])

    def test_missing_and_extra_identities(self):
        a = [{"ID": "SYN_A", "NES": 1}]; b = [{"ID": "SYN_B", "NES": 1}]
        self.assertEqual({d["status"] for d in module.compare_records(a, b, ["ID"], ["NES"])["differences"]}, {"EXTRA_COMPUTED", "MISSING_COMPUTED"})

    def test_missing_is_not_numeric_zero(self):
        a = [{"ID": "SYN_A", "value": None}]; b = [{"ID": "SYN_A", "value": 0}]
        self.assertEqual(len(module.compare_records(a, b, ["ID"], ["value"])["differences"]), 1)

    def test_case_alias_is_required(self):
        with self.assertRaisesRegex(ValueError, "Missing comparison identity"):
            module.compare_records([{"case": "NA"}], [{"case": "NA"}], ["case"], [])

    def test_contrast_direction_is_checked(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_contrast(directory, reverse=True)
            with self.assertRaisesRegex(ValueError, "direction"):
                module.check_wts_direction({"expected_numerator": "SYN_POS", "expected_denominator": "SYN_NEG"}, directory, 1)

    def test_group_semantics_cannot_be_inferred_from_labels(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            write_contrast(directory)
            rows = module.read_tsv(directory / "metadata_used.tsv")
            rows[0]["type"] = "Tumor"
            write_tsv(directory / "metadata_used.tsv", rows)
            with self.assertRaisesRegex(ValueError, "metadata conflicts"):
                module.check_wts_direction({"expected_numerator": "SYN_POS", "expected_denominator": "SYN_NEG"}, directory, 1)

    @unittest.skipUnless(importlib.util.find_spec("openpyxl"), "Optional openpyxl unavailable")
    def test_artificial_workbook_cli_all_modules_and_no_overwrite(self):
        import openpyxl
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            workbook = openpyxl.Workbook()
            sheet = workbook.active
            sheet.title = "Table S13"
            wes = {"case_number": "SYN_CASE_A", "gene_symbol": "SYN_GENE_A", "variant_function": "SYN_MISSENSE",
                   "protein_position": "1", "amino_acids": "A/G", "cosmic96_coding": "SYNTHETIC_ONLY", "INFO": "SYNTHETIC_ONLY",
                   "Chr": "SYN_CHR", "Start": 101, "End": 101, "Ref": "A", "Alt": "G"}
            sheet.append(list(wes)); sheet.append(list(wes.values()))
            wts = {"GENE": "SYN_GENE_A", "baseMean": 120, "log2FoldChange": 1.2, "pvalue": 0.01, "padj": 0.02}
            sheet = workbook.create_sheet("Table S1")
            sheet.append(list(wts)); sheet.append(list(wts.values()))
            sheet.append(["Table S1. Artificial test caption only"])
            workbook_path = directory / "synthetic_only.xlsx"
            workbook.save(workbook_path)
            write_tsv(directory / "variants.tsv", [wes])
            wts_computed = dict(wts); wts_computed["pvalue"] = 0.7
            write_tsv(directory / "DEG_selected.tsv", [wts_computed])
            write_contrast(directory)
            micro = {"taxon_id": "SYN_TAXON", "taxon_name": "Synthetic species", "sample_id": "SYN_A", "count": 10, "relative_abundance": 0.1}
            write_tsv(directory / "micro_reported.tsv", [micro]); write_tsv(directory / "micro_computed.tsv", [micro])
            manifest = {
                "wes": {"variants": str(directory / "variants.tsv"), "columns": {key: key for key in wes}},
                "wts": [{"table_number": 1, "results": str(directory / "DEG_selected.tsv"), "expected_numerator": "SYN_POS",
                         "expected_denominator": "SYN_NEG", "absolute_tolerance": 0, "relative_tolerance": 0}],
                "microbiome": {"reported_tsv": str(directory / "micro_reported.tsv"), "computed_tsv": str(directory / "micro_computed.tsv"),
                               "metrics": ["count", "relative_abundance"], "absolute_tolerance": 0, "relative_tolerance": 0},
            }
            manifest_path = directory / "manifest.json"; manifest_path.write_text(json.dumps(manifest))
            output = directory / "comparison"
            command = [sys.executable, str(path), "--supplement", str(workbook_path), "--manifest", str(manifest_path), "--output", str(output)]
            completed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            report = json.loads((output / "comparison.json").read_text())
            self.assertFalse(report["manuscript_reproduced"])
            self.assertEqual(report["modules"]["WES"]["differences"], [])
            self.assertEqual(report["modules"]["Microbiome"]["differences"], [])
            self.assertEqual(report["modules"]["WTS_S1"]["differences"][0]["field"], "pvalue")
            self.assertEqual(len(report["modules"]["WTS_S1"]["differences"]), 1)
            previous = (output / "comparison.json").read_bytes()
            repeated = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(repeated.returncode, 0)
            self.assertEqual(previous, (output / "comparison.json").read_bytes())


if __name__ == "__main__":
    unittest.main()

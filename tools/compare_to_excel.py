#!/usr/bin/env python3
"""Read-only identity/numeric comparison. Requires user-supplied private workbooks.

No filter, pvalue, adjusted value, alias or comparison direction is rewritten.
The output directory must not already exist. Results are private by default.
"""
import argparse
import csv
import json
import math
from pathlib import Path


def read_tsv(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if len(reader.fieldnames or []) != len(set(reader.fieldnames or [])):
            raise ValueError("Duplicate TSV column names")
        return list(reader)


def key_index(rows, keys, context):
    result = {}
    for row in rows:
        key = tuple(str(row[k]) for k in keys)
        if any(x in ("", "None", "NA") for x in key):
            raise ValueError(f"Missing comparison identity in {context}")
        if key in result:
            raise ValueError(f"Ambiguous duplicate comparison identity in {context}; no first-match selection applied")
        result[key] = row
    return result


def number(x):
    if x is None or str(x) in ("", "NA", "NaN", "nan"):
        return None
    value = float(x)
    if not math.isfinite(value):
        raise ValueError("Nonfinite numeric comparison value")
    return value


def compare_records(expected, observed, keys, numeric_fields, text_fields=(), atol=1e-8, rtol=1e-6):
    e = key_index(expected, keys, "reported table")
    o = key_index(observed, keys, "computed output")
    differences = []
    for key in sorted(set(e) | set(o)):
        if key not in e or key not in o:
            differences.append({"identity": key, "field": "record", "status": "EXTRA_COMPUTED" if key not in e else "MISSING_COMPUTED",
                                "reported": None, "computed": None})
            continue
        for field in numeric_fields:
            a = number(e[key][field]); b = number(o[key][field])
            same = (a is None and b is None) or (a is not None and b is not None and math.isclose(a, b, abs_tol=atol, rel_tol=rtol))
            if not same:
                differences.append({"identity": key, "field": field, "status": "VALUE_DIFF", "reported": a, "computed": b})
        for field in text_fields:
            a = str(e[key][field]); b = str(o[key][field])
            if a != b:
                differences.append({"identity": key, "field": field, "status": "VALUE_DIFF", "reported": a, "computed": b})
    return {"reported_records": len(e), "computed_records": len(o), "matched_identities": len(set(e) & set(o)),
            "absolute_tolerance": atol, "relative_tolerance": rtol, "differences": differences}


def worksheet_records(workbook, name):
    sheet = workbook[name]
    values = list(sheet.values)
    headers = [str(v) if v is not None else "" for v in values[0]]
    if len(headers) != len(set(headers)):
        raise ValueError(f"Duplicate Excel headers: {name}")
    # Final captions are prose rows, not analytical records.
    return [dict(zip(headers, row)) for row in values[1:] if row[0] is not None and not str(row[0]).startswith("Table S")]


def check_wts_direction(part, result_dir, table_number):
    state = json.loads((result_dir / "contrast_used.json").read_text())
    expected_numerator = part["expected_numerator"]
    expected_denominator = part["expected_denominator"]
    if state["numerator"] != expected_numerator or state["denominator"] != expected_denominator:
        raise ValueError("Computed contrast direction differs from the explicit comparison manifest")
    # Verify actual metadata semantics, rather than assuming the meaning of a group name.
    metadata = read_tsv(result_dir / "metadata_used.tsv")
    roles = {
        1: (("Normal", "Positive"), ("Normal", "Negative")),
        2: (("Tumor", "Positive"), ("Normal", "Positive")),
        3: (("Tumor", "Positive"), ("Tumor", "Negative")),
        4: (("Tumor", "Negative"), ("Normal", "Negative")),
    }
    factor = state["factor"]
    for label, semantics in zip((expected_numerator, expected_denominator), roles[(table_number - 1) % 4 + 1]):
        selected = [r for r in metadata if r[factor] == label]
        if not selected or any((r["type"], r["hp_status"]) != semantics for r in selected):
            raise ValueError("Computed group metadata conflicts with the final sheet comparison caption")
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--supplement", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    import openpyxl  # Read-only dependency; never used to create or edit a workbook.
    manifest = json.loads(args.manifest.read_text())
    workbook = openpyxl.load_workbook(args.supplement, data_only=True, read_only=True)
    if args.output.exists():
        raise ValueError("Refusing to overwrite comparison output directory")
    args.output.mkdir(parents=True)
    report = {"status": "COMPARISON_OF_SUPPLIED_OUTPUTS_ONLY", "manuscript_reproduced": False, "modules": {}}
    try:
        if manifest.get("wes"):
            part = manifest["wes"]
            observed = read_tsv(part["variants"])
            mapping = part["columns"]
            expected = worksheet_records(workbook, "Table S13")
            transformed = [{k: row[v] for k, v in mapping.items()} for row in observed]
            # Case aliases are required. No matching of aliases from numeric suffixes.
            report["modules"]["WES"] = compare_records(expected, transformed,
                ["case_number", "Chr", "Start", "End", "Ref", "Alt"], [],
                ["gene_symbol", "variant_function", "protein_position", "amino_acids", "cosmic96_coding", "INFO"], 0, 0)
        for part in manifest.get("wts", []):
            n = int(part["table_number"])
            if not 1 <= n <= 12:
                raise ValueError("WTS table_number must be 1-12")
            expected = worksheet_records(workbook, f"Table S{n}")
            computed = read_tsv(part["results"])
            result_dir = Path(part["results"]).parent
            direction = check_wts_direction(part, result_dir, n)
            identity = "GENE" if n <= 4 else "ID"
            fields = ["baseMean", "log2FoldChange", "pvalue", "padj"] if n <= 4 else ["NES", "pvalue", "p.adjust"]
            columns = part.get("columns", {k: k for k in [identity] + fields})
            observed = [{k: row[v] for k, v in columns.items()} for row in computed]
            comparison = compare_records(expected, observed, [identity], fields,
                atol=float(part["absolute_tolerance"]), rtol=float(part["relative_tolerance"]))
            comparison["contrast_verified"] = direction
            if n <= 4:
                comparison["identifier_limitation"] = "Final sheets provide symbols, not Ensembl feature IDs. A unique verified symbol mapping is required; ambiguous symbols cause an error."
            report["modules"][f"WTS_S{n}"] = comparison
        if manifest.get("microbiome"):
            part = manifest["microbiome"]
            # No microbiome Excel sheet was supplied. A separately authorized reported TSV is required.
            expected = read_tsv(part["reported_tsv"]); computed = read_tsv(part["computed_tsv"])
            report["modules"]["Microbiome"] = compare_records(expected, computed,
                ["taxon_id", "sample_id"], part["metrics"], ["taxon_name"],
                atol=float(part["absolute_tolerance"]), rtol=float(part["relative_tolerance"]))
        else:
            report["modules"]["Microbiome"] = {"status": "NOT_RUN", "reason": "No corresponding final microbiome table provided"}
        if not manifest.get("wes") and not manifest.get("wts") and not manifest.get("microbiome"):
            raise ValueError("No computed results supplied. The comparison utility does not perform or imply analysis.")
        (args.output / "comparison.json").write_text(json.dumps(report, indent=2))
        print("Comparison saved. Differences were recorded without changing inputs or statistics.")
    except Exception as error:
        (args.output / "FAILED.json").write_text(json.dumps({"status": "FAILED", "error": str(error)}, indent=2))
        raise


if __name__ == "__main__":
    main()

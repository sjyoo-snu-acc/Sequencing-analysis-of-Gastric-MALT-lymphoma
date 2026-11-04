#!/usr/bin/env python3
"""Static public-payload scan. Reports locations/rules, never matched secret values."""
import argparse
import ast
import hashlib
import json
import re
from pathlib import Path

RULES = {
    "personal_absolute_path": re.compile(r"(?:[A-Za-z]:[/\\]Users[/\\]|/Users/|/home/[^/\s]+/|/workspace/scratch/)", re.I),
    "study_specimen_identifier": re.compile(r"\b\d{2}[-.]\d{4}[-.][NT][-.]W(?:E|T)S(?:[-.]\d+)?\b"),
    "reported_case_alias": re.compile(r"\bcase_\d{2}\b"),
    "credential_assignment": re.compile(r"\b(?:idpw|password|passwd|userpwd|api_key|access_token)\s*(?:=|:|<-)\s*['\"][^'\"]+['\"]", re.I),
    "private_endpoint": re.compile(r"(?:sftp://|bio-ui\d+\.|/bio/scratch/)", re.I),
    "credential_token": re.compile(r"(?:gh[pousr]_[A-Za-z0-9]{25,}|AKIA[0-9A-Z]{16}|sk-[A-Za-z0-9]{32,})"),
}
BLOCKED = {".xlsx", ".xls", ".zip", ".bam", ".cram", ".vcf", ".pem", ".key", ".rds", ".rdata"}
BASE_R_PACKAGES = {"base", "utils", "stats", "grDevices", "graphics", "methods", "tools"}
DECLARED = BASE_R_PACKAGES | {"yaml", "jsonlite", "digest", "DESeq2", "maftools", "clusterProfiler", "biomaRt", "AnnotationDbi", "BiocManager", "codetools"}


def audit(root):
    findings = []; dependencies = set(); files = []
    ignored = {"data", "results", "test-output", "internal-validation", ".git", "__pycache__"}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root)
        if not path.is_file() or any(part in ignored for part in relative.parts):
            continue
        files.append({"path": str(relative), "bytes": path.stat().st_size, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        if path.suffix.lower() in BLOCKED or path.name.startswith((".env", ".Renviron")):
            findings.append({"path": str(relative), "rule": "restricted_file_type"})
        if path.suffix in {".R", ".py", ".md", ".yaml", ".json", ".tsv", ".txt"}:
            text = path.read_text(encoding="utf-8-sig")
            if path.suffix == ".py":
                ast.parse(text, filename=str(relative))
            if path.suffix == ".R":
                dependencies.update(re.findall(r"\b([A-Za-z][A-Za-z0-9._]+)::", text))
                # Namespace resolution is explicit; external files sourced by runners exist.
                if relative.parts[0] in {"WES", "WTS", "Microbiome"}:
                    for source in re.findall(r'file.path\(root, "R", "([^"]+)"\)', text):
                        if not (root / "R" / source).is_file():
                            findings.append({"path": str(relative), "rule": "missing_source_file"})
            if path.name == "audit_public_release.py":
                continue  # Do not detect the scanner's own regex literals as actual secrets.
            for row, line in enumerate(text.splitlines(), 1):
                for rule, pattern in RULES.items():
                    if pattern.search(line):
                        findings.append({"path": str(relative), "line": row, "rule": rule})
    for package in sorted(dependencies - DECLARED):
        findings.append({"package": package, "rule": "undeclared_namespace"})
    return {"status": "PASS" if not findings else "FAIL", "files_inspected": len(files),
            "namespace_packages": sorted(dependencies), "findings": findings,
            "limits": "Heuristic static scan and code review do not establish license rights or correctness of future private inputs.",
            "files": files}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default=".", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    result = audit(args.root.resolve())
    if args.report:
        if args.report.exists():
            raise ValueError("Refusing to overwrite static audit report")
        args.report.write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != "files"}, indent=2))
    raise SystemExit(0 if result["status"] == "PASS" else 1)

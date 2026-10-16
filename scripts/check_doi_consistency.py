#!/usr/bin/env python3
"""Keep all repository model/ML citations consistent with the 46 DOI registry.

This offline test catches wrong but syntactically valid DOIs in training cohorts,
historic benchmarks and regenerated artifact metadata. DOI identity against
Crossref is checked separately by scripts/verify_model_dois.py.
"""
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CATALOG = {m["id"]: m for m in json.loads((ROOT / "models/catalog.json").read_text())["models"]}
PATTERN = re.compile(r"doi:\s*(10\.\d{4,9}/[^\s,;]+)", re.I)
errors = []
checked = 0


def assert_reference(model_id, value, location):
    global checked
    found = PATTERN.search(value or "")
    if not found:
        return
    if model_id not in CATALOG:
        return  # A historical record may concern a retired model.
    checked += 1
    expected = CATALOG[model_id]["doi"].casefold()
    observed = found.group(1).rstrip(".").casefold()
    if observed != expected:
        errors.append(f"{location}: {model_id}: {observed} != {expected}")


metadata = json.loads((ROOT / "models/metadata.json").read_text())["models"]
for id_, row in CATALOG.items():
    if row["doi"].casefold() != metadata[id_]["doi"].casefold():
        errors.append(f"models/metadata.json: {id_} DOI differs from catalog")

training = json.loads((ROOT / "ml/training-populations.json").read_text())["models"]
for id_, row in training.items():
    assert_reference(id_, row.get("source", ""), "ml/training-populations.json")

registry = json.loads((ROOT / "ml/registry.json").read_text())["artifacts"]
for artifact in registry:
    assert_reference(
        artifact["baseModelId"],
        artifact.get("methodology", {}).get("sourcePopulation", ""),
        "ml/registry.json"
    )

for path in ("ml/validation/source-audit-retraining.csv",
             "ml/validation/c0-c0post-all-models-2026-09-24.csv"):
    with (ROOT / path).open(newline="", encoding="utf-8") as file:
        for row in csv.DictReader(file):
            assert_reference(row.get("base_model", ""), row.get("population_source", ""), path)

if errors:
    raise SystemExit("\n".join(errors))
print(f"DOI provenance coherence: 46 models; {checked} DOI-tagged ML references checked; no mismatches.")

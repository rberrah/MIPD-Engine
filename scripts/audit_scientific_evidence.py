#!/usr/bin/env python3
"""Audit scientific MIPD evidence without conflating levels of validation.

No network or third-party packages. This script generates a complete inventory
from the checked-in catalogue and independently authored evidence records.
Actual numerical claims are checked by tests/source_equation_reproduction.R.
"""

import argparse
import json
import os
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEVELS = {"PARTIAL_SOURCE_EQUATIONS", "FULL_SOURCE_REPRODUCTION"}
REQUIRED = {
    "modelId", "level", "primaryDoi", "primaryUrl", "sourceLocation",
    "population", "sourceValues", "mathematicalOracle", "testPath",
    "scenarios", "verifiedScope", "exclusions", "reviewedOn"
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        default="scientific-validation/coverage.generated.json",
        help="Output path relative to repository root",
    )
    args = parser.parse_args()
    catalogue = json.loads((ROOT / "models/catalog.json").read_text(encoding="utf-8"))["models"]
    manifest = json.loads(
        (ROOT / "scientific-validation/source-evidence.json").read_text(encoding="utf-8")
    )
    if manifest["schemaVersion"] != 1:
        raise SystemExit("Unsupported scientific evidence manifest schema")
    ids = [m["id"] for m in catalogue]
    if len(ids) != 46 or len(set(ids)) != len(ids):
        raise SystemExit(f"Expected 46 unique shipped models, found {len(ids)}")
    by_id = {m["id"]: m for m in catalogue}

    evidence = {}
    for item in manifest["entries"]:
        missing = REQUIRED.difference(item)
        if missing:
            raise SystemExit(f"{item.get('modelId')}: missing evidence: {sorted(missing)}")
        mid = item["modelId"]
        if mid in evidence or mid not in by_id:
            raise SystemExit(f"Duplicate or unknown evidence entry: {mid}")
        if item["level"] not in LEVELS:
            raise SystemExit(f"{mid}: invalid evidence level {item['level']}")
        if item["primaryDoi"].lower() != by_id[mid]["doi"].lower():
            raise SystemExit(f"{mid}: primary source DOI does not match published catalogue")
        if not item["primaryUrl"].startswith("https://"):
            raise SystemExit(f"{mid}: primary source must use an HTTPS URL")
        if not isinstance(item["sourceValues"], dict) or not item["sourceValues"]:
            raise SystemExit(f"{mid}: independently transcribed numeric values required")
        for field in ("scenarios", "verifiedScope", "exclusions"):
            if not isinstance(item[field], list) or not item[field]:
                raise SystemExit(f"{mid}: {field} must be a nonempty list")
        if not (ROOT / item["testPath"]).is_file():
            raise SystemExit(f"{mid}: claimed source reproduction lacks runnable test")
        try:
            date.fromisoformat(item["reviewedOn"])
        except ValueError as error:
            raise SystemExit(f"{mid}: invalid review date") from error
        if item["level"] == "FULL_SOURCE_REPRODUCTION":
            # A future promotion requires independently reviewed entire-model
            # scope; zero currently qualify. Separate adjudication is required.
            if item.get("independentReviewer") in (None, "") or item["exclusions"]:
                raise SystemExit(f"{mid}: cannot claim FULL with unreviewed scope/exclusions")
        evidence[mid] = item

    report = []
    for model in catalogue:
        mid = model["id"]
        ev = evidence.get(mid)
        verification = model.get("verificationStatus", "NOT_TESTED")
        if verification == "NUMERICALLY_REPRODUCED" and (
            ev is None or ev["level"] != "FULL_SOURCE_REPRODUCTION"
        ):
            raise SystemExit(
                f"{mid}: catalog promotes NUMERICALLY_REPRODUCED without full evidence"
            )
        report.append({
            "id": mid,
            "primaryDoi": model["doi"],
            "implementationStatus": model["implementationStatus"],
            "analysisEligible": bool(model["analysisEligible"]),
            "softwareVerification": verification,
            "sourceReproduction": ev["level"] if ev else "NOT_REPRODUCED",
            "evidencePath": "scientific-validation/source-evidence.json" if ev else None,
            "testPath": ev["testPath"] if ev else None,
            "unverifiedScope": ev["exclusions"] if ev else ["Independent source reproduction not yet performed"],
            "externalValidation": model.get("validationStatus", "NONE"),
        })

    summary = {
        "catalogueModels": len(report),
        "softwareTested": sum(r["softwareVerification"] != "NOT_TESTED" for r in report),
        "independentSourceEquationsPartial": sum(
            r["sourceReproduction"] == "PARTIAL_SOURCE_EQUATIONS" for r in report
        ),
        "fullSourceReproductions": sum(
            r["sourceReproduction"] == "FULL_SOURCE_REPRODUCTION" for r in report
        ),
        "pendingSourceReproduction": sum(
            r["sourceReproduction"] == "NOT_REPRODUCED" for r in report
        ),
        "mapEligibleTotal": sum(r["analysisEligible"] for r in report),
        "mapEligibleWithPartialOrFullSourceEvidence": sum(
            r["analysisEligible"] and r["sourceReproduction"] != "NOT_REPRODUCED"
            for r in report
        ),
        "externallyValidated": sum(r["externalValidation"] == "EXTERNAL" for r in report),
    }
    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps({"schemaVersion": 1, "summary": summary, "models": report}, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2))
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as f:
            f.write(
                "\n### Scientific qualification (not clinical validation)\n\n"
                "| Scope | Count |\n|---|---:|\n"
                f"| Shipped library models | {summary['catalogueModels']} |\n"
                f"| Independently checked, partial source-equation scope | {summary['independentSourceEquationsPartial']} |\n"
                f"| Full source reproductions | {summary['fullSourceReproductions']} |\n"
                f"| Pending independent reproduction | {summary['pendingSourceReproduction']} |\n"
                f"| MAP-eligible models with source-equation evidence | {summary['mapEligibleWithPartialOrFullSourceEvidence']}/{summary['mapEligibleTotal']} |\n"
                f"| Externally validated | {summary['externallyValidated']} |\n\n"
                "Source-equation comparisons are not external patient-level validation.\n"
            )


if __name__ == "__main__":
    main()

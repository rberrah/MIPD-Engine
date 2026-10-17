#!/usr/bin/env python3
"""Audit ML evidence without treating source compatibility as validation.

Reads public aggregate synthetic metrics only; never loads patient data or
executes/rewrites RDS. Fail-closed policy is inherited from the MIPD catalogue.
"""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
source = json.loads((ROOT / "ml/catalog-summary.generated.json").read_text())
catalogue = {
    row["id"]: row
    for row in json.loads((ROOT / "models/catalog.json").read_text())["models"]
}
rows = []
errors = []
for artifact in source["artifacts"]:
    mid = artifact["modelId"]
    v = artifact.get("validation") or {}
    thresholds = v.get("thresholds") or {}
    cv = v.get("repeatedCv") or {}
    holdout = v.get("untouchedHoldout") or {}
    compatible = bool(artifact["hashMatches"])
    available = bool(artifact["availableInTdm"])
    eligible = bool(catalogue[mid]["analysisEligible"])
    if available != (eligible and compatible):
        errors.append(f"{mid}: incorrect software source-compatibility gate")
    if not compatible and available:
        errors.append(f"{mid}: obsolete source SHA artifact exposed to TDM")

    def within_spec(metrics: dict) -> bool:
        if not all(isinstance(metrics.get(key), (int, float))
                   for key in ("relativeRmsePct", "relativeBiasPct", "within20Pct")):
            return False
        return (
            metrics["relativeRmsePct"] <= thresholds.get("relativeRmsePct", 15)
            and abs(metrics["relativeBiasPct"]) <= thresholds.get("absoluteBiasPct", 10)
            and metrics["within20Pct"] >= thresholds.get("within20Pct", 80)
        )
    for label, result in (("repeatedCv", cv), ("untouchedHoldout", holdout)):
        if result.get("passed") is True and not within_spec(result):
            errors.append(f"{mid}: {label} is marked passed despite violating a numeric threshold")
    internal_pass = cv.get("passed") is True and holdout.get("passed") is True
    clinical = v.get("realPatient", {}).get("status") == "validated"
    if clinical:
        errors.append(f"{mid}: new external patient validation claim requires independent review")
    rows.append({
        "modelId": mid,
        "articleImplementationStatus": catalogue[mid]["implementationStatus"],
        "sourceHashCompatible": compatible,
        "mapEligible": eligible,
        "technicallyAvailableInTdm": available,
        "repeatedCvPassed": cv.get("passed") is True,
        "untouchedHoldoutPassed": holdout.get("passed") is True,
        "internalResearchPass": internal_pass,
        "holdoutRelativeRmsePct": holdout.get("relativeRmsePct"),
        "holdoutP20Pct": holdout.get("within20Pct"),
        "alternatePopPkPassed": (v.get("alternatePopPk") or {}).get("passed") is True,
        "externalPatientValidation": clinical,
    })
stale = [r["modelId"] for r in rows if not r["sourceHashCompatible"]]
available = [r for r in rows if r["technicallyAvailableInTdm"]]
internal_pass = [r["modelId"] for r in available if r["internalResearchPass"]]
underperforming = [r["modelId"] for r in available if not r["internalResearchPass"]]
if len(rows) != 46 or len(stale) != 12 or len(available) != 12 or len(internal_pass) != 3:
    errors.append("ML status baseline changed; audit updated models/metrics before changing the qualified baseline")
goti = next(r for r in rows if r["modelId"] == "vanco_goti")
if not goti["sourceHashCompatible"] or goti["internalResearchPass"] or goti["externalPatientValidation"]:
    errors.append("Goti retraining is source compatible, but does not pass internal or external validation")
result = {
    "schema": "mipd-ml-evidence.v1",
    "policy": "Technical source compatibility, synthetic CV/holdout performance, alternate PopPK transportability and real-patient validation are distinct. No patient data used.",
    "artifactCount": len(rows),
    "staleSourceModels": stale,
    "technicallyAvailableCount": len(available),
    "internalResearchPassModels": internal_pass,
    "belowInternalThresholdModels": underperforming,
    "externallyValidatedCount": sum(r["externalPatientValidation"] for r in rows),
    "errors": errors,
    "models": rows
}
out = ROOT / "scientific-validation/ml-evidence.generated.json"
out.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n")
if errors:
    raise SystemExit("\n".join(errors))
print(
    f"ML evidence: {len(rows)} artifacts; {len(stale)} stale/blocked; "
    f"{len(available)} technically available; {len(internal_pass)} pass internal "
    f"CV and holdout; {len(underperforming)} under thresholds; zero external clinical validations."
)

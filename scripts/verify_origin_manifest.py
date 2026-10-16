#!/usr/bin/env python3
import json, pathlib, subprocess, sys

root = pathlib.Path(__file__).resolve().parents[1]
manifest = json.loads((root / "migration" / "source-manifest.json").read_text(encoding="utf-8"))

failures = []
if manifest.get("copiedCount") != len(manifest["copied"]):
    failures.append("copiedCount does not match copied entries")
if manifest.get("adaptedCount") != len(manifest.get("adapted", [])):
    failures.append("adaptedCount does not match adapted entries")

copied_targets = {item["target"] for item in manifest["copied"]}
adapted_targets = {item["target"] for item in manifest.get("adapted", [])}
if copied_targets & adapted_targets:
    failures.append("targets cannot be both copied and adapted")

for item in manifest["copied"]:
    path = root / item["target"]
    if not path.is_file():
        failures.append(f"missing: {item['target']}")
        continue
    digest = subprocess.check_output(
        ["git", "hash-object", f"--path={item['target']}", item["target"]],
        cwd=root,
        text=True,
    ).strip()
    if digest != item["gitBlob"]:
        failures.append(f"hash mismatch: {item['target']} expected {item['gitBlob']} got {digest}")

for item in manifest.get("adapted", []):
    path = root / item["target"]
    if not path.is_file():
        failures.append(f"missing adapted target: {item['target']}")
        continue
    baseline = subprocess.run(
        ["git", "cat-file", "-e", f"{item['gitBlob']}^{{blob}}"],
        cwd=root,
        capture_output=True,
    )
    if baseline.returncode:
        failures.append(f"missing origin blob for adapted target: {item['target']} ({item['gitBlob']})")

if failures:
    print("\n".join(failures), file=sys.stderr)
    raise SystemExit(1)

print(f"Origin parity OK: {len(manifest['copied'])} copied files match {manifest['sourceRepository']}@{manifest['sourceCommit']}.")
print(f"Origin baselines retained for {len(manifest.get('adapted', []))} intentionally adapted files.")
print(f"Pending migration entries: {len(manifest['pending'])}.")

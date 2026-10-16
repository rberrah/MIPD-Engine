#!/usr/bin/env python3
import json, pathlib, subprocess, sys

root = pathlib.Path(__file__).resolve().parents[1]
manifest = json.loads((root / "migration" / "source-manifest.json").read_text(encoding="utf-8"))

failures = []
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

if failures:
    print("\n".join(failures), file=sys.stderr)
    raise SystemExit(1)

print(f"Origin parity OK: {len(manifest['copied'])} copied files match {manifest['sourceRepository']}@{manifest['sourceCommit']}.")
print(f"Pending migration entries: {len(manifest['pending'])}.")

#!/usr/bin/env python3
"""Offline contract test for the post-deploy read-only public HTTP verifier."""
import contextlib
import importlib.util
import io
import json
import pathlib
import sys

root = pathlib.Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("live_check", root / "scripts" / "check_live_deployment.py")
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)
sha = "a" * 40
valid = {
    "schema": "mipd-public-release.v1",
    "target": "tdmhub/MIPD_Engine",
    "appId": "17776160",
    "sourceSha": sha,
    "modelCount": 46,
    "allowCustomModels": False,
    "enablePmetrics": False,
}

def exercise(profile, html=b"<html><body>Shiny app</body></html>"):
    def fake_fetch(url, timeout):
        if "deployment-release.json" in url:
            return json.dumps(profile).encode("utf-8"), "application/json"
        return html, "text/html"
    check.fetch = fake_fetch
    sys.argv = ["check_live_deployment.py", "--expected-sha", sha, "--attempts", "1"]
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        return check.main()

assert exercise(valid) == 0
for key, value in [
    ("sourceSha", "b" * 40),
    ("enablePmetrics", True),
    ("allowCustomModels", True),
    ("modelCount", 45),
    ("target", "someone-else/mirror"),
    ("schema", "unknown"),
]:
    bad = {**valid, key: value}
    assert exercise(bad) == 1, f"Unsafe or stale public release accepted: {key}"
assert exercise(valid, html=b"service unavailable") == 1
print("Live-release checker contract PASS: exact source SHA, public safety, stale and invalid responses.")

#!/usr/bin/env python3
"""Read-only public Shiny deployment verification. No patient data are sent."""
import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

DEFAULT = "https://tdmhub.shinyapps.io/MIPD_Engine/"
MAX_BYTES = 4_000_000


def fetch(url, timeout):
    request = urllib.request.Request(
        url, headers={"User-Agent": "MIPD-Engine-release-verification/1.0",
                      "Accept": "text/html,application/json"}
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        payload = response.read(MAX_BYTES + 1)
        if len(payload) > MAX_BYTES:
            raise ValueError("Oversized public response")
        if response.status != 200:
            raise ValueError(f"Unexpected HTTP {response.status} from {url}")
        return payload, response.headers.get("Content-Type", "")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--base", default=DEFAULT)
    p.add_argument("--expected-sha", required=True)
    p.add_argument("--attempts", type=int, default=6)
    p.add_argument("--timeout", type=int, default=30)
    args = p.parse_args()
    if not re.fullmatch(r"[a-f0-9]{40}", args.expected_sha):
        p.error("--expected-sha must be a 40-character immutable commit SHA")
    base = args.base.rstrip("/") + "/"
    if not base.startswith("https://"):
        p.error("The public app must use HTTPS")
    errors = []
    # A cold Shinyapps instance can require activation before serving www files.
    for attempt in range(1, args.attempts + 1):
        try:
            html, _ = fetch(base + "?lang=en", args.timeout)
            decoded = html.decode("utf-8", errors="replace")
            if not re.search(r"<html\\b", decoded, re.IGNORECASE):
                raise ValueError("Public app did not respond with an HTML document")
            if not re.search(r"(shiny|shinyapps)", decoded, re.IGNORECASE):
                raise ValueError("Public response is not recognizable as a Shiny app")
            release_data, _ = fetch(
                urllib.parse.urljoin(base, "deployment-release.json"),
                args.timeout
            )
            release = json.loads(release_data)
            expected = {
                "schema": "mipd-public-release.v1",
                "target": "tdmhub/MIPD_Engine",
                "appId": "17776160",
                "sourceSha": args.expected_sha,
                "modelCount": 46,
                "allowCustomModels": False,
                "enablePmetrics": False,
            }
            for key, value in expected.items():
                if release.get(key) != value:
                    raise ValueError(f"Public release does not match {key}: "
                                     f"expected {value!r}, received {release.get(key)!r}")
            print("PUBLIC MIPD RELEASE PASS", release["sourceSha"])
            print("Verified active app HTML and immutable bundle fingerprint "
                  "(46 models, custom C++ and Pmetrics disabled).")
            return 0
        except (ValueError, OSError, json.JSONDecodeError, UnicodeError) as exc:
            errors.append(f"Attempt {attempt}: {exc}")
            if attempt < args.attempts:
                time.sleep(6)
    print("PUBLIC MIPD RELEASE FAIL:\n" + "\n".join(errors), file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())

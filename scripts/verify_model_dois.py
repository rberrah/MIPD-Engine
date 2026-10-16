#!/usr/bin/env python3
"""Resolve and verify every shipped model DOI against authoritative Crossref records.

Only bibliographic alignment is assessed here. A valid DOI is not proof that the
implemented PK equations faithfully reproduce the cited article. No PDFs or
patient data are downloaded, committed or inspected.
"""
import argparse
import json
import re
import sys
import time
import unicodedata
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
API = "https://api.crossref.org/works/"
HEADERS = {
    "Accept": "application/json",
    "User-Agent": "MIPD-Engine-DOI-audit/1.0 (https://github.com/rberrah/MIPD-Engine)",
}


def normalized(value):
    plain = unicodedata.normalize("NFKD", value or "")
    return re.sub(r"[^a-z0-9]", "", "".join(c for c in plain if not unicodedata.combining(c)).lower())


def crossref(doi, retries=4):
    url = API + quote(doi, safe="")
    for attempt in range(retries):
        try:
            request = Request(url, headers=HEADERS)
            with urlopen(request, timeout=25) as response:
                return json.load(response)["message"], None
        except HTTPError as err:
            if err.code == 404:
                return None, "Crossref DOI not found"
            if err.code not in (429, 500, 502, 503, 504) or attempt == retries - 1:
                return None, f"Crossref HTTP {err.code}"
        except (URLError, TimeoutError, OSError, json.JSONDecodeError) as err:
            if attempt == retries - 1:
                return None, f"Crossref unavailable ({type(err).__name__})"
        time.sleep(1.5 * (attempt + 1))
    return None, "Crossref unavailable"


def source_year(entry):
    match = re.search(r"\b(19\d{2}|20\d{2})\s*;", entry["citation"])
    if match is None:
        raise ValueError(f"Missing year: {entry['id']}")
    return int(match.group(1))


def published_year(message):
    for field in ("published-print", "published", "published-online", "issued", "created"):
        parts = message.get(field, {}).get("date-parts", [[]])
        if parts and parts[0] and isinstance(parts[0][0], int):
            return parts[0][0]
    return None


def check_entry(entry, record, error):
    expected_author = entry["model"]
    status = "UNVERIFIED"
    reasons = []
    got_author = None
    got_year = None
    title = None
    if error:
        reasons.append(error)
    elif record:
        authors = record.get("author", [])
        got_author = authors[0].get("family", "") if authors else ""
        got_year = published_year(record)
        titles = record.get("title", [])
        title = titles[0] if titles else ""
        actual_doi = record.get("DOI", "")
        if normalized(actual_doi) != normalized(entry["doi"]):
            reasons.append(f"DOI mismatch: {actual_doi}")
        if normalized(got_author) != normalized(expected_author):
            reasons.append(f"First author mismatch: expected {expected_author}, got {got_author}")
        year = source_year(entry)
        if got_year is None or abs(year - got_year) > 1:
            reasons.append(f"Year mismatch: expected {year}, got {got_year}")
        if not title:
            reasons.append("Missing registered article title")
        if reasons:
            status = "MISMATCH"
        else:
            status = "VERIFIED"
    return {
        "id": entry["id"],
        "doi": entry["doi"],
        "modelAuthor": expected_author,
        "citation": entry["citation"],
        "registeredFirstAuthor": got_author,
        "registeredPublicationYear": got_year,
        "registeredTitle": title,
        "registeredJournal": (record or {}).get("container-title", [None])[0]
            if (record or {}).get("container-title") else None,
        "status": status,
        "issues": reasons,
        "resolver": "Crossref REST API",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="doi-audit-report.json")
    parser.add_argument("--strict", action="store_true", help="Fail if any unmatched or unresolved model DOI")
    args = parser.parse_args()
    catalog = json.loads((ROOT / "models/catalog.json").read_text(encoding="utf-8"))["models"]
    if len(catalog) != 46:
        raise SystemExit(f"Expected exactly 46 active models, got {len(catalog)}")
    cache = {}
    result = []
    for i, model in enumerate(catalog, 1):
        doi = model["doi"].lower().strip()
        if doi not in cache:
            cache[doi] = crossref(doi)
        record, error = cache[doi]
        checked = check_entry(model, record, error)
        result.append(checked)
        print(f"[{i:02d}/46] {checked['status']:10s} {model['id']}  {doi}", flush=True)
        if checked["issues"]:
            print("   " + "; ".join(checked["issues"]), flush=True)
    counts = {s: sum(r["status"] == s for r in result)
              for s in ("VERIFIED", "MISMATCH", "UNVERIFIED")}
    summary = {"catalogueModels": len(result), "uniqueDois": len(cache), **counts}
    output = Path(args.output)
    output.write_text(json.dumps({"schemaVersion": 1, "summary": summary, "models": result},
                                 indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2), flush=True)
    if args.strict and (counts["MISMATCH"] or counts["UNVERIFIED"]):
        raise SystemExit("DOI audit found bibliographic mismatches or unresolved records")


if __name__ == "__main__":
    main()

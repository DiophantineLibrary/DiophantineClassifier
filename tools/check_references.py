r"""Reference checking/annotation pipeline for the Diophantine Library.

Validates ``diophantine_classifier/data/references.bib`` against the family
registry and against locally downloaded papers, and writes a status report.

Checks
------

structural (always):
  - the BibTeX file parses and entries carry author/title/year;
  - every reference key used by a family exists in the bibliography and
    carries a nonempty ``why`` annotation;
  - ``doi``/``url``/``eprint`` fields are well-formed;
  - unused bibliography entries are reported (warning only).

local documents (``references/pdf/<key>.pdf``):
  - when a PDF has been downloaded for an entry, its text (extracted with
    ``pdftotext`` if available) is checked against the entry: most long
    title words and at least one author surname must appear.  This catches
    wrong downloads and mismatched metadata.

online (``--online``):
  - DOIs must resolve at ``https://doi.org/`` and ``url`` fields must be
    reachable.

Usage::

    sage -python tools/check_references.py [--online] [--quiet]

Writes ``references/REPORT.md`` and exits nonzero on hard errors (unknown
keys, malformed entries, empty annotations).
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from diophantine_classifier.references import (bibliography, _delatex)  # noqa: E402
from diophantine_classifier.registry import families  # noqa: E402

PDF_DIR = os.path.join(ROOT, "references", "pdf")
REPORT = os.path.join(ROOT, "references", "REPORT.md")

REQUIRED_FIELDS = ("author", "title", "year")
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
ARXIV_RE = re.compile(r"^(\d{4}\.\d{4,5}|[a-z-]+/\d{7})(v\d+)?$")


def usage_map():
    """Map bib key -> list of (family slug, why annotation)."""
    used = {}
    for fam in families().values():
        for ref in fam.references:
            used.setdefault(ref["key"], []).append((fam.slug,
                                                    ref.get("why", "")))
    return used


def surnames(author_field):
    """Rough surname list from a BibTeX author field."""
    names = []
    for person in _delatex(author_field).split(" and "):
        person = person.strip().rstrip(",")
        if not person:
            continue
        # drop suffixes like Jr.
        parts = [p for p in re.split(r"[\s,]+", person)
                 if p and not p.endswith(".")]
        if parts:
            names.append(parts[-1].lower())
    return names


def check_structure(bib, used):
    """Structural validation; returns (errors, warnings)."""
    errors, warnings = [], []
    for key, annotations in used.items():
        if key not in bib:
            errors.append(f"unknown reference key {key!r} used by "
                          f"{[slug for slug, _ in annotations]}")
        for slug, why in annotations:
            if not why.strip():
                errors.append(f"{slug}: reference {key} lacks a 'why' "
                              "annotation")
    for key, entry in bib.items():
        fields = entry["fields"]
        for req in REQUIRED_FIELDS:
            if req not in fields:
                errors.append(f"{key}: missing required field {req!r}")
        doi = fields.get("doi")
        if doi and not DOI_RE.match(doi):
            errors.append(f"{key}: malformed doi {doi!r}")
        url = fields.get("url")
        if url and not url.startswith(("http://", "https://")):
            errors.append(f"{key}: malformed url {url!r}")
        eprint = fields.get("eprint")
        if eprint and not ARXIV_RE.match(eprint):
            errors.append(f"{key}: malformed arXiv id {eprint!r}")
        if key not in used:
            warnings.append(f"{key}: not referenced by any family")
    return errors, warnings


def pdf_status(key, entry):
    """Check a locally downloaded PDF against its bibliography entry."""
    path = os.path.join(PDF_DIR, f"{key}.pdf")
    if not os.path.exists(path):
        return "missing", "no local PDF"
    if os.path.getsize(path) < 10_000:
        return "suspect", "file smaller than 10KB"
    if shutil.which("pdftotext") is None:
        return "present", "downloaded (install pdftotext for text checks)"
    try:
        out = subprocess.run(
            ["pdftotext", "-l", "3", "-q", path, "-"],
            capture_output=True, timeout=60)
        text = out.stdout.decode("utf-8", errors="ignore").lower()
    except Exception as err:
        return "suspect", f"pdftotext failed: {err}"
    if not text.strip():
        return "suspect", "no extractable text (scanned image?)"
    fields = entry["fields"]
    title_words = [w for w in re.findall(r"[a-z]{5,}",
                                         _delatex(fields.get("title", "")).lower())]
    hits = sum(1 for w in set(title_words) if w in text)
    coverage = hits / max(1, len(set(title_words)))
    author_hit = any(s in text for s in surnames(fields.get("author", "")))
    if coverage >= 0.6 and author_hit:
        return "verified", f"title-word coverage {coverage:.0%}, author found"
    return ("mismatch",
            f"title-word coverage {coverage:.0%}, author "
            f"{'found' if author_hit else 'NOT found'} — check the download")


def online_status(entry):
    """Check that doi/url resolve (network access required)."""
    results = []
    fields = entry["fields"]
    for label, target in (("doi", "https://doi.org/" + fields["doi"]
                           if "doi" in fields else None),
                          ("url", fields.get("url"))):
        if target is None:
            continue
        req = urllib.request.Request(target, method="HEAD",
                                     headers={"User-Agent": "diophlib-check"})
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                results.append((label, "ok" if resp.status < 400
                                else f"HTTP {resp.status}"))
        except Exception as err:
            results.append((label, f"unreachable ({err})"))
    return results


def write_report(bib, used, pdf, online, errors, warnings):
    lines = ["# Reference status report", "",
             "Generated by `tools/check_references.py`. Download papers to "
             "`references/pdf/<key>.pdf` (legally!) and re-run to verify "
             "them against the bibliography.", ""]
    n_doi = sum(1 for e in bib.values() if "doi" in e["fields"])
    n_url = sum(1 for e in bib.values() if "url" in e["fields"])
    n_ver = sum(1 for s, _ in pdf.values() if s == "verified")
    lines += [f"- entries: {len(bib)}; used by families: {len(used)}",
              f"- with DOI: {n_doi}; with free URL: {n_url}; "
              f"local PDFs verified: {n_ver}", ""]
    if errors:
        lines += ["## Errors", ""] + [f"- {e}" for e in errors] + [""]
    if warnings:
        lines += ["## Warnings", ""] + [f"- {w}" for w in warnings] + [""]
    lines += ["## Entries", "",
              "| key | families | doi | free url | local pdf |",
              "|-----|----------|-----|----------|-----------|"]
    for key in sorted(bib):
        fields = bib[key]["fields"]
        fams = ", ".join(slug for slug, _ in used.get(key, [])) or "—"
        doi = "yes" if "doi" in fields else "TODO"
        url = "yes" if "url" in fields else "TODO"
        status, note = pdf.get(key, ("missing", ""))
        cell = status if status == "missing" else f"{status} ({note})"
        if online and key in online:
            cell += "; online: " + "; ".join(f"{l}: {s}"
                                             for l, s in online[key])
        lines.append(f"| {key} | {fams} | {doi} | {url} | {cell} |")
    todo_doi = sorted(k for k, e in bib.items() if "doi" not in e["fields"])
    todo_url = sorted(k for k, e in bib.items() if "url" not in e["fields"])
    todo_pdf = sorted(k for k in bib if pdf.get(k, ("missing",))[0] == "missing")
    lines += ["", "## TODO", "",
              f"- add DOIs ({len(todo_doi)}): " + ", ".join(todo_doi),
              f"- find legally free URLs ({len(todo_url)}): "
              + ", ".join(todo_url),
              f"- download PDFs ({len(todo_pdf)}): " + ", ".join(todo_pdf),
              ""]
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    with open(REPORT, "w") as fobj:
        fobj.write("\n".join(lines))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--online", action="store_true",
                        help="also check that DOIs and URLs resolve")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    bib = bibliography()
    used = usage_map()
    errors, warnings = check_structure(bib, used)
    pdf = {key: pdf_status(key, entry) for key, entry in bib.items()}
    online = {}
    if args.online:
        for key, entry in bib.items():
            res = online_status(entry)
            if res:
                online[key] = res
                bad = [f"{key} {l}: {s}" for l, s in res if s != "ok"]
                errors.extend(bad)
    write_report(bib, used, pdf, online, errors, warnings)

    if not args.quiet:
        n_ver = sum(1 for s, _ in pdf.values() if s == "verified")
        n_mis = sum(1 for s, _ in pdf.values() if s == "mismatch")
        print(f"{len(bib)} entries; {len(errors)} errors, "
              f"{len(warnings)} warnings; PDFs: {n_ver} verified, "
              f"{n_mis} mismatched, "
              f"{sum(1 for s, _ in pdf.values() if s == 'missing')} missing")
        for e in errors:
            print(f"ERROR: {e}")
        print(f"report: {os.path.relpath(REPORT, os.getcwd())}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

r"""
Reference checking/annotation pipeline for the Diophantine Library.

Validates ``diophantine_classifier/data/references.bib`` against the family
registry and against locally downloaded papers, and writes a status report.

Because different people run this on different machines, verification status
is **monotone**: successes are recorded in a committed ledger
(``references/status.yaml``) and are never downgraded by a run on a machine
that lacks the file — a missing local PDF reports as "missing here" while the
ledger still shows where and when it was verified.  A *local* mismatch never
erases a ledger verification either (your download may simply be the wrong
file); it is surfaced as a warning instead.

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
    title words and at least one author surname must appear.  A success is
    recorded in the ledger (write-once).

online (``--online``):
  - DOIs must resolve at ``https://doi.org/`` and ``url`` fields must be
    reachable; successes are recorded in the ledger (write-once).

Usage::

    sage -python tools/check_references.py [--online] [--quiet]

Writes ``references/REPORT.md``, updates ``references/status.yaml`` when new
verifications happen, and exits nonzero on hard errors (unknown keys,
malformed entries, empty annotations).
"""

import argparse
import datetime
import os
import re
import shutil
import subprocess
import sys
import urllib.request

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from diophantine_classifier.references import (bibliography, _delatex)  # noqa: E402
from diophantine_classifier.registry import families  # noqa: E402

PDF_DIR = os.path.join(ROOT, "references", "pdf")
REPORT = os.path.join(ROOT, "references", "REPORT.md")
LEDGER = os.path.join(ROOT, "references", "status.yaml")

REQUIRED_FIELDS = ("author", "title", "year")
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
ARXIV_RE = re.compile(r"^(\d{4}\.\d{4,5}|[a-z-]+/\d{7})(v\d+)?$")

#: ledger fields (write-once dates); statuses only ever upgrade
LEDGER_FIELDS = ("pdf_verified", "doi_ok", "url_ok")


def usage_map():
    r"""
    Map bib key -> list of ``(family slug, why annotation)``.
    """
    used = {}
    for fam in families().values():
        for ref in fam.references:
            used.setdefault(ref["key"], []).append((fam.slug,
                                                    ref.get("why", "")))
    return used


def load_ledger():
    r"""
    Read the committed verification ledger (empty dict if absent).
    """
    if os.path.exists(LEDGER):
        with open(LEDGER) as fobj:
            return yaml.safe_load(fobj) or {}
    return {}


def record(ledger, key, field):
    r"""
    Record a success in the ledger, write-once (monotone upgrades).

    Returns ``True`` if the ledger changed.  An existing record — including
    its original date — is never overwritten or removed.
    """
    if field not in LEDGER_FIELDS:
        raise ValueError(f"unknown ledger field {field!r}")
    entry = ledger.setdefault(key, {})
    if field in entry:
        return False
    entry[field] = datetime.date.today().isoformat()
    return True


def save_ledger(ledger):
    r"""
    Write the ledger with stable ordering (clean diffs).
    """
    os.makedirs(os.path.dirname(LEDGER), exist_ok=True)
    with open(LEDGER, "w") as fobj:
        fobj.write("# Verification ledger: successes recorded by "
                   "tools/check_references.py.\n"
                   "# Monotone by design: entries are write-once and never "
                   "downgraded by a run\n"
                   "# on a machine without the file. Commit this file.\n")
        yaml.safe_dump({k: dict(sorted(v.items()))
                        for k, v in sorted(ledger.items())},
                       fobj, sort_keys=True)


def surnames(author_field):
    r"""
    Rough surname list from a BibTeX author field.
    """
    names = []
    for person in _delatex(author_field).split(" and "):
        person = person.strip().rstrip(",")
        if not person:
            continue
        parts = [p for p in re.split(r"[\s,]+", person)
                 if p and not p.endswith(".")]
        if parts:
            names.append(parts[-1].lower())
    return names


def check_structure(bib, used):
    r"""
    Structural validation; returns ``(errors, warnings)``.
    """
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
    r"""
    Check a locally downloaded PDF against its bibliography entry.

    Returns ``(status, note)`` with status one of ``missing``, ``present``
    (no text tool), ``suspect``, ``mismatch``, ``verified``.  Only
    ``verified`` enters the ledger.
    """
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
    title_words = re.findall(r"[a-z]{5,}",
                             _delatex(fields.get("title", "")).lower())
    hits = sum(1 for w in set(title_words) if w in text)
    coverage = hits / max(1, len(set(title_words)))
    author_hit = any(s in text for s in surnames(fields.get("author", "")))
    if coverage >= 0.6 and author_hit:
        return "verified", f"title-word coverage {coverage:.0%}, author found"
    return ("mismatch",
            f"title-word coverage {coverage:.0%}, author "
            f"{'found' if author_hit else 'NOT found'} — check the download")


def online_status(entry):
    r"""
    Check that doi/url resolve (network access required).
    """
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


def _pdf_cell(key, status, note, ledger):
    r"""
    Render the local-PDF column, merging in the ledger's knowledge.
    """
    recorded = ledger.get(key, {}).get("pdf_verified")
    if status == "verified":
        return f"verified ({note})"
    cell = status if status == "missing" else f"{status} ({note})"
    if recorded:
        cell += f"; ledger: verified {recorded}"
    return cell


def write_report(bib, used, pdf, online, ledger, errors, warnings):
    r"""
    Write ``references/REPORT.md`` from this run merged with the ledger.
    """
    lines = ["# Reference status report", "",
             "Generated by `tools/check_references.py`. Download papers to "
             "`references/pdf/<key>.pdf` (legally!) and re-run to verify "
             "them against the bibliography. Verification status is "
             "monotone: successes are recorded in `status.yaml` and never "
             "downgraded by a machine that lacks the file.", ""]
    n_doi = sum(1 for e in bib.values() if "doi" in e["fields"])
    n_url = sum(1 for e in bib.values() if "url" in e["fields"])
    n_ver = sum(1 for k in bib
                if pdf.get(k, ("missing",))[0] == "verified"
                or "pdf_verified" in ledger.get(k, {}))
    lines += [f"- entries: {len(bib)}; used by families: {len(used)}",
              f"- with DOI: {n_doi}; with free URL: {n_url}; "
              f"PDFs verified (here or in ledger): {n_ver}", ""]
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
        rec = ledger.get(key, {})
        doi = "yes" if "doi" in fields else "TODO"
        if "doi_ok" in rec:
            doi += f" (ok {rec['doi_ok']})"
        url = "yes" if "url" in fields else "TODO"
        if "url_ok" in rec:
            url += f" (ok {rec['url_ok']})"
        status, note = pdf.get(key, ("missing", ""))
        cell = _pdf_cell(key, status, note, ledger)
        if online and key in online:
            cell += "; online: " + "; ".join(f"{l}: {s}"
                                             for l, s in online[key])
        lines.append(f"| {key} | {fams} | {doi} | {url} | {cell} |")
    todo_doi = sorted(k for k, e in bib.items() if "doi" not in e["fields"])
    todo_url = sorted(k for k, e in bib.items() if "url" not in e["fields"])
    todo_pdf = sorted(k for k in bib
                      if pdf.get(k, ("missing",))[0] == "missing"
                      and "pdf_verified" not in ledger.get(k, {}))
    lines += ["", "## TODO", "",
              f"- add DOIs ({len(todo_doi)}): " + ", ".join(todo_doi),
              f"- find legally free URLs ({len(todo_url)}): "
              + ", ".join(todo_url),
              f"- download and verify PDFs ({len(todo_pdf)}): "
              + ", ".join(todo_pdf),
              ""]
    os.makedirs(os.path.dirname(REPORT), exist_ok=True)
    with open(REPORT, "w") as fobj:
        fobj.write("\n".join(lines))


def main(argv=None):
    r"""
    Run the pipeline; returns the process exit code.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--online", action="store_true",
                        help="also check that DOIs and URLs resolve")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)

    bib = bibliography()
    used = usage_map()
    ledger = load_ledger()
    errors, warnings = check_structure(bib, used)
    pdf = {key: pdf_status(key, entry) for key, entry in bib.items()}

    changed = False
    for key, (status, _) in pdf.items():
        if status == "verified":
            changed |= record(ledger, key, "pdf_verified")
        elif status == "mismatch" and "pdf_verified" in ledger.get(key, {}):
            warnings.append(
                f"{key}: local file mismatches, but the ledger records a "
                f"verification ({ledger[key]['pdf_verified']}) — your local "
                "copy may be the wrong paper")

    online = {}
    if args.online:
        for key, entry in bib.items():
            res = online_status(entry)
            if res:
                online[key] = res
                for label, s in res:
                    if s == "ok":
                        changed |= record(ledger, key, f"{label}_ok")
                    else:
                        warnings.append(f"{key} {label}: {s}")

    if changed:
        save_ledger(ledger)
    write_report(bib, used, pdf, online, ledger, errors, warnings)

    if not args.quiet:
        n_ver = sum(1 for s, _ in pdf.values() if s == "verified")
        n_led = sum(1 for k in bib
                    if "pdf_verified" in ledger.get(k, {})
                    and pdf.get(k, ("missing",))[0] != "verified")
        n_mis = sum(1 for s, _ in pdf.values() if s == "mismatch")
        print(f"{len(bib)} entries; {len(errors)} errors, "
              f"{len(warnings)} warnings; PDFs: {n_ver} verified here, "
              f"{n_led} verified elsewhere (ledger), {n_mis} mismatched, "
              f"{sum(1 for s, _ in pdf.values() if s == 'missing')} missing "
              "locally")
        for e in errors:
            print(f"ERROR: {e}")
        print(f"report: {os.path.relpath(REPORT, os.getcwd())}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

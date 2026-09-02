r"""
Reference checking/annotation pipeline for the Diophantine Library.

Validates ``diophantine_classifier/data/references.bib`` against the family
registry and against locally downloaded papers, and writes a status report.

Because different people run this on different machines, the ledger
(``references/status.yaml``) is **monotone**: successes are appended and
never removed, so a missing local PDF reports as "missing here" while the
ledger still shows where and when it was verified.

Monotone history is not the same as a standing claim, though, so each record
carries a fingerprint of *what was verified*: the SHA-256 of the PDF and of
the entry's bibliographic metadata, or the exact DOI/URL string that
resolved.  A record counts as a verification of the entry **as it stands
today** only when those fingerprints still match.  Change a title, swap the
PDF, or edit a DOI and the old record stays in the ledger but is reported as
stale — never as a verification of the new content.

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
    recorded together with the hash of the file and of the metadata it was
    checked against.

online (``--online``):
  - DOIs must resolve at ``https://doi.org/`` and ``url`` fields must be
    reachable; successes are recorded together with the value that
    resolved.

Usage::

    sage -python tools/check_references.py [--online] [--quiet]

Writes ``references/REPORT.md``, updates ``references/status.yaml`` when new
verifications happen, and exits nonzero on hard errors (unknown keys,
malformed entries, empty annotations).
"""

import argparse
import datetime
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from dataclasses import dataclass

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from diophantine_classifier.references import (bibliography, _delatex)  # noqa: E402
from diophantine_classifier.registry import families  # noqa: E402

REQUIRED_FIELDS = ("author", "title", "year")
DOI_RE = re.compile(r"^10\.\d{4,9}/\S+$")
ARXIV_RE = re.compile(r"^(\d{4}\.\d{4,5}|[a-z-]+/\d{7})(v\d+)?$")

#: the kinds of verification the ledger records
LEDGER_KINDS = ("pdf", "doi", "url")


@dataclass(frozen=True)
class Paths:
    r"""
    Where the pipeline reads and writes.

    Passed around rather than hard-coded so a test can point the whole run
    at a temporary directory without touching committed files.
    """
    pdf_dir: str
    report: str
    ledger: str


def default_paths(root=ROOT):
    r"""
    The committed locations, under the repository root.
    """
    return Paths(pdf_dir=os.path.join(root, "references", "pdf"),
                 report=os.path.join(root, "references", "REPORT.md"),
                 ledger=os.path.join(root, "references", "status.yaml"))


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


def load_ledger(paths):
    r"""
    Read the committed verification ledger (empty dict if absent).
    """
    if os.path.exists(paths.ledger):
        with open(paths.ledger) as fobj:
            return yaml.safe_load(fobj) or {}
    return {}


def metadata_fingerprint(entry):
    r"""
    A stable hash of the bibliographic data a text check is made against.

    Normalized JSON over author, title, year and venue, so the fingerprint
    is insensitive to LaTeX markup and whitespace but changes the moment any
    of those fields is edited.
    """
    fields = entry["fields"]

    def norm(name):
        return " ".join(_delatex(fields.get(name, "")).split()).lower()

    payload = {"author": norm("author"), "title": norm("title"),
               "year": norm("year"),
               "venue": (norm("journal") or norm("booktitle")
                         or norm("publisher"))}
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def file_fingerprint(path):
    r"""
    SHA-256 of a file, or ``None`` when it is not there.
    """
    if not os.path.exists(path):
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as fobj:
        for chunk in iter(lambda: fobj.read(1 << 16), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(ledger, key, kind, fingerprint, today):
    r"""
    Append a verification of specific content to the ledger.

    ``fingerprint`` is the dict of fields identifying what was verified
    (``pdf_sha256``/``metadata_sha256``, or ``value``).  A record for the
    same content is not duplicated; records for *other* content are kept, so
    the history stays monotone while only matching records count as current.

    Returns ``True`` if the ledger changed.
    """
    if kind not in LEDGER_KINDS:
        raise ValueError(f"unknown ledger kind {kind!r}")
    records = ledger.setdefault(key, {}).setdefault(kind, [])
    if any(all(rec.get(k) == v for k, v in fingerprint.items())
           for rec in records):
        return False
    records.append(dict(fingerprint, verified_on=str(today)))
    return True


def _latest(records):
    r"""
    The most recent ``verified_on`` among some records.
    """
    return max((rec.get("verified_on", "") for rec in records), default="")


def verification(ledger, key, kind, fingerprint):
    r"""
    How the ledger stands against the entry's *current* content.

    Returns ``(status, date)`` with status ``"current"`` (a record matches
    the content as it is now), ``"stale"`` (records exist, but for other
    content), or ``"none"``.
    """
    records = ledger.get(key, {}).get(kind) or []
    if not records:
        return "none", None
    matching = [rec for rec in records
                if all(rec.get(k) == v for k, v in fingerprint.items())]
    if matching:
        return "current", _latest(matching)
    return "stale", _latest(records)


def pdf_verification(ledger, key, metadata_sha256, pdf_sha256):
    r"""
    The ledger's verdict on a PDF, given the current content.

    Adds one status to :func:`verification`: ``"unconfirmed"``, for a record
    made against exactly this metadata when the file itself is not on this
    machine to hash.  That is the monotone case — nobody's missing download
    invalidates someone else's verification — while an edited title or a
    replaced file both come back ``"stale"``.
    """
    records = ledger.get(key, {}).get("pdf") or []
    if not records:
        return "none", None
    same_metadata = [rec for rec in records
                     if rec.get("metadata_sha256") == metadata_sha256]
    if not same_metadata:
        return "stale", _latest(records)
    if pdf_sha256 is None:
        return "unconfirmed", _latest(same_metadata)
    exact = [rec for rec in same_metadata
             if rec.get("pdf_sha256") == pdf_sha256]
    if exact:
        return "current", _latest(exact)
    return "stale", _latest(records)


def save_ledger(ledger, paths):
    r"""
    Write the ledger with stable ordering (clean diffs).
    """
    os.makedirs(os.path.dirname(paths.ledger), exist_ok=True)
    with open(paths.ledger, "w") as fobj:
        fobj.write("# Verification ledger: successes recorded by "
                   "tools/check_references.py.\n"
                   "# Monotone: records are appended, never removed. Each "
                   "one fingerprints what was\n"
                   "# verified, so it stops counting as current when that "
                   "content changes.\n")
        yaml.safe_dump({k: {kind: list(recs)
                            for kind, recs in sorted(v.items())}
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


def pdf_status(key, entry, paths):
    r"""
    Check a locally downloaded PDF against its bibliography entry.

    Returns ``(status, note)`` with status one of ``missing``, ``present``
    (no text tool), ``suspect``, ``mismatch``, ``verified``.  Only
    ``verified`` enters the ledger.
    """
    path = os.path.join(paths.pdf_dir, f"{key}.pdf")
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


def ledger_verdicts(bib, ledger, paths):
    r"""
    For every entry, how the ledger stands against its current content.

    Returns ``{key: {"pdf": (status, date), "doi": ..., "url": ...}}``.
    This is the single place that decides "verified" versus "stale", so the
    report, the warnings and the tests all agree.
    """
    out = {}
    for key, entry in bib.items():
        fields = entry["fields"]
        pdf_sha = file_fingerprint(os.path.join(paths.pdf_dir, f"{key}.pdf"))
        verdict = {"pdf": pdf_verification(ledger, key,
                                           metadata_fingerprint(entry),
                                           pdf_sha)}
        for kind in ("doi", "url"):
            value = fields.get(kind)
            verdict[kind] = (("none", None) if value is None else
                             verification(ledger, key, kind,
                                          {"value": value}))
        out[key] = verdict
    return out


def _pdf_cell(status, note, ledger_status, ledger_date):
    r"""
    Render the local-PDF column, merging in the ledger's verdict.
    """
    if status == "verified":
        cell = f"verified ({note})"
    else:
        cell = status if status == "missing" else f"{status} ({note})"
    if ledger_status == "current":
        cell += f"; ledger: verified {ledger_date}"
    elif ledger_status == "unconfirmed":
        cell += (f"; ledger: verified {ledger_date} elsewhere (file not "
                 "here to hash)")
    elif ledger_status == "stale":
        cell += (f"; ledger: STALE — the {ledger_date} record was for "
                 "different content")
    return cell


def _link_cell(present, ledger_status, ledger_date):
    r"""
    Render a doi/url column against the ledger's verdict.
    """
    cell = "yes" if present else "TODO"
    if ledger_status == "current":
        cell += f" (ok {ledger_date})"
    elif ledger_status == "stale":
        cell += f" (STALE: {ledger_date} record was for another value)"
    return cell


def write_report(bib, used, pdf, online, ledger, errors, warnings, paths):
    r"""
    Write ``references/REPORT.md`` from this run merged with the ledger.
    """
    lines = ["# Reference status report", "",
             "Generated by `tools/check_references.py`. Download papers to "
             "`references/pdf/<key>.pdf` (legally!) and re-run to verify "
             "them against the bibliography. The ledger in `status.yaml` is "
             "monotone — records are appended, never removed — but each one "
             "fingerprints what was verified, so editing an entry or "
             "replacing a PDF reports the old record as stale rather than "
             "as a verification of the new content.", ""]
    verdicts = ledger_verdicts(bib, ledger, paths)
    n_doi = sum(1 for e in bib.values() if "doi" in e["fields"])
    n_url = sum(1 for e in bib.values() if "url" in e["fields"])
    n_ver = sum(1 for k in bib
                if pdf.get(k, ("missing",))[0] == "verified"
                or verdicts[k]["pdf"][0] in ("current", "unconfirmed"))
    n_stale = sum(1 for k in bib
                  if any(verdicts[k][kind][0] == "stale"
                         for kind in LEDGER_KINDS))
    lines += [f"- entries: {len(bib)}; used by families: {len(used)}",
              f"- with DOI: {n_doi}; with free URL: {n_url}; "
              f"PDFs verified (here or in ledger): {n_ver}",
              f"- entries whose ledger records are stale (content changed "
              f"since): {n_stale}", ""]
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
        verdict = verdicts[key]
        doi = _link_cell("doi" in fields, *verdict["doi"])
        url = _link_cell("url" in fields, *verdict["url"])
        status, note = pdf.get(key, ("missing", ""))
        cell = _pdf_cell(status, note, *verdict["pdf"])
        if online and key in online:
            cell += "; online: " + "; ".join(f"{l}: {s}"
                                             for l, s in online[key])
        lines.append(f"| {key} | {fams} | {doi} | {url} | {cell} |")
    todo_doi = sorted(k for k, e in bib.items() if "doi" not in e["fields"])
    todo_url = sorted(k for k, e in bib.items() if "url" not in e["fields"])
    todo_pdf = sorted(k for k in bib
                      if pdf.get(k, ("missing",))[0] == "missing"
                      and verdicts[k]["pdf"][0] not in ("current",
                                                        "unconfirmed"))
    lines += ["", "## TODO", "",
              f"- add DOIs ({len(todo_doi)}): " + ", ".join(todo_doi),
              f"- find legally free URLs ({len(todo_url)}): "
              + ", ".join(todo_url),
              f"- download and verify PDFs ({len(todo_pdf)}): "
              + ", ".join(todo_pdf),
              ""]
    os.makedirs(os.path.dirname(paths.report), exist_ok=True)
    with open(paths.report, "w") as fobj:
        fobj.write("\n".join(lines))


def main(argv=None, paths=None, today=None):
    r"""
    Run the pipeline; returns the process exit code.

    ``paths`` and ``today`` are injectable so a test can run the whole thing
    against a temporary directory on a fixed date.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    parser.add_argument("--online", action="store_true",
                        help="also check that DOIs and URLs resolve")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args(argv)
    paths = paths or default_paths()
    today = today or datetime.date.today().isoformat()

    bib = bibliography()
    used = usage_map()
    ledger = load_ledger(paths)
    errors, warnings = check_structure(bib, used)
    pdf = {key: pdf_status(key, entry, paths) for key, entry in bib.items()}

    changed = False
    for key, (status, _) in pdf.items():
        if status != "verified":
            continue
        changed |= record(
            ledger, key, "pdf",
            {"pdf_sha256": file_fingerprint(
                os.path.join(paths.pdf_dir, f"{key}.pdf")),
             "metadata_sha256": metadata_fingerprint(bib[key])},
            today)

    online = {}
    if args.online:
        for key, entry in bib.items():
            res = online_status(entry)
            if res:
                online[key] = res
                for label, status in res:
                    if status == "ok":
                        changed |= record(
                            ledger, key, label,
                            {"value": entry["fields"][label]}, today)
                    else:
                        warnings.append(f"{key} {label}: {status}")

    verdicts = ledger_verdicts(bib, ledger, paths)
    for key, verdict in sorted(verdicts.items()):
        for kind in LEDGER_KINDS:
            if verdict[kind][0] == "stale":
                warnings.append(
                    f"{key}: the {kind} verification of {verdict[kind][1]} "
                    "is stale — the entry no longer holds the content that "
                    "was checked")
        if pdf[key][0] == "mismatch" and verdict["pdf"][0] != "none":
            warnings.append(
                f"{key}: the local file does not match the entry, though "
                f"the ledger has a {verdict['pdf'][0]} record — your copy "
                "may be the wrong paper")

    if changed:
        save_ledger(ledger, paths)
    write_report(bib, used, pdf, online, ledger, errors, warnings, paths)

    if not args.quiet:
        n_ver = sum(1 for s, _ in pdf.values() if s == "verified")
        n_led = sum(1 for k in bib
                    if verdicts[k]["pdf"][0] == "unconfirmed")
        n_stale = sum(1 for k in bib
                      if any(verdicts[k][kind][0] == "stale"
                             for kind in LEDGER_KINDS))
        n_mis = sum(1 for s, _ in pdf.values() if s == "mismatch")
        print(f"{len(bib)} entries; {len(errors)} errors, "
              f"{len(warnings)} warnings; PDFs: {n_ver} verified here, "
              f"{n_led} verified elsewhere (ledger), {n_mis} mismatched, "
              f"{sum(1 for s, _ in pdf.values() if s == 'missing')} missing "
              f"locally; {n_stale} stale ledger record(s)")
        for e in errors:
            print(f"ERROR: {e}")
        print(f"report: {os.path.relpath(paths.report, os.getcwd())}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

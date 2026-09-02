"""Bibliography and the reference-checking pipeline."""

import importlib.util
import os

import pytest

from diophantine_classifier import families
from diophantine_classifier.references import (bibliography, check_key,
                                               format_reference)

TOOLS = os.path.join(os.path.dirname(__file__), "..", "tools",
                     "check_references.py")


def _load_tool():
    spec = importlib.util.spec_from_file_location("check_references", TOOLS)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_bibliography_parses():
    bib = bibliography()
    assert len(bib) > 80
    assert bib["Wiles1995"]["fields"]["journal"].startswith("Annals")


def test_every_family_reference_resolves():
    for fam in families().values():
        for ref in fam.references:
            assert check_key(ref["key"]), f"{fam.slug}: {ref['key']}"
            assert ref["why"].strip(), f"{fam.slug}: {ref['key']} needs a why"


def test_every_family_has_references():
    missing = [fam.slug for fam in families().values() if not fam.references]
    assert not missing, f"families without references: {missing}"


def test_formatting_all_used_keys():
    used = {ref["key"] for fam in families().values()
            for ref in fam.references}
    for key in used:
        line = format_reference(key)
        assert bibliography()[key]["fields"]["year"] in line, key


def test_free_urls_are_plausible():
    for key, entry in bibliography().items():
        url = entry["fields"].get("url")
        if url:
            assert url.startswith("https://"), key


def test_pipeline_structural(tmp_path):
    mod = _load_tool()
    bib = bibliography()
    used = mod.usage_map()
    errors, warnings = mod.check_structure(bib, used)
    assert errors == []


def test_pipeline_report(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    rc = mod.main(["--quiet"], paths=paths, today="2026-08-09")
    assert rc == 0
    text = open(paths.report).read()
    assert "TODO" in text and "Lenstra2002" in text


def _paths(mod, tmp_path):
    """A run confined to a temporary directory."""
    return mod.Paths(pdf_dir=str(tmp_path / "pdf"),
                     report=str(tmp_path / "REPORT.md"),
                     ledger=str(tmp_path / "status.yaml"))


def _fake_pdf(path, key):
    path.parent.mkdir(parents=True, exist_ok=True)
    (path.parent / f"{key}.pdf").write_bytes(
        b"%PDF-1.4\n" + b"unrelated content " * 1000)


def test_pdf_mismatch_detection(tmp_path):
    """A bogus local PDF must not verify."""
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    _fake_pdf(tmp_path / "pdf" / "x", "Wiles1995")
    status, note = mod.pdf_status("Wiles1995", bibliography()["Wiles1995"],
                                  paths)
    assert status in ("mismatch", "suspect", "present")
    assert status != "verified"


# --- verification is tied to the content that was verified (brief 4.2) ---

def _verified_run(mod, paths, key, today="2026-08-09"):
    """Run the pipeline with `key`'s PDF check forced to succeed."""
    real = mod.pdf_status

    def fake(k, entry, p):
        return ("verified", "simulated") if k == key else real(k, entry, p)

    mod.pdf_status = fake
    try:
        return mod.main(["--quiet"], paths=paths, today=today)
    finally:
        mod.pdf_status = real


def test_matching_pdf_and_metadata_is_a_current_verification(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    _fake_pdf(tmp_path / "pdf" / "x", "Lenstra2002")
    assert _verified_run(mod, paths, "Lenstra2002") == 0
    ledger = mod.load_ledger(paths)
    record, = ledger["Lenstra2002"]["pdf"]
    assert record["verified_on"] == "2026-08-09"
    assert record["pdf_sha256"] and record["metadata_sha256"]
    verdict = mod.ledger_verdicts(bibliography(), ledger, paths)
    assert verdict["Lenstra2002"]["pdf"] == ("current", "2026-08-09")


def test_changing_the_title_makes_the_verification_stale(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    _fake_pdf(tmp_path / "pdf" / "x", "Lenstra2002")
    _verified_run(mod, paths, "Lenstra2002")
    ledger = mod.load_ledger(paths)

    bib = bibliography()
    edited = {k: {**v, "fields": dict(v["fields"])} for k, v in bib.items()}
    edited["Lenstra2002"]["fields"]["title"] = "A completely different paper"
    verdict = mod.ledger_verdicts(edited, ledger, paths)
    assert verdict["Lenstra2002"]["pdf"][0] == "stale"


def test_replacing_the_pdf_makes_the_verification_stale(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    _fake_pdf(tmp_path / "pdf" / "x", "Lenstra2002")
    _verified_run(mod, paths, "Lenstra2002")
    ledger = mod.load_ledger(paths)
    assert mod.ledger_verdicts(bibliography(), ledger,
                               paths)["Lenstra2002"]["pdf"][0] == "current"

    (tmp_path / "pdf" / "Lenstra2002.pdf").write_bytes(
        b"%PDF-1.4\n" + b"a different file entirely " * 1000)
    verdict = mod.ledger_verdicts(bibliography(), ledger, paths)
    assert verdict["Lenstra2002"]["pdf"][0] == "stale"


def test_a_missing_local_pdf_does_not_downgrade_the_ledger(tmp_path):
    """Monotone: someone else's verification survives a bare checkout."""
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    _fake_pdf(tmp_path / "pdf" / "x", "Lenstra2002")
    _verified_run(mod, paths, "Lenstra2002")
    ledger = mod.load_ledger(paths)

    (tmp_path / "pdf" / "Lenstra2002.pdf").unlink()
    verdict = mod.ledger_verdicts(bibliography(), ledger, paths)
    assert verdict["Lenstra2002"]["pdf"] == ("unconfirmed", "2026-08-09")


def test_changing_a_doi_or_url_makes_that_check_stale(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    bib = bibliography()
    key = next(k for k, e in bib.items() if "doi" in e["fields"])
    ledger = {}
    mod.record(ledger, key, "doi", {"value": bib[key]["fields"]["doi"]},
               "2026-08-09")
    assert mod.ledger_verdicts(bib, ledger, paths)[key]["doi"] == (
        "current", "2026-08-09")

    edited = {k: {**v, "fields": dict(v["fields"])} for k, v in bib.items()}
    edited[key]["fields"]["doi"] = "10.1234/somewhere-else"
    assert mod.ledger_verdicts(edited, ledger, paths)[key]["doi"][0] == "stale"


def test_history_is_kept_but_only_matching_records_are_current(tmp_path):
    """An old record stays in the ledger; it just stops counting."""
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    bib = bibliography()
    key = next(k for k, e in bib.items() if "url" in e["fields"])
    ledger = {}
    mod.record(ledger, key, "url", {"value": "https://old.example/paper"},
               "2020-01-01")
    mod.record(ledger, key, "url", {"value": bib[key]["fields"]["url"]},
               "2026-08-09")
    assert len(ledger[key]["url"]) == 2                # history preserved
    assert mod.ledger_verdicts(bib, ledger, paths)[key]["url"] == (
        "current", "2026-08-09")

    stale_only = {key: {"url": [ledger[key]["url"][0]]}}
    status, date = mod.ledger_verdicts(bib, stale_only, paths)[key]["url"]
    assert (status, date) == ("stale", "2020-01-01")


def test_record_does_not_duplicate_the_same_content(tmp_path):
    mod = _load_tool()
    ledger = {}
    fingerprint = {"value": "https://example.org/paper"}
    assert mod.record(ledger, "Key", "url", fingerprint, "2026-08-09")
    assert not mod.record(ledger, "Key", "url", fingerprint, "2026-08-10")
    assert len(ledger["Key"]["url"]) == 1


def test_stale_records_are_reported_as_warnings(tmp_path):
    mod = _load_tool()
    paths = _paths(mod, tmp_path)
    bib = bibliography()
    key = next(k for k, e in bib.items() if "doi" in e["fields"])
    os.makedirs(os.path.dirname(paths.ledger), exist_ok=True)
    with open(paths.ledger, "w") as fobj:
        fobj.write(f"{key}:\n  doi:\n  - value: 10.9999/not-this-one\n"
                   "    verified_on: '2020-01-01'\n")
    rc = mod.main([], paths=paths, today="2026-08-09")
    assert rc == 0                      # stale is a warning, not an error
    report = open(paths.report).read()
    assert "STALE" in report


def test_metadata_fingerprint_ignores_markup_but_not_content():
    mod = _load_tool()
    entry = {"fields": {"author": "A. Wiles", "title": "Modular  elliptic",
                        "year": "1995", "journal": "Annals"}}
    same = {"fields": {"author": "A. Wiles", "title": "Modular elliptic",
                       "year": "1995", "journal": "Annals"}}
    other = {"fields": {"author": "A. Wiles", "title": "Modular elliptic",
                        "year": "1996", "journal": "Annals"}}
    assert mod.metadata_fingerprint(entry) == mod.metadata_fingerprint(same)
    assert mod.metadata_fingerprint(entry) != mod.metadata_fingerprint(other)

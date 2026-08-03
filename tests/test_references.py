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


def test_pipeline_report(tmp_path, monkeypatch):
    mod = _load_tool()
    monkeypatch.setattr(mod, "REPORT", str(tmp_path / "REPORT.md"))
    rc = mod.main(["--quiet"])
    assert rc == 0
    text = (tmp_path / "REPORT.md").read_text()
    assert "TODO" in text and "Lenstra2002" in text


def test_pdf_mismatch_detection(tmp_path, monkeypatch):
    """A bogus local PDF must not verify."""
    mod = _load_tool()
    monkeypatch.setattr(mod, "PDF_DIR", str(tmp_path))
    fake = tmp_path / "Wiles1995.pdf"
    fake.write_bytes(b"%PDF-1.4\n" + b"unrelated content " * 1000)
    status, note = mod.pdf_status("Wiles1995", bibliography()["Wiles1995"])
    assert status in ("mismatch", "suspect", "present")
    assert status != "verified"


def test_ledger_never_downgrades(tmp_path, monkeypatch):
    """A verification recorded in the ledger survives a run on a machine
    without the file (status is monotone)."""
    mod = _load_tool()
    monkeypatch.setattr(mod, "PDF_DIR", str(tmp_path / "pdf"))  # empty
    monkeypatch.setattr(mod, "REPORT", str(tmp_path / "REPORT.md"))
    monkeypatch.setattr(mod, "LEDGER", str(tmp_path / "status.yaml"))
    (tmp_path / "status.yaml").write_text(
        "Wiles1995:\n  pdf_verified: '2026-01-01'\n")
    rc = mod.main(["--quiet"])
    assert rc == 0
    # ledger untouched (no changes -> not rewritten, content preserved)
    assert "2026-01-01" in (tmp_path / "status.yaml").read_text()
    report = (tmp_path / "REPORT.md").read_text()
    assert "ledger: verified 2026-01-01" in report
    # the ledger entry keeps Wiles1995 off the download TODO list
    todo = [l for l in report.splitlines() if "download and verify" in l]
    assert todo
    names = todo[0].split(": ", 1)[1].split(", ")
    assert "Wiles1995" not in names and "TaylorWiles1995" in names


def test_ledger_records_upgrade(tmp_path, monkeypatch):
    """A local verification is written to the ledger, write-once."""
    mod = _load_tool()
    monkeypatch.setattr(mod, "PDF_DIR", str(tmp_path / "pdf"))
    monkeypatch.setattr(mod, "REPORT", str(tmp_path / "REPORT.md"))
    monkeypatch.setattr(mod, "LEDGER", str(tmp_path / "status.yaml"))
    real_pdf_status = mod.pdf_status

    def fake_status(key, entry):
        if key == "Lenstra2002":
            return "verified", "simulated"
        return real_pdf_status(key, entry)

    monkeypatch.setattr(mod, "pdf_status", fake_status)
    rc = mod.main(["--quiet"])
    assert rc == 0
    ledger = (tmp_path / "status.yaml").read_text()
    assert "Lenstra2002" in ledger and "pdf_verified" in ledger


def test_ledger_mismatch_keeps_verification(tmp_path, monkeypatch):
    """A local mismatch warns but never erases a recorded verification."""
    mod = _load_tool()
    monkeypatch.setattr(mod, "PDF_DIR", str(tmp_path))
    monkeypatch.setattr(mod, "REPORT", str(tmp_path / "REPORT.md"))
    monkeypatch.setattr(mod, "LEDGER", str(tmp_path / "status.yaml"))
    (tmp_path / "status.yaml").write_text(
        "Wiles1995:\n  pdf_verified: '2026-01-01'\n")
    (tmp_path / "Wiles1995.pdf").write_bytes(
        b"%PDF-1.4\n" + b"unrelated content " * 1000)
    rc = mod.main(["--quiet"])
    assert rc == 0                      # a mismatch is a warning, not an error
    assert "2026-01-01" in (tmp_path / "status.yaml").read_text()

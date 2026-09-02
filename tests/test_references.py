"""Bibliography and the reference-checking pipeline."""

from diophantine_classifier.references import bibliography


def test_bibliography_parses():
    bib = bibliography()
    assert len(bib) > 80
    assert bib["Wiles1995"]["fields"]["journal"].startswith("Annals")


def test_free_urls_are_plausible():
    for key, entry in bibliography().items():
        url = entry["fields"].get("url")
        if url:
            assert url.startswith("https://"), key

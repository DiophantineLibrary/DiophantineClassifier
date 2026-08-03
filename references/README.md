# Local reference library

Companion to `diophantine_classifier/data/references.bib`.

- Download papers you have legal access to as `pdf/<bibtex-key>.pdf`
  (e.g. `pdf/Lenstra2002.pdf`). PDFs are **not** committed to the repo
  (see `.gitignore`); this directory holds the verification report and
  conventions.
- Run the pipeline to validate the bibliography and check each downloaded
  PDF against its entry (title words + author surnames must appear in the
  extracted text):

  ```bash
  sage -python tools/check_references.py           # structural + local PDFs
  sage -python tools/check_references.py --online  # also resolve DOIs/URLs
  ```

- The result is written to [REPORT.md](REPORT.md), including TODO lists of
  entries missing DOIs, missing legally-free URLs, and missing local PDFs.

Policy: the `url` field in the bibliography must point to a *legally free*
copy only — arXiv, an open journal archive (AMS Notices, Math. Comp. open
archive), the author's own page, or a public-domain scan. When in doubt,
leave it out and let the report track it as a TODO.

Text extraction uses `pdftotext` (poppler) when installed:
`brew install poppler`.

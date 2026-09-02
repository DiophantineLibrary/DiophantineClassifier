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
- Verification history is **monotone** across machines: successes (PDF
  verified, DOI/URL resolved) are appended to [status.yaml](status.yaml),
  which is committed, and never removed. A run on a machine without a given
  PDF reports it as "missing locally" but never erases the recorded
  verification; a *local* mismatch warns (your copy may be the wrong file)
  without downgrading the ledger.
- Each record fingerprints **what** was verified: the SHA-256 of the PDF
  and of the entry's author/title/year/venue, or the exact DOI or URL
  string that resolved. A record counts as a verification of the entry *as
  it stands today* only while those fingerprints still match, so editing a
  title, replacing a PDF or changing a DOI reports the old record as stale
  instead of quietly carrying its success over to different content. Only
  entries with a *current* record stay off the download TODO list.

Policy: the `url` field in the bibliography must point to a *legally free*
copy only — arXiv, an open journal archive (AMS Notices, Math. Comp. open
archive), the author's own page, or a public-domain scan. When in doubt,
leave it out and let the report track it as a TODO.

Text extraction uses `pdftotext` (poppler) when installed:
`brew install poppler`.

r"""
Bibliography support: parse ``data/references.bib`` and format entries.

Families refer to references by BibTeX key (see ``data/families/<slug>.yaml``),
each use carrying a ``why`` annotation.  This module provides the parser (a
small, tolerant scanner sufficient for our own controlled file — not a general
BibTeX implementation) and a plain-text formatter used by
:meth:`~diophantine_classifier.classify.Classification.explain`.

EXAMPLES::

    sage: from diophantine_classifier.references import bibliography, format_reference
    sage: entry = bibliography()["Lenstra2002"]
    sage: entry["fields"]["year"]
    '2002'
    sage: fmt = format_reference("Lenstra2002")
    sage: "Solving the Pell equation" in fmt and "182-192" in fmt
    True
"""

import os
import re
import unicodedata
from functools import lru_cache

BIB_PATH = os.path.join(os.path.dirname(__file__), "data", "references.bib")

#: LaTeX accent commands -> Unicode combining characters
_ACCENTS = {
    '"': "\u0308", "'": "\u0301", "`": "\u0300", "^": "\u0302",
    "~": "\u0303", "H": "\u030b", "u": "\u0306", "v": "\u030c",
    "c": "\u0327", "=": "\u0304", ".": "\u0307",
}


class BibError(ValueError):
    """
    Raised when ``references.bib`` cannot be parsed.
    """


def _skip_whitespace(text, i):
    r"""
    Return the first index ``>= i`` of a non-whitespace character.

    INPUT:

    - ``text`` -- string
    - ``i`` -- starting index

    OUTPUT: an integer index (``i`` itself if it is already past the end)

    EXAMPLES::

        sage: from diophantine_classifier.references import _skip_whitespace
        sage: _skip_whitespace("a   b", 1)
        4
        sage: _skip_whitespace("ab", 5)
        5
    """
    while i < len(text) and text[i].isspace():
        i += 1
    return i


def _read_braced(text, i):
    r"""
    Read a ``{...}``-balanced group starting at index ``i``.

    INPUT:

    - ``text`` -- string with ``text[i] == "{"``
    - ``i`` -- index of the opening brace

    OUTPUT: pair ``(contents, j)`` where ``contents`` excludes the outer
    braces and ``j`` is the index just past the closing brace

    EXAMPLES::

        sage: from diophantine_classifier.references import _read_braced
        sage: _read_braced("{a{b}c} rest", 0)
        ('a{b}c', 7)
    """
    depth = 0
    j = i
    while j < len(text):
        if text[j] == "{":
            depth += 1
        elif text[j] == "}":
            depth -= 1
            if depth == 0:
                return text[i + 1:j], j + 1
        j += 1
    raise BibError(f"unbalanced braces starting at index {i}")


def _parse_fields(body):
    r"""
    Parse the ``name = value`` fields of one BibTeX entry body.

    INPUT:

    - ``body`` -- the text between the entry's outer braces, *after* the key
      and its trailing comma

    OUTPUT: a dict mapping lowercased field names to raw string values

    EXAMPLES::

        sage: from diophantine_classifier.references import _parse_fields
        sage: _parse_fields('author = {A. B.},\n year = {1999}')
        {'author': 'A. B.', 'year': '1999'}
        sage: _parse_fields('pages = "1--5"')
        {'pages': '1--5'}
    """
    fields = {}
    i = 0
    while True:
        i = _skip_whitespace(body, i)
        if i >= len(body):
            break
        m = re.match(r"([A-Za-z]+)\s*=\s*", body[i:])
        if not m:
            raise BibError(f"expected 'name =' at ...{body[i:i + 30]!r}")
        name = m.group(1).lower()
        i += m.end()
        if i < len(body) and body[i] == "{":
            value, i = _read_braced(body, i)
        elif i < len(body) and body[i] == '"':
            j = body.index('"', i + 1)
            value, i = body[i + 1:j], j + 1
        else:
            m2 = re.match(r"[^,\s]+", body[i:])
            value = m2.group(0)
            i += m2.end()
        fields[name] = re.sub(r"\s+", " ", value.strip())
        i = _skip_whitespace(body, i)
        if i < len(body) and body[i] == ",":
            i += 1
    return fields


@lru_cache(maxsize=None)
def bibliography(path=None):
    r"""
    Parse the bibliography file into a dict keyed by BibTeX key.

    INPUT:

    - ``path`` -- (default: the package's ``data/references.bib``) file to parse

    OUTPUT: dict mapping key to ``{"type": ..., "fields": {...}}``

    EXAMPLES::

        sage: from diophantine_classifier.references import bibliography
        sage: bib = bibliography()
        sage: bib["Wiles1995"]["type"]
        'article'
        sage: bib["Wiles1995"]["fields"]["volume"]
        '141'
        sage: len(bib) > 80
        True
    """
    path = path or BIB_PATH
    text = "\n".join(line for line in open(path).read().splitlines()
                     if not line.lstrip().startswith("%"))
    entries = {}
    i = 0
    while True:
        at = text.find("@", i)
        if at < 0:
            break
        m = re.match(r"@(\w+)\s*", text[at:])
        if not m:
            raise BibError(f"malformed entry at index {at}")
        etype = m.group(1).lower()
        j = at + m.end()
        if j >= len(text) or text[j] != "{":
            raise BibError(f"expected '{{' after @{etype}")
        body, i = _read_braced(text, j)
        key, _, rest = body.partition(",")
        key = key.strip()
        if not key:
            raise BibError(f"missing key in @{etype} entry")
        if key in entries:
            raise BibError(f"duplicate key {key!r}")
        entries[key] = {"type": etype, "fields": _parse_fields(rest)}
    return entries


def _delatex(s):
    r"""
    Convert the LaTeX markup used in ``references.bib`` to plain text.

    Handles the accent commands appearing in the file (``{\"u}``, ``{\'e}``,
    ``{\H o}``, ...), strips braces, dollar signs and remaining backslashes,
    and rewrites ``--`` as ``-``.

    EXAMPLES::

        sage: from diophantine_classifier.references import _delatex
        sage: _delatex(r'Journal f{\"u}r die reine und angewandte Mathematik')
        'Journal für die reine und angewandte Mathematik'
        sage: _delatex(r'{$x^2 + 7 = 2^n$}')
        'x^2 + 7 = 2^n'
        sage: _delatex('1--5')
        '1-5'
    """
    def accent(m):
        mark, letter = m.group(1), m.group(2)
        return unicodedata.normalize("NFC", letter + _ACCENTS[mark])

    s = re.sub(r"\\([\"'`^~Huvc=.])\s*\{?([A-Za-z])\}?", accent, s)
    s = s.replace("--", "-")
    s = re.sub(r"[{}$]", "", s)
    s = s.replace("\\", "")
    return re.sub(r"\s+", " ", s).strip()


def format_reference(key, bib=None):
    r"""
    Format one bibliography entry as a single plain-text line.

    INPUT:

    - ``key`` -- BibTeX key present in the bibliography
    - ``bib`` -- (default: the package bibliography) parsed dict, as returned
      by :func:`bibliography`

    OUTPUT: string of the shape ``authors, title, venue (year), pages.
    doi:... [free url]`` — pieces are omitted when absent

    EXAMPLES::

        sage: from diophantine_classifier.references import format_reference
        sage: fmt = format_reference("Thue1909")
        sage: fmt.startswith('Axel Thue, Über Annäherungswerte')
        True
        sage: 'doi:10.1515/crll.1909.135.284' in fmt
        True
        sage: format_reference("nope")
        Traceback (most recent call last):
        ...
        KeyError: 'nope'
    """
    if bib is None:
        bib = bibliography()
    entry = bib[key]
    f = entry["fields"]
    authors = _delatex(f.get("author", "")).replace(" and ", ", ")
    title = _delatex(f.get("title", ""))
    parts = [p for p in (authors, title) if p]
    venue = ""
    if "journal" in f:
        venue = _delatex(f["journal"])
        if "volume" in f:
            venue += f" {f['volume']}"
    elif "booktitle" in f:
        venue = "in " + _delatex(f["booktitle"])
    elif "publisher" in f:
        venue = _delatex(f["publisher"])
    elif "note" in f:
        venue = _delatex(f["note"])
    if venue:
        if "year" in f:
            venue += f" ({f['year']})"
        parts.append(venue)
    elif "year" in f:
        parts.append(f"({f['year']})")
    if "pages" in f and "journal" in f:
        parts.append(_delatex(f["pages"]))
    line = ", ".join(parts) + "."
    if "doi" in f:
        line += f" doi:{f['doi']}"
    if "url" in f:
        line += f" [{f['url']}]"
    return line


def check_key(key):
    r"""
    Return ``True`` if ``key`` exists in the package bibliography.

    EXAMPLES::

        sage: from diophantine_classifier.references import check_key
        sage: check_key("Faltings1983")
        True
        sage: check_key("Fermat1637")
        False
    """
    return key in bibliography()

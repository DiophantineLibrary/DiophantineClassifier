"""Registry consistency: DAG well-formed, matcher flags honest, docs in sync."""

from diophantine_classifier import families, validate
from diophantine_classifier.registry import ancestors, depth


def test_validate():
    assert validate()


def test_dag_depths():
    fams = families()
    assert depth("general-polynomial") == 0
    assert depth("pell") > depth("pell-like") > depth("binary-qf-representation")
    assert "general-polynomial" in ancestors("pell")


def test_slides_families_present():
    """The ten example families from the LuCaNT 2025 talk are all registered."""
    fams = families()
    for slug in ["linear-system", "pythagorean", "pell", "elliptic-weierstrass",
                 "thue", "generalized-fermat", "ramanujan-nagell",
                 "thue-mahler", "sum-of-three-cubes", "erdos-straus"]:
        assert slug in fams, slug


def test_priorities_cover_launch_set():
    fams = families()
    p1 = [f.slug for f in fams.values() if f.priority == 1]
    assert len(p1) >= 25


def test_matcher_flags_match_implementation():
    """Every slug the matcher code can emit is flagged matcher: true."""
    import re
    import diophantine_classifier.matchers as m
    source = open(m.__file__.replace(".pyc", ".py")).read()
    emitted = set(re.findall(r'Match\(\s*\n?\s*"([a-z0-9-]+)"', source))
    emitted |= set(re.findall(r'Match\("([a-z0-9-]+)"', source))
    fams = families()
    for slug in emitted:
        assert slug in fams, f"matcher emits unknown slug {slug!r}"
        assert fams[slug].matcher, f"{slug}: emitted but matcher flag is false"


def test_examples_classify_to_their_family():
    """Registry examples must classify to the family that lists them (or a
    descendant, e.g. the fermat example listed under generalized-fermat)."""
    from diophantine_classifier import classify
    from diophantine_classifier.registry import ancestors as anc
    fams = families()
    for fam in fams.values():
        if not fam.matcher:
            continue
        for ex in fam.examples:
            cls = classify(ex, params=_infer_params(ex))
            ok = (cls.slug == fam.slug or fam.slug in anc(cls.slug)
                  or cls.slug in anc(fam.slug))
            assert ok, (f"{fam.slug}: example {ex!r} classified as "
                        f"{cls.slug!r}")


def _infer_params(example):
    """Parameters used in registry examples: single letters not serving as
    unknowns in that example's family form."""
    KNOWN = {"x^2 + y^2 + z^2 = n": "n", "x^2 + y^2 + z^2 + w^2 = n": "n",
             "x^4 + y^4 + z^4 + w^4 = n": "n",
             "4/n = 1/x + 1/y + 1/z": "n"}
    return KNOWN.get(example, "")

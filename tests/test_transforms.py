"""Coordinate maps: every normalization a matcher performs must be
executable and invertible, not merely described in prose.

A solver works in a family's standard coordinates; the transform is what
carries its answer back to the variables the user wrote.
"""

import json

import pytest

from diophantine_classifier.matchers import run
from diophantine_classifier.parsing import parse
from diophantine_classifier.transforms import (CoordinateTransform, identity,
                                               negate, rename)


def match_for(equation, slug, **kwargs):
    for m in run(parse(equation, **kwargs)):
        if m.slug == slug:
            return m
    raise AssertionError(f"{equation!r} did not match {slug!r}")


# --- the map itself -------------------------------------------------------

def test_identity_round_trips_and_knows_it_is_one():
    t = identity(("x", "y"))
    assert t.is_identity
    assert t.push_forward({"x": 1, "y": 2}) == {"x": 1, "y": 2}
    assert t.pull_back({"x": 1, "y": 2}) == {"x": 1, "y": 2}


def test_equation_level_operations_are_not_coordinate_changes():
    """Multiplying an equation by -1 changes no coordinate."""
    t = identity(("x", "y"), operations=("multiplied by -1",))
    assert t.is_identity
    assert bool(t)                      # something did happen, though
    assert str(t) == "multiplied by -1"
    assert t.pull_back({"x": 3, "y": 4}) == {"x": 3, "y": 4}


def test_rename_is_invertible():
    t = rename([("x", "b"), ("y", "a")])
    assert t.push_forward({"a": 2, "b": 3}) == {"x": 3, "y": 2}
    assert t.pull_back(t.push_forward({"a": 2, "b": 3})) == {"a": 2, "b": 3}


def test_negate_is_its_own_inverse():
    t = negate(("x", "y", "z"), ("z",))
    source = {"x": 1, "y": 2, "z": 3}
    assert t.push_forward(source) == {"x": 1, "y": 2, "z": -3}
    assert t.pull_back(t.push_forward(source)) == source


def test_missing_coordinate_is_an_error_not_a_guess():
    t = rename([("x", "a"), ("y", "b")])
    with pytest.raises(KeyError):
        t.push_forward({"a": 1})


def test_rational_coordinates_survive_transport():
    from sage.all import QQ
    t = negate(("x",), ("x",))
    assert t.pull_back({"x": QQ(1) / 2}) == {"x": QQ(-1) / 2}


def test_as_dict_is_json_serializable():
    t = rename([("x", "b"), ("y", "a")], description="swapped variables",
               operations=("multiplied by -1",))
    blob = json.dumps(t.as_dict())
    assert "swapped variables" in blob
    assert json.loads(blob)["to_source"] == {"a": "y", "b": "x"}


def test_default_transform_is_the_identity():
    assert CoordinateTransform().is_identity


# --- composition (brief 5.1) ---------------------------------------------

def test_transform_composition_round_trip():
    first = negate(("u", "v"), ("v",))            # source -> intermediate
    second = rename([("x", "v"), ("y", "u")])     # intermediate -> standard
    combined = first.then(second)
    source = {"u": 2, "v": 3}
    assert combined.pull_back(combined.push_forward(source)) == source


def test_composition_agrees_with_applying_the_steps_in_turn():
    first = negate(("u", "v"), ("v",))
    second = rename([("x", "v"), ("y", "u")])
    combined = first.then(second)
    source = {"u": 5, "v": -7}
    assert combined.push_forward(source) == second.push_forward(
        first.push_forward(source))
    standard = {"x": 11, "y": 13}
    assert combined.pull_back(standard) == first.pull_back(
        second.pull_back(standard))


def test_composition_tracks_the_variable_lists():
    first = rename([("a", "u"), ("b", "v")])
    second = rename([("x", "a"), ("y", "b")])
    combined = first.then(second)
    assert combined.source_variables == ("u", "v")
    assert combined.normalized_variables == ("x", "y")


def test_composition_with_identity_changes_nothing():
    t = rename([("x", "v"), ("y", "u")])
    assert identity(("u", "v")).then(t).to_normalized == t.to_normalized
    assert t.then(identity(("x", "y"))).to_source == t.to_source


def test_composition_keeps_operations_in_order():
    first = identity(("u",), operations=("multiplied by -1",))
    second = identity(("u",), operations=("divided by the content",))
    assert first.then(second).operations == ("multiplied by -1",
                                             "divided by the content")


def test_composition_transports_conditions():
    """A condition of the second map is stated in *its* source coordinates;
    composing must pull it back so the whole map's conditions are checkable
    on the original assignment."""
    from diophantine_classifier.conditions import NonzeroCondition
    first = negate(("u", "v"), ("v",))
    second = rename([("x", "v"), ("y", "u")],
                    conditions=(NonzeroCondition("v", ("v",)),))
    combined = first.then(second)
    assert [str(c) for c in combined.conditions] == ["-v != 0"]
    assert combined.conditions_hold({"u": 1, "v": 3}) is True
    assert combined.conditions_hold({"u": 1, "v": 0}) is False


def test_composed_conditions_serialize():
    from diophantine_classifier.conditions import NonzeroCondition
    first = rename([("a", "u"), ("b", "v")],
                   conditions=(NonzeroCondition("u", ("u",)),))
    second = rename([("x", "a"), ("y", "b")],
                    conditions=(NonzeroCondition("b", ("b",)),))
    blob = json.dumps(first.then(second).as_dict())
    assert "conditions" in blob


# --- generalized Fermat: data and map agree (brief 5.2) ------------------

def test_generalized_fermat_sorting_updates_the_transform():
    m = match_for("u^5 + v^3 = w^7", "generalized-fermat")
    assert m.data["signature"] == "(3, 5, 7)"
    assert m.transform.to_normalized["x"] == "v"
    assert m.transform.to_normalized["y"] == "u"
    assert m.transform.to_normalized["z"] == "w"


def test_generalized_fermat_coefficients_follow_the_same_variables():
    """a is the coefficient of standard x, which is the user's v."""
    m = match_for("2*u^5 + 3*v^3 = 5*w^7", "generalized-fermat")
    assert m.data["signature"] == "(3, 5, 7)"
    assert (m.data["a"], m.data["b"], m.data["c"]) == ("3", "2", "5")
    assert m.transform.to_normalized["x"] == "v"


def test_generalized_fermat_sign_flip_and_permutation_compose():
    m = match_for("u^5 + v^3 + w^7 = 0", "generalized-fermat")
    source = {"u": 2, "v": 3, "w": 4}
    assert m.transform.pull_back(m.transform.push_forward(source)) == source
    negated = [expr for expr in m.transform.to_normalized.values()
               if expr.startswith("-")]
    assert len(negated) == 1


def test_global_orientation_records_equation_operation():
    m = match_for("-x^3 - y^3 = -z^3", "generalized-fermat")
    assert m.transform.operations == ("multiplied by -1",)


def test_symbolic_signature_does_not_claim_a_false_bijection():
    """x^p + y^p = z^q reuses an exponent, so there is no invertible rename
    onto independent canonical p, q, r; the structure is roles instead."""
    m = match_for("x^p + y^p = z^q", "generalized-fermat")
    assert m.transform.is_identity
    assert m.transform.roles["exponents"] == ["p", "p", "q"]
    assert m.transform.roles["bases"] == ["x", "y", "z"]


# --- the maps the matchers actually emit ---------------------------------

def test_every_match_carries_a_map_over_the_real_unknowns():
    for equation in ["3*x + 5*y = 1", "x^2 - 61*y^2 = 1", "y^2 = x^3 - 2",
                     "x^2 + y^2 = z^2", "x^2 + 7 = 2^n", "x^p - y^q = 1",
                     "x^3 + y^3 + z^3 = 0"]:
        pe = parse(equation)
        for m in run(pe):
            assert set(m.transform.source_variables) <= set(pe.unknowns)
            assert set(m.transform.to_source) <= set(pe.unknowns)


def test_pell_swap_transform_round_trip():
    """5*x^2 - y^2 = 1 is Pell only after reading the user's y as x."""
    m = match_for("5*x^2 - y^2 = 1", "pell")
    assert not m.transform.is_identity
    assert m.transform.to_normalized == {"x": "y", "y": "x"}
    source = {"x": 2, "y": 3}
    assert m.transform.pull_back(m.transform.push_forward(source)) == source


def test_pell_without_swap_is_the_identity():
    m = match_for("x^2 - 61*y^2 = 1", "pell")
    assert m.transform.is_identity
    assert m.transform.source_variables == ("x", "y")


def test_sign_flip_multiplication_is_an_operation_not_a_swap():
    """-x^2 + 5*y^2 = -1 is oriented by multiplying through by -1."""
    m = match_for("-x^2 + 5*y^2 = -1", "pell")
    assert m.transform.operations == ("multiplied by -1",)
    assert m.transform.is_identity


def test_odd_exponent_sign_flip_round_trip():
    """x^3 + y^3 + z^3 = 0 becomes generalized Fermat by a sign flip, then
    a permutation onto the canonical (x, y, z); the two compose into one
    invertible map."""
    m = match_for("x^3 + y^3 + z^3 = 0", "generalized-fermat")
    assert not m.transform.is_identity
    negated = [expr for expr in m.transform.to_normalized.values()
               if expr.startswith("-")]
    assert len(negated) == 1
    source = {"x": 3, "y": 4, "z": 5}
    assert m.transform.pull_back(m.transform.push_forward(source)) == source


def test_curve_role_reversal_is_executable():
    """x^3 = y^2 - 2 is Mordell with the user's y as the standard x."""
    m = match_for("x^3 = y^2 - 2", "mordell")
    assert m.transform.to_normalized == {"x": "x", "y": "y"}
    m2 = match_for("y^2 = x^3 - 2", "mordell")
    assert m2.transform.pull_back({"x": 3, "y": 5}) == {"x": 3, "y": 5}


def test_curve_roles_follow_the_user_names():
    """u^2 = v^3 - 2: the standard x is the user's v."""
    m = match_for("u^2 = v^3 - 2", "mordell")
    assert m.transform.to_normalized == {"x": "v", "y": "u"}
    assert m.transform.pull_back({"x": 3, "y": 5}) == {"v": 3, "u": 5}


def test_roles_live_on_the_transform_not_in_the_data():
    m = match_for("x^2 + y^2 = z^2", "pythagorean")
    assert m.transform.roles["hypotenuse"] == "z"
    assert sorted(m.transform.roles["legs"]) == ["x", "y"]
    assert "roles" not in m.data


def test_catalan_roles_are_a_coordinate_map():
    m = match_for("x^p - y^q = 1", "catalan")
    assert m.transform.pull_back({"x": 3, "p": 2, "y": 2, "q": 3}) == {
        "x": 3, "p": 2, "y": 2, "q": 3}


def test_catalan_orientation_is_invertible():
    """1 = y^q - x^p is the same equation with the roles the other way."""
    m = match_for("y^q - x^p = 1", "catalan")
    pulled = m.transform.pull_back({"x": 3, "p": 2, "y": 2, "q": 3})
    assert pulled == {"y": 3, "q": 2, "x": 2, "p": 3}


def test_ramanujan_nagell_roles_are_a_coordinate_map():
    m = match_for("u^2 + 7 = 2^k", "ramanujan-nagell")
    assert m.transform.pull_back({"x": 181, "n": 15}) == {"u": 181, "k": 15}

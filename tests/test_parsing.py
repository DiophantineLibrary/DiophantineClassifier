"""The parser: side conditions, the source problem, and exact validation."""

import pytest

from diophantine_classifier import parse
from diophantine_classifier.parsing import ParseError


def test_denominator_expression_is_preserved():
    """A condition names the denominator, not just its variables."""
    pe = parse("1/(x - 1) = 1/(y - 1)")
    assert {str(c) for c in pe.conditions} == {"x - 1 != 0", "y - 1 != 0"}


def test_condition_only_variable_survives_cancellation():
    """x/x = y normalizes to y = 1, but x is still an unknown of the problem."""
    pe = parse("x/x = y")
    assert "x" in pe.unknowns
    assert any(c.expression == "x" for c in pe.conditions)


def test_original_equation_validation_rejects_pole():
    pe = parse("1/(x - 1) = 1/(y - 1)")
    assert not pe.is_solution({"x": 1, "y": 1})
    assert pe.is_solution({"x": 2, "y": 2})


def test_numeric_denominator_adds_no_condition():
    pe = parse("x/2 = 1")
    assert not pe.conditions


def test_nested_denominator():
    pe = parse("1/(x*(y + 1)) = 2")
    assert [str(c) for c in pe.conditions] == ["x*(y + 1) != 0"]
    assert pe.conditions[0].evaluate({"x": 1, "y": -1}) is False
    assert pe.conditions[0].evaluate({"x": 1, "y": 1}) is True


def test_negative_power_is_a_denominator():
    pe = parse("x^-2 + y = 1")
    assert [str(c) for c in pe.conditions] == ["x != 0"]


def test_conditions_are_in_source_order():
    pe = parse("4/n = 1/x + 1/y + 1/z", params="n")
    assert [str(c) for c in pe.conditions] == [
        "n != 0", "x != 0", "y != 0", "z != 0"]


def test_conditions_are_three_valued():
    pe = parse("1/(x - 1) = 1/(y - 1)")
    assert pe.conditions_hold({"x": 2, "y": 2}) is True
    assert pe.conditions_hold({"x": 1, "y": 2}) is False
    assert pe.conditions_hold({"x": 2}) is None


def test_condition_json_is_structured():
    import json
    pe = parse("1/(x - 1) = 1/(y - 1)")
    payload = [c.as_dict() for c in pe.conditions]
    assert payload[0] == {"type": "nonzero", "expression": "x - 1",
                          "variables": ["x"], "source": "denominator"}
    json.dumps(payload)


def test_source_equation_is_kept_separate_from_the_model():
    """The cleared polynomial is the working model, not the problem."""
    pe = parse("1/(x - 1) = 1/(y - 1)")
    assert pe.original == "1/(x - 1) = 1/(y - 1)"
    assert pe.poly == pe.poly.parent()("-x + y")
    assert pe.satisfies_original({"x": 1, "y": 1}) is False


def test_validation_is_exact_over_QQ():
    pe = parse("2*x - 1 = 0", domain="QQ")
    from sage.all import QQ
    assert pe.is_solution({"x": QQ(1) / 2})
    assert not pe.is_solution({"x": 0})


def test_validation_is_undetermined_without_parameter_values():
    pe = parse("y^2 = x^3 + k", params="k")
    assert pe.satisfies_original({"x": 3, "y": 5}) is None
    assert pe.is_solution({"x": 3, "y": 5}) is False
    assert pe.is_solution({"x": 3, "y": 5, "k": -2})


def test_parse_still_rejects_nonsense():
    with pytest.raises(ParseError):
        parse("x = x")


# --- conditional identities (brief 3) -------------------------------------

def test_conditional_identity_survives_cancellation():
    """x/x = 1 is true for every x != 0 -- not the unrestricted tautology."""
    pe = parse("x/x = 1")
    assert pe.is_conditional_identity
    assert pe.unknowns == ("x",)
    assert [str(c) for c in pe.conditions] == ["x != 0"]
    assert pe.is_solution({"x": 1})
    assert not pe.is_solution({"x": 0})


def test_unrestricted_identity_is_still_rejected():
    with pytest.raises(ParseError):
        parse("x = x")


def test_conditional_identity_has_no_fabricated_polynomial():
    pe = parse("x/x = 1")
    assert pe.terms == []
    assert pe.poly is None
    assert not pe.is_polynomial


def test_conditional_identity_keeps_every_source_identifier():
    pe = parse("(x - y)/(x - y) = 1")
    assert pe.is_conditional_identity
    assert pe.unknowns == ("x", "y")
    assert [str(c) for c in pe.conditions] == ["x - y != 0"]
    assert pe.is_solution({"x": 3, "y": 1})
    assert not pe.is_solution({"x": 2, "y": 2})


def test_conditional_identity_keeps_domain_and_parameters():
    pe = parse("k*x/(k*x) = 1", params="k", domain="QQ")
    assert pe.is_conditional_identity
    assert pe.domain == "QQ" and pe.params == ("k",)
    assert pe.unknowns == ("x",)


def test_ordinary_equation_is_not_a_conditional_identity():
    assert not parse("x^2 - 61*y^2 = 1").is_conditional_identity
    assert not parse("x/x = y").is_conditional_identity

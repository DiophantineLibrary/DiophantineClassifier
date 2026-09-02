"""Solver layer: concrete instances actually solved via Sage/PARI."""

import pytest

from diophantine_classifier import classify, parse, solve, SolverUnavailable


def assert_valid_solutions(equation, solutions, *, params=(), domain="ZZ"):
    """Every tuple must solve the problem as submitted, exactly.

    Catches failures in denominator handling and in transporting solutions
    back to the user's coordinates.
    """
    pe = parse(equation, params=params, domain=domain)
    for solution in solutions:
        assignment = dict(zip(pe.unknowns, solution))
        assert pe.is_solution(assignment), (equation, solution)


def test_solver_does_not_emit_denominator_pole():
    """Clearing 1/(x-1) = 1/(y-1) gives x = y, whose (1, 1) is a pole."""
    S = solve("1/(x - 1) = 1/(y - 1)")
    assert (1, 1) not in S.first(30)
    assert all(x != 1 and y != 1 for x, y in S.first(30))
    assert_valid_solutions("1/(x - 1) = 1/(y - 1)", S.first(30))


def test_linear():
    s = solve("3*x + 5*y = 1")
    ((x, y),) = s.solutions
    assert 3 * x + 5 * y == 1


def test_linear_empty():
    s = solve("6*x + 9*y = 5")
    assert s.kind == "empty" and s.complete


def test_univariate():
    s = solve("x^2 - 5*x + 6 = 0")
    assert [t[0] for t in s.solutions] == [2, 3]


def test_linear_over_qq_is_not_integer_gcd_problem():
    S = solve("2*x + 4*y = 1", domain="QQ")
    assert S.kind != "empty"
    assert all(2 * x + 4 * y == 1 for x, y in S.first(5))
    assert_valid_solutions("2*x + 4*y = 1", S.first(5), domain="QQ")


def test_linear_over_qq_reaches_genuinely_rational_points():
    """The integer lattice is a proper subset; iteration must leave it."""
    from sage.all import ZZ
    S = solve("2*x + 4*y = 1", domain="QQ")
    sols = S.first(20)
    assert len(set(sols)) == len(sols)
    assert any(x not in ZZ or y not in ZZ for x, y in sols)


def test_linear_nn_never_returns_negative_coordinates():
    try:
        S = solve("3*x + 5*y = 1", domain="NN")
    except SolverUnavailable:
        return
    assert all(x >= 0 and y >= 0 for x, y in S.first(20))


def test_linear_nn_is_complete_where_it_answers():
    S = solve("3*x + 5*y = 47", domain="NN")
    assert S.complete and S.kind == "finite-complete"
    assert S.solutions == [(4, 7), (9, 4), (14, 1)]
    assert_valid_solutions("3*x + 5*y = 47", S.solutions, domain="NN")


def test_linear_nn_infinite_line_stays_nonnegative():
    S = solve("x - y = 1", domain="NN")
    sols = S.first(20)
    assert len(sols) == 20
    assert all(x >= 0 and y >= 0 for x, y in sols)


def test_linear_nn_declines_rather_than_guess():
    """Two-dimensional mixed signs are an unbounded lattice-point count."""
    with pytest.raises(SolverUnavailable):
        solve("x + y - z = 1", domain="NN")


def test_univariate_domains_differ():
    from sage.all import QQ
    assert solve("2*x - 1 = 0").kind == "empty"
    root, = solve("2*x - 1 = 0", domain="QQ").solutions
    assert root == (QQ(1) / 2,) and root[0].parent() is QQ  # exact, not 0.5
    assert solve("x^2 - 4 = 0", domain="NN").solutions == [(2,)]


def test_parametric_univariate_declines_cleanly():
    with pytest.raises(SolverUnavailable):
        solve("x^2 - k = 0", params="k")


def test_high_rank_linear_set_is_iterable():
    eq = " + ".join(f"x{i}" for i in range(10)) + " = 0"
    S = solve(eq)
    assert S.kind == "infinite"
    assert len(S.first(2)) == 2
    assert S.first(2)[0] != S.first(2)[1]
    assert_valid_solutions(eq, S.first(5))


def test_infinite_without_a_stream_is_rejected():
    from diophantine_classifier.solvers import SolutionSet
    with pytest.raises(ValueError):
        SolutionSet(("x",), [(1,)], "infinite", complete=True)


def test_empty_must_be_empty_and_complete():
    from diophantine_classifier.solvers import SolutionSet
    with pytest.raises(ValueError):
        SolutionSet(("x",), [(1,)], "empty", complete=True)
    with pytest.raises(ValueError):
        SolutionSet(("x",), [], "empty", complete=False)


def test_stored_tuples_match_the_variables():
    from diophantine_classifier.solvers import SolutionSet
    with pytest.raises(ValueError):
        SolutionSet(("x", "y"), [(1,)], "finite-complete", complete=True)


# --- exact serialization (brief 7.5) -------------------------------------

def test_solution_set_json_preserves_rationals():
    import json
    d = solve("2*x - 1 = 0", domain="QQ").as_dict()
    assert d["solutions"] == [["1/2"]]
    json.dumps(d)


def test_solution_set_json_keeps_integers_as_integers():
    d = solve("x^2 - 5*x + 6 = 0").as_dict()
    assert d["solutions"] == [[2], [3]]


def test_infinite_solution_set_json_does_not_consume_the_stream():
    d = solve("3*x + 5*y = 1").as_dict()
    assert d["kind"] == "infinite" and d["infinite"] is True
    assert d["solutions"] == [[2, -1]]


# --- dispatch uses the match it is solving (brief 6.1 / 7.1) -------------

def test_solver_receives_its_own_match(monkeypatch):
    """A fallback solver gets the fallback match, not the primary's data."""
    from diophantine_classifier import solvers
    seen = {}

    def capture(cls, match):
        seen["slug"] = match.slug
        seen["data"] = dict(match.data)
        return solvers.SolutionSet(match.transform.normalized_variables, [],
                                   "empty", "", complete=True)

    monkeypatch.setitem(solvers.SOLVERS, "general-polynomial", capture)
    monkeypatch.delitem(solvers.SOLVERS, "linear")
    solve("3*x + 5*y = 1")
    assert seen["slug"] == "general-polynomial"
    assert "degree" in seen["data"]        # its own data ...
    assert "coeffs" not in seen["data"]    # ... not the linear match's


def test_component_of_a_rational_equation_keeps_the_parent_conditions():
    """(x - 2)(x - y)/y = 0 factors; y != 0 still binds each component."""
    cls = classify("(x - 2)*(x - y)/y = 0")
    component, = [c for c in cls.components if c.slug == "linear"]
    S = solve(component)
    sols = S.first(10)
    assert sols
    assert all(y != 0 for _, y in sols)
    assert (0, 0) not in sols
    for x, y in sols:
        assert cls.parsed.is_solution({"x": x, "y": y})


def test_a_wrong_normalization_is_an_error_not_an_empty_stream(monkeypatch):
    """Refuting every distinguished solution means the map is wrong.

    Silently returning a stream that filters everything out would hang the
    first caller who asked for a solution.
    """
    from diophantine_classifier import solvers
    from diophantine_classifier.transforms import negate

    cls = classify("3*x + 5*y = 1")
    match = cls.match_for("linear")
    # a transform that is not an equivalence of this equation
    monkeypatch.setattr(match, "transform", negate(("x", "y"), ("x",)))
    with pytest.raises(AssertionError) as err:
        solve(cls)
    assert "does not invert onto this problem" in str(err.value)


def test_every_wired_solver_returns_valid_solutions():
    """The generic invariant: whatever comes out solves what went in.

    Families whose solver has not landed on this branch are skipped, so the
    check grows with the registry instead of pinning it.
    """
    corpus = [
        ("3*x + 5*y = 1", "ZZ"), ("2*x + 4*y = 1", "QQ"),
        ("x^2 - 5*x + 6 = 0", "ZZ"), ("2*x - 1 = 0", "QQ"),
        ("x^2 - 61*y^2 = 1", "ZZ"), ("5*x^2 - y^2 = 1", "ZZ"),
        ("x^2 - 2*y^2 = 7", "ZZ"), ("x^2 + y^2 = 610", "ZZ"),
        ("x^2 + y^2 + z^2 = 62", "ZZ"), ("3*x^2 + 7*y^2 = 19", "ZZ"),
        ("y^2 = x^3 - 2", "ZZ"), ("u^2 = v^3 - 2", "ZZ"),
        ("x^3 + 2*y^3 = 11", "ZZ"), ("a^2 + b^2 = c^2", "ZZ"),
        # the same shapes with the unknowns in a different order: a tuple
        # is only an answer if its coordinates land on the right variables
        ("7*b^2 = a^2 + 3*c^2", "ZZ"), ("3*z^2 = x^2 + 7*y^2", "ZZ"),
        ("y^2 + x^2 = 2*z^2", "ZZ"), ("610 = y^2 + x^2", "ZZ"),
        ("1/x + 1/y + 1/z = 1", "ZZ"), ("x^2 + 7 = 2^n", "ZZ"),
        # rational inputs whose cleared model has a pole at its stored
        # witness, a conditional identity, and a reduced repeated factor
        ("1/(x - 1) = 1/(y - 1)", "ZZ"), ("1/x = 1/y", "ZZ"),
        ("x/x = y", "ZZ"), ("x/x = 1", "ZZ"), ("x/x = 1", "QQ"),
        ("(x + y)^2 = 0", "ZZ"), ("(x + y)^2 = 0", "QQ"),
    ]
    for equation, domain in corpus:
        try:
            S = solve(equation, domain=domain)
        except SolverUnavailable:
            continue
        assert_valid_solutions(equation, S.first(8), domain=domain)


# --- expected poles are filtered, not asserted (brief 7.1) ---------------

def test_zero_particular_solution_can_be_a_pole():
    """The homogeneous linear model's particular solution is (0, 0), which
    is a pole of the original -- that is a filtered point, not a wrong map."""
    S = solve("1/x = 1/y")
    sols = S.first(20)
    assert sols
    assert all(x == y and x != 0 for x, y in sols)
    assert_valid_solutions("1/x = 1/y", sols)


def test_condition_only_free_variable_is_filtered_not_asserted():
    S = solve("x/x = y")
    sols = S.first(20)
    assert sols
    assert all(x != 0 and y == 1 for x, y in sols)
    assert_valid_solutions("x/x = y", sols)


def test_an_infinite_set_keeps_its_stream_when_the_witness_is_a_pole():
    S = solve("1/x = 1/y")
    assert S.kind == "infinite"
    assert S.complete
    assert len(S.first(5)) == 5


def test_bad_coordinate_map_still_raises(monkeypatch):
    """The synthetic wrong-normalization check must survive the reordering."""
    from diophantine_classifier.transforms import negate
    cls = classify("3*x + 5*y = 1")
    match = cls.match_for("linear")
    monkeypatch.setattr(match, "transform", negate(("x", "y"), ("x",)))
    with pytest.raises(AssertionError) as err:
        solve(cls)
    assert "does not invert onto this problem" in str(err.value)


# --- the effective transform is what solvers travel back through (7.2) ---

def test_finalizer_uses_the_reduction_as_well_as_the_match_transform():
    """(x + y)^2 = 0 reduces, then the linear match maps coordinates; a
    solution has to come back through both."""
    S = solve("(x + y)^2 = 0")
    sols = S.first(6)
    assert sols
    assert all(x + y == 0 for x, y in sols)
    assert_valid_solutions("(x + y)^2 = 0", sols)


def test_reduction_and_rename_compose_for_the_solver(monkeypatch):
    from diophantine_classifier.transforms import rename
    cls = classify("(x + y)^2 = 0")
    match = cls.match_for("linear")
    # a nontrivial matcher rename on top of the nontrivial reduction
    monkeypatch.setattr(match, "transform",
                        rename([("a", "x"), ("b", "y")]))
    S = solve(cls)
    assert S.variables == ("x", "y")
    assert all(x + y == 0 for x, y in S.first(5))


# --- transform conditions are enforced (brief 7.3) -----------------------

def test_transform_condition_filters_solver_output(monkeypatch):
    """A map defined only away from a divisor must not hand back points on
    it, even when they solve the equation."""
    from diophantine_classifier.conditions import NonzeroCondition
    from diophantine_classifier.transforms import identity
    cls = classify("3*x + 5*y = 1")
    match = cls.match_for("linear")
    monkeypatch.setattr(match, "transform",
                        identity(("x", "y"),
                                 conditions=(NonzeroCondition("x - 2",
                                                              ("x",)),)))
    sols = solve(cls).first(6)
    assert sols
    assert all(x != 2 for x, y in sols)          # (2, -1) is on the divisor
    assert all(3 * x + 5 * y == 1 for x, y in sols)


# --- component solvers are safe (brief 7.4) ------------------------------

def test_component_solver_never_returns_the_lower_dimensional_point_only():
    cls = classify("(x - 2)*(x - y)/y = 0")
    comp = next(c for c in cls.components if str(c.working.poly) == "x - 2")
    with pytest.raises(SolverUnavailable) as err:
        solve(comp)
    assert "free" in str(err.value) and "y" in str(err.value)


def test_component_without_free_variables_still_solves():
    cls = classify("(x^2 - 2)*(x - y) = 0")
    comp = next(c for c in cls.components if not c.free_variables)
    S = solve(comp)
    assert_valid_solutions("(x^2 - 2)*(x - y) = 0", S.first(5))


# --- conditional identities (brief 3.4) ----------------------------------

def test_conditional_identity_is_solved_by_its_conditions():
    S = solve("x/x = 1")
    assert S.kind == "infinite" and S.complete
    sols = S.first(10)
    assert len(sols) == 10
    assert all(x != 0 for x, in sols)
    assert_valid_solutions("x/x = 1", sols)


def test_conditional_identity_respects_the_domain():
    assert all(x > 0 for x, in solve("x/x = 1", domain="NN").first(5))
    rationals = solve("x/x = 1", domain="QQ").first(8)
    assert len(set(rationals)) == 8
    assert_valid_solutions("x/x = 1", rationals, domain="QQ")


def test_conditional_identity_in_two_variables():
    S = solve("(x - y)/(x - y) = 1")
    sols = S.first(10)
    assert all(x != y for x, y in sols)
    assert_valid_solutions("(x - y)/(x - y) = 1", sols)



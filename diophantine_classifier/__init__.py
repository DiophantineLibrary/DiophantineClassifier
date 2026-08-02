"""DiophantineClassifier: classify Diophantine equations into named families.

Part of the Diophantine Library (https://github.com/DiophantineLibrary).
Runs on top of SageMath; some solver suggestions target Magma or PARI.

Basic usage (inside Sage)::

    from diophantine_classifier import classify, solve

    classify("x^2 - 61*y^2 = 1")           # Pell equation, D = 61
    classify("y^2 = x^3 + k", params="k")  # Mordell equation
    classify("x^2 + 7 = 2^n")              # Ramanujan-Nagell
    solve("x^3 + 2*y^3 = 11")              # [(3, -2)] via PARI's thue
"""

__version__ = "0.1.0"

from .parsing import parse, ParseError, UnsupportedEquationError
from .registry import families, family, Family, ancestors, validate
from .classify import classify, Classification
from .solvers import solve, SolutionSet, SolverUnavailable

__all__ = [
    "classify", "Classification", "solve", "SolutionSet", "SolverUnavailable",
    "parse", "ParseError", "UnsupportedEquationError",
    "families", "family", "Family", "ancestors", "validate", "__version__",
]

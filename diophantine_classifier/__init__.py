"""
DiophantineClassifier: classify Diophantine equations into named families.

Part of the Diophantine Library (https://github.com/DiophantineLibrary).
Runs on top of SageMath; some solver suggestions target Magma or PARI.

Basic usage (inside Sage)::

    from diophantine_classifier import classify

    classify("3*x + 5*y = 1")              # linear Diophantine equation
    classify("x^2 - 5*x + 6 = 0")          # univariate
"""

__version__ = "0.1.0"

from .parsing import parse, ParseError, UnsupportedEquationError
from .registry import (families, family, Family, ancestors, lineage_paths,
                       lineage_graph, validate)
from .classify import classify, Classification

__all__ = [
    "classify", "Classification",
    "parse", "ParseError", "UnsupportedEquationError",
    "families", "family", "Family", "ancestors", "lineage_paths",
    "lineage_graph", "validate", "__version__",
]

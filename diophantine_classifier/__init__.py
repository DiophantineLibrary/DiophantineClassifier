"""
DiophantineClassifier: classify Diophantine equations into named families.

Part of the Diophantine Library (https://github.com/DiophantineLibrary).
Runs on top of SageMath; some solver suggestions target Magma or PARI.
"""

__version__ = "0.1.0"

from .parsing import parse, ParseError, UnsupportedEquationError
from .registry import (families, family, Family, ancestors, lineage_paths,
                       lineage_graph, validate)

__all__ = [
    "parse", "ParseError", "UnsupportedEquationError",
    "families", "family", "Family", "ancestors", "lineage_paths",
    "lineage_graph", "validate", "__version__",
]

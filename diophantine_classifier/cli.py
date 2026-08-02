"""Command-line interface: ``dioclassify 'x^2 - 61*y^2 = 1' [--solve]``.

Must run under Sage's Python (``sage -python -m diophantine_classifier.cli ...``
or the ``dioclassify`` entry point after ``sage -pip install -e .``).
"""

import argparse
import json
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="dioclassify",
        description="Classify a Diophantine equation into a named family "
                    "(Diophantine Library).",
    )
    parser.add_argument("equation", help="e.g. 'x^2 - 61*y^2 = 1'")
    parser.add_argument("--params", default="",
                        help="comma-separated parameter names, e.g. 'k' or 'd,n'")
    parser.add_argument("--domain", default="ZZ", choices=["ZZ", "NN", "QQ"])
    parser.add_argument("--solve", action="store_true",
                        help="also attempt to solve via Sage/PARI")
    parser.add_argument("--json", action="store_true",
                        help="machine-readable output")
    args = parser.parse_args(argv)

    from .classify import classify
    from .solvers import solve, SolverUnavailable

    cls = classify(args.equation, params=args.params, domain=args.domain)
    if args.json:
        out = cls.as_dict()
        if args.solve:
            try:
                sols = solve(cls)
                out["solutions"] = {
                    "kind": sols.kind,
                    "complete": sols.complete,
                    "variables": list(sols.variables),
                    "solutions": [[int(x) for x in s] for s in sols.solutions],
                    "description": sols.description,
                }
            except SolverUnavailable as err:
                out["solutions"] = {"error": str(err)}
        print(json.dumps(out, indent=2))
        return 0
    print(cls.explain())
    if args.solve:
        try:
            print()
            print(solve(cls))
        except SolverUnavailable as err:
            print(f"\nsolver unavailable: {err}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

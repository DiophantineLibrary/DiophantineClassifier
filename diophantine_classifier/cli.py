r"""Command-line interface: ``dioclassify 'x^2 - 61*y^2 = 1' [--solve]``.

Must run under Sage's Python (``sage -python -m diophantine_classifier.cli
...`` or the ``dioclassify`` entry point after ``sage -pip install -e .``).

EXAMPLES::

    sage: from diophantine_classifier.cli import main
    sage: main(["3*x + 5*y = 1"])                    # doctest: +ELLIPSIS
    3*x + 5*y = 1
      family: linear — Linear Diophantine equation   [P1, solved]
    ...
    0
"""

import argparse
import json
import sys


def build_parser():
    r"""Build the argument parser (separated out for testing).

    OUTPUT: an :class:`argparse.ArgumentParser`

    EXAMPLES::

        sage: from diophantine_classifier.cli import build_parser
        sage: args = build_parser().parse_args(["x^2 = 2", "--solve"])
        sage: args.solve, args.domain
        (True, 'ZZ')
    """
    parser = argparse.ArgumentParser(
        prog="dioclassify",
        description="Classify a Diophantine equation into a named family "
                    "(Diophantine Library).",
    )
    parser.add_argument("equation", help="e.g. 'x^2 - 61*y^2 = 1'")
    parser.add_argument("--params", default="",
                        help="comma-separated parameter names, e.g. 'k' or "
                             "'d,n'")
    parser.add_argument("--domain", default="ZZ", choices=["ZZ", "NN", "QQ"])
    parser.add_argument("--solve", action="store_true",
                        help="also attempt to solve via Sage/PARI")
    parser.add_argument("--json", action="store_true",
                        help="machine-readable output")
    return parser


def main(argv=None):
    r"""Entry point: classify (and optionally solve) one equation.

    INPUT:

    - ``argv`` -- (default: ``sys.argv[1:]``) list of command-line arguments

    OUTPUT: process exit code (0)

    EXAMPLES::

        sage: from diophantine_classifier.cli import main
        sage: rc = main(["x^2 - 5*x + 6 = 0", "--solve"])
        x^2 - 5*x + 6 = 0
        ...
        SolutionSet(finite-complete, vars=['x'])
          solutions: (2,), (3,)
          all roots in the domain
        sage: rc
        0

    JSON mode emits the ``as_dict`` payload::

        sage: import json, io, contextlib
        sage: buf = io.StringIO()
        sage: with contextlib.redirect_stdout(buf):
        ....:     rc = main(["x^2 - 61*y^2 = 1", "--json"])
        sage: json.loads(buf.getvalue())["family"]
        'pell'
    """
    args = build_parser().parse_args(argv)

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
                    "solutions": [[int(x) for x in s]
                                  for s in sols.solutions],
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

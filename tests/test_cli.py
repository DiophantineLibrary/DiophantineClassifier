"""The command-line interface, and the JSON payload the website consumes."""

import json

from diophantine_classifier.cli import main


def payload(argv, capsys):
    assert main(argv) == 0
    return json.loads(capsys.readouterr().out)


def test_cli_json_preserves_rational_root(capsys):
    """int(x) used to turn the only root, 1/2, into 0."""
    out = payload(["2*x - 1 = 0", "--domain", "QQ", "--solve", "--json"],
                  capsys)
    assert out["solutions"]["solutions"] == [["1/2"]]
    assert out["solutions"]["complete"] is True


def test_cli_json_keeps_integers_as_integers(capsys):
    out = payload(["x^2 - 5*x + 6 = 0", "--solve", "--json"], capsys)
    assert out["solutions"]["solutions"] == [[2], [3]]


def test_cli_json_is_serializable_and_structured(capsys):
    out = payload(["1/x = y", "--json"], capsys)
    assert out["equation"] == "1/x = y"
    assert out["conditions"][0]["expression"] == "x"
    assert out["transform"]["identity"] is True
    assert isinstance(out["lineage_paths"], list)


def test_cli_json_reports_an_unavailable_solver(capsys):
    out = payload(["x + y - z = 1", "--domain", "NN", "--solve", "--json"],
                  capsys)
    assert "error" in out["solutions"]


def test_cli_text_mode_runs(capsys):
    assert main(["3*x + 5*y = 1", "--solve"]) == 0
    text = capsys.readouterr().out
    assert "family: linear" in text
    assert "SolutionSet" in text


def test_cli_domain_is_threaded_through(capsys):
    out = payload(["x^2 - 5*x + 6 = 0", "--domain", "NN", "--json"], capsys)
    assert out["domain"] == "NN"

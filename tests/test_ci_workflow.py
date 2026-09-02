"""The CI coverage gate must fail when coverage does.

The step is shell, so the only honest test is to run it: the workflow's own
script is extracted and executed against a stub `sage` whose output we
control.  Piping the coverage command into `tee` (or letting a missing SCORE
line pass vacuously) is exactly the regression this catches.
"""

import os
import stat
import subprocess

import pytest
import yaml

WORKFLOW = os.path.join(os.path.dirname(__file__), "..", ".github",
                        "workflows", "ci.yml")

FULL_SCORE = "\n".join(
    f"SCORE diophantine_classifier/{name}.py: 100.0% (9 of 9)"
    for name in ["classify", "parsing", "solvers"])


def coverage_script():
    """The `run:` body of the workflow's docstring-coverage step."""
    with open(WORKFLOW) as fobj:
        workflow = yaml.safe_load(fobj)
    for step in workflow["jobs"]["test"]["steps"]:
        if "coverage" in step.get("name", "").lower():
            return step["run"]
    raise AssertionError("no docstring-coverage step in the workflow")


def run_gate(tmp_path, output, exit_code=0):
    """Run the gate with a stub `sage` printing `output` and exiting."""
    stub = tmp_path / "sage"
    stub.write_text("#!/bin/sh\ncat <<'EOF'\n%s\nEOF\nexit %d\n"
                    % (output, exit_code))
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    env = dict(os.environ, PATH=f"{tmp_path}{os.pathsep}{os.environ['PATH']}")
    return subprocess.run(["sh", "-c", coverage_script()], cwd=tmp_path,
                          env=env, capture_output=True, text=True).returncode


def test_the_gate_is_not_piped_into_tee():
    """A pipeline hands the step tee's exit status, not sage's."""
    script = coverage_script()
    assert "sage --coverage" in script
    command = next(line for line in script.splitlines()
                   if "sage --coverage" in line)
    assert "|" not in command, command


def test_a_full_score_passes(tmp_path):
    assert run_gate(tmp_path, FULL_SCORE) == 0


def test_a_failing_coverage_command_fails_the_job(tmp_path):
    assert run_gate(tmp_path, FULL_SCORE, exit_code=2) != 0


def test_missing_score_output_fails_the_job(tmp_path):
    assert run_gate(tmp_path, "sage: something went wrong") != 0


def test_a_score_below_100_fails_the_job(tmp_path):
    degraded = FULL_SCORE + (
        "\nSCORE diophantine_classifier/matchers.py: 87.5% (7 of 8)")
    assert run_gate(tmp_path, degraded) != 0


def test_one_perfect_file_does_not_excuse_the_rest(tmp_path):
    """`case "$score" in *': 100.0%'*` would pass this; it must not."""
    mixed = ("SCORE diophantine_classifier/a.py: 100.0% (9 of 9)\n"
             "SCORE diophantine_classifier/b.py: 50.0% (1 of 2)")
    assert run_gate(tmp_path, mixed) != 0

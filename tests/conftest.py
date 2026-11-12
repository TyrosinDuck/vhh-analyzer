"""
What the suite needs on the machine, checked once, before anything runs.

ANARCI shells out to HMMER, which is a binary and not a Python package:
Community Cloud installs it from `packages.txt`, and nothing in
`requirements-dev.txt` can bring it in. Without it the app still works, the
regex annotator takes over, but the smoke test drives the Comparison view,
which calls ANARCI, and failed with

    Could not compare: [Errno 2] No such file or directory: 'hmmscan'

buried in an assertion about a Streamlit element. The error is accurate and
tells a contributor nothing about what to install, which is the point of this
file: the missing binary is reported once, at the top of the run, with the
command that fixes it.

This warns rather than skips. A suite that quietly skipped the annotator tests
on a machine without HMMER would pass on the one machine where the app's
default annotator is never exercised.

`HAS_HMMER` is here rather than recomputed in each test file so that the
warning above and the tests that branch on it cannot disagree. A test that
branches on it must assert something in *both* branches: without HMMER the
reference comparison is expected to say what it needs, and a test that merely
skipped there would hide the day that message stops being drawn.
"""

from __future__ import annotations

import shutil

HAS_HMMER = shutil.which("hmmscan") is not None


def pytest_configure(config):
    if HAS_HMMER:
        return
    message = (
        "HMMER is not installed (no `hmmscan` on PATH), so ANARCI cannot run "
        "and the tests that annotate through it will fail. Install it with "
        "`apt-get install hmmer`, `brew install hmmer` or "
        "`conda install -c bioconda hmmer`. The deployed app gets it from "
        "packages.txt.")
    reporter = config.pluginmanager.get_plugin("terminalreporter")
    if reporter is not None:
        reporter.write_line("")
        reporter.write_line(message, yellow=True, bold=True)
    else:
        print(message)

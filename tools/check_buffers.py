"""Check data/buffers/*.csv before committing them.

    python3 tools/check_buffers.py            # from the repository root
    python3 tools/check_buffers.py path/to/buffers

Runs the loader the app uses, `buffers.load`, so a file this accepts is a
file the app shows in full, and the two cannot disagree. Prints every problem
with its file and line, and exits non-zero if there is any. It checks
structure, not chemistry: see data/buffers/README.md.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import buffers  # noqa: E402


def main(argv: list[str]) -> int:
    folder = Path(argv[1]) if len(argv) > 1 else buffers.DATA_DIR
    book, problems = buffers.load(folder)
    for problem in problems:
        print(f"ERROR    {problem}")
    if problems:
        print(f"\n{len(problems)} problem(s).")
        return 1
    print(f"{len(book['recipes'])} recipe(s), {len(book['compounds'])} "
          "compound(s), no problem in the structure. Formulas and amounts "
          "still need checking against your bottles and protocols.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

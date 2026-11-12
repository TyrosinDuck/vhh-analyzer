"""Regenerate the section table in README.md from `shell.SECTIONS`.

The table used to be typed into the README, which is one more list of the
sections to keep in step with the one in `shell.py`. It drifted, as a copied
list does: by the time this script was written the README's own file table was
wrong about eleven files out of thirteen, and its section descriptions were a
paraphrase of `SECTIONS` from a few releases earlier.

So the README now holds generated text between two markers, and
`tests/test_metadata.py` fails if what is between them is not what this script
would write. Regenerating is one command and the test says when it is needed:

    python tools/readme_sections.py

Nothing here is imported by the app.

The description in the table is `SECTIONS`' own, which is also what the home
page card and the sidebar show, so a reader of the README and a user of the
app are told the same thing in the same words.
"""

from __future__ import annotations

import sys
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

import shell  # noqa: E402

README = ROOT / "README.md"
START = "<!-- generated: sections, by tools/readme_sections.py -->"
END = "<!-- /generated: sections -->"


def table() -> str:
    """The Markdown table, exactly as it should appear between the markers."""
    lines = ["| Group | Section | What it does |", "|---|---|---|"]
    for group, keys in shell.sections_by_group().items():
        for n, key in enumerate(keys):
            icon, label, description, _group = shell.SECTIONS[key]
            # The group cell is written once per group, as in a printed table.
            cell = f"**{group}**" if n == 0 else ""
            flat = " ".join(description.split())
            lines.append(f"| {cell} | **{icon} {label}** | {flat} |")
    return "\n".join(lines)


def block() -> str:
    """The table with its markers, which is what the README must contain."""
    return f"{START}\n\n{table()}\n\n{END}"


def current(text: str) -> str:
    """What sits between the markers in `text`, markers included."""
    start = text.index(START)
    end = text.index(END) + len(END)
    return text[start:end]


def main() -> int:
    text = README.read_text(encoding="utf-8")
    wanted = block()
    if current(text) == wanted:
        print("README.md section table is already up to date")
        return 0
    README.write_text(text.replace(current(text), wanted), encoding="utf-8")
    print(textwrap.dedent(f"""\
        README.md section table rewritten from shell.SECTIONS
        ({len(shell.SECTIONS)} sections in {len(shell.sections_by_group())} groups)"""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

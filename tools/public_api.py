"""Write the snapshot of vhh_core's public API that tests/test_public_api.py
holds every later 1.x release to.

    python3 tools/public_api.py

Run it only when the promise itself moves: after a MAJOR release, whose whole
point is that the old signatures no longer hold, or to record functions added
in a MINOR release so that they are protected from then on. Running it to make
a failing test pass on a MINOR or PATCH release defeats the test: the failure
was the test saying that the change breaks somebody's notebook.

The snapshot is deliberately a copy, which this project otherwise avoids. A
contract test has to compare the code against something the code cannot
change by itself; a list derived from the module at test time would always
agree with it. What keeps the copy from going stale is that it is written by
this script from the module, never edited by hand, and that the test refuses a
snapshot taken under a different MAJOR number.
"""

from __future__ import annotations

import inspect
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SNAPSHOT = ROOT / "tests" / "compat" / "vhh_core_api.json"


def describe(module) -> dict:
    """Each name in ``module.__all__``, with the parameters of the callables:
    name, kind and whether it has a default. Defaults' values are not kept,
    since changing one is not a break of the signature."""
    names = {}
    for name in sorted(module.__all__):
        obj = getattr(module, name)
        if not callable(obj) or name == "__version__":
            names[name] = {"callable": False}
            continue
        params = []
        for p in inspect.signature(obj).parameters.values():
            params.append({"name": p.name, "kind": p.kind.name,
                           "default": p.default is not inspect.Parameter.empty})
        names[name] = {"callable": True, "parameters": params}
    return names


def main() -> int:
    sys.path.insert(0, str(ROOT))
    import shell
    import vhh_core

    snapshot = {
        "release_major": int(shell.RELEASE.split(".")[0]),
        "taken_at_release": shell.RELEASE,
        "names": describe(vhh_core),
    }
    SNAPSHOT.parent.mkdir(parents=True, exist_ok=True)
    SNAPSHOT.write_text(json.dumps(snapshot, indent=1) + "\n", encoding="utf-8")
    print(f"wrote {SNAPSHOT.relative_to(ROOT)}: "
          f"{len(snapshot['names'])} names at release {shell.RELEASE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

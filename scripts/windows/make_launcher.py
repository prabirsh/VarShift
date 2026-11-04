"""Make varshift.exe, the same kind of small launcher pip makes for its own commands.

A .cmd file runs through Command Prompt, which reads > as "write to a file" unless it is in quotes, and
PowerShell removes the quotes from any argument without a space. So NM_004333.6:c.1799T>A would lose its
ending and leave a stray file called A. An .exe launcher never goes through Command Prompt.

install.ps1 falls back on a .cmd file if this can't run. Usage: make_launcher.py BIN_FOLDER SOURCE_FOLDER
"""

import sys

_SCRIPT = """import sys
sys.path.insert(0, {source!r})
from %(module)s import %(import_name)s
if __name__ == "__main__":
    sys.exit(%(func)s())
"""


def main() -> int:
    from pip._vendor.distlib.scripts import ScriptMaker

    bin_folder, source_folder = sys.argv[1:3]
    maker = ScriptMaker(None, bin_folder, add_launchers=True)
    maker.clobber = True
    maker.variants = {""}
    maker.executable = sys.executable
    maker.script_template = _SCRIPT.format(source=source_folder)
    made = maker.make("varshift = varshift.cli:main")
    print(made[0] if made else "")
    return 0 if made else 1


if __name__ == "__main__":
    sys.exit(main())

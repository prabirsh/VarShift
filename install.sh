#!/bin/sh
# Install the VarSight command. It needs only Python 3.10 or newer: no pip, no sudo, no packages.
#
#   ./install.sh                                    command only, then run the demo
#   ./install.sh --setup --data-dir ~/varsight-data  also download and build all the data
#
set -eu

usage() {
    cat <<'EOF'
Usage: ./install.sh [options]

  --bin-dir DIR     where to put the varsight command (default: ~/.local/bin)
  --data-dir DIR    data folder (default: ~/varsight-data)
  --build NAME      hg19, hg38 or both (default: both)
  --setup           also download and build the data now (about 8 GB transfer, see README)
  --no-demo         skip the bundled demo
  --test            run the unit tests as well
  -h, --help        show this help
EOF
}

bin_dir="${HOME}/.local/bin"
data_dir=""
build="both"
run_setup=0
run_demo=1
run_tests=0

while [ $# -gt 0 ]; do
    case "$1" in
        --bin-dir) bin_dir="${2:?--bin-dir needs a folder}"; shift 2 ;;
        --data-dir) data_dir="${2:?--data-dir needs a folder}"; shift 2 ;;
        --build) build="${2:?--build needs hg19, hg38 or both}"; shift 2 ;;
        --setup) run_setup=1; shift ;;
        --no-demo) run_demo=0; shift ;;
        --test) run_tests=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "unknown option: $1" >&2; usage >&2; exit 2 ;;
    esac
done

case "$build" in
    hg19|hg38|both) ;;
    *) echo "--build must be hg19, hg38 or both" >&2; exit 2 ;;
esac

repo_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd -P)

python_bin=""
for candidate in python3 python; do
    if command -v "$candidate" >/dev/null 2>&1 \
        && "$candidate" -c 'import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)' 2>/dev/null; then
        python_bin=$(command -v "$candidate")
        break
    fi
done
if [ -z "$python_bin" ]; then
    echo "VarSight needs Python 3.10 or newer, and none was found." >&2
    echo "On Ubuntu or WSL: sudo apt update && sudo apt install python3" >&2
    exit 1
fi
if ! "$python_bin" -c 'import sqlite3, gzip, mmap, ssl' 2>/dev/null; then
    echo "$python_bin is missing a standard module (sqlite3, gzip, mmap or ssl)." >&2
    exit 1
fi
echo "Python: $("$python_bin" --version) at $python_bin"

mkdir -p "$bin_dir"
launcher="$bin_dir/varsight"
# Quote the paths for the shell, so a folder name with spaces or quotes stays safe.
quote() { printf "'%s'" "$(printf '%s' "$1" | sed "s/'/'\\\\''/g")"; }
{
    echo '#!/bin/sh'
    echo "PYTHONPATH=$(quote "$repo_dir/src") exec $(quote "$python_bin") -m varsight \"\$@\""
} > "$launcher.tmp"
chmod 755 "$launcher.tmp"
mv "$launcher.tmp" "$launcher"
echo "Command: $launcher"

"$launcher" --version

if [ "$run_tests" -eq 1 ]; then
    echo
    echo "Running the unit tests"
    (cd "$repo_dir" && PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 "$python_bin" -m unittest discover -s tests -t .)
fi

if [ "$run_demo" -eq 1 ]; then
    echo
    echo "Running the demo (no download)"
    "$launcher" demo --demo-dir "$repo_dir/demo" --out "${PWD}/demo_output" --force
fi

if [ "$run_setup" -eq 1 ]; then
    echo
    target="${data_dir:-$HOME/varsight-data}"
    "$launcher" setup --data-dir "$target" --build "$build" --yes
fi

case ":${PATH}:" in
    *":$bin_dir:"*) ;;
    *)
        echo
        echo "Add the command to your PATH by running this once:"
        echo "  echo 'export PATH=\"$bin_dir:\$PATH\"' >> ~/.bashrc && . ~/.bashrc"
        ;;
esac

echo
echo "Installed. Next steps:"
if [ "$run_setup" -eq 0 ]; then
    echo "  varsight setup --data-dir ~/varsight-data    download and build the data (about 8 GB transfer)"
fi
echo "  varsight annotate \"BRAF V600E\" --build hg19    annotate a variant"
echo "  varsight --help                                 see every command"

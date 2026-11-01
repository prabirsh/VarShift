#!/bin/sh
# Install the VarSight command on Linux or WSL. It needs only Python 3.10 or newer: no pip, no sudo.
#
#   sh install.sh               install the command, run the demo, then offer to set up the databases
#   sh install.sh setup=no      install the command and the demo only
#
set -eu

usage() {
    cat <<'EOF'
Usage: sh install.sh [option=value ...]

  bin=DIR      where to put the varsight command (default: ~/.local/bin)
  demo=no      skip the bundled demo
  setup=no     don't offer to set up the databases
  test=yes     run the unit tests as well
EOF
}

bin_dir="${HOME}/.local/bin"
run_demo=yes
run_setup=ask
run_tests=no

for argument in "$@"; do
    case "$argument" in
        bin=*) bin_dir="${argument#bin=}" ;;
        demo=*) run_demo="${argument#demo=}" ;;
        setup=*) run_setup="${argument#setup=}" ;;
        test=*) run_tests="${argument#test=}" ;;
        -h|help) usage; exit 0 ;;
        *) echo "unknown option: $argument" >&2; usage >&2; exit 2 ;;
    esac
done

repo_dir=$(CDPATH='' cd "$(dirname "$0")" && pwd -P)

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
echo "Python: $("$python_bin" -V) at $python_bin"

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
"$launcher" version

if [ "$run_tests" = "yes" ]; then
    echo
    echo "Running the unit tests"
    (cd "$repo_dir" && PYTHONPATH=src PYTHONDONTWRITEBYTECODE=1 "$python_bin" -m unittest discover -s tests -t .)
fi

if [ "$run_demo" = "yes" ]; then
    echo
    echo "Running the demo (no download)"
    "$launcher" demo demodir="$repo_dir/demo" out="${PWD}/demo_output" replace
fi

case ":${PATH}:" in
    *":$bin_dir:"*) ;;
    *)
        echo
        echo "To use the varsight command from anywhere, run this once:"
        echo "  echo 'export PATH=\"$bin_dir:\$PATH\"' >> ~/.bashrc && . ~/.bashrc"
        ;;
esac

if [ "$run_setup" = "ask" ] && [ -t 0 ]; then
    echo
    printf "Set up the databases now? It asks whether you want hg19, hg38 or both. yes or no [no] "
    read -r answer || answer=no
    case "$answer" in
        y|Y|yes|Yes|YES) "$launcher" setup ;;
    esac
fi

echo
echo "Installed. Next steps:"
echo "  varsight setup                               choose hg19, hg38 or both, and download"
echo "  varsight annotate \"BRAF V600E\" build=hg19    annotate a variant"
echo "  varsight help                                see every command"

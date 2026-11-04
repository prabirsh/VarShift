#!/bin/sh
# One data update, for cron or Windows Task Scheduler with WSL. Output goes to a dated log in the data
# folder. Updates that can change names (transcripts) are reported but wait for your approval.
#
# Usage: weekly-update.sh [DATA_FOLDER] [VARSHIFT_COMMAND]
set -u

data_dir="${1:-${VARSHIFT_DATA:-$HOME/varshift-data}}"
command="${2:-$HOME/.local/bin/varshift}"

log_dir="$data_dir/reports"
mkdir -p "$log_dir"
log_file="$log_dir/weekly-$(date -u +%Y%m%dT%H%M%SZ).log"

status=0
"$command" update data="$data_dir" >"$log_file" 2>&1 || status=$?
if [ "$status" -ne 0 ]; then
    echo "VarShift weekly update reported a problem; see $log_file" >&2
fi
exit "$status"

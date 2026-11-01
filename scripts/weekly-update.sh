#!/bin/sh
# One weekly data update, for cron or Windows Task Scheduler. Output goes to a dated log in the data
# folder. Results-changing sources (transcripts) are reported but wait for your approval.
#
# Usage: weekly-update.sh [DATA_DIR] [VARSIGHT_COMMAND]
set -u

data_dir="${1:-${VARSIGHT_DATA:-$HOME/varsight-data}}"
command="${2:-$HOME/.local/bin/varsight}"

log_dir="$data_dir/reports"
mkdir -p "$log_dir"
log_file="$log_dir/weekly-$(date -u +%Y%m%dT%H%M%SZ).log"

status=0
"$command" data update --data-dir "$data_dir" >"$log_file" 2>&1 || status=$?
if [ "$status" -ne 0 ]; then
    echo "VarSight weekly update reported a problem; see $log_file" >&2
fi
exit "$status"

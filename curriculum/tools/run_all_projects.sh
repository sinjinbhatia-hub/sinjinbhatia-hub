#!/usr/bin/env bash
# Run every project.py headless and report failures. Usage: bash curriculum/tools/run_all_projects.sh
set -u
root="$(cd "$(dirname "$0")/.." && pwd)"
fail=0
while IFS= read -r f; do
  d="$(dirname "$f")"
  if (cd "$d" && MPLBACKEND=Agg timeout 300 python3 project.py > /tmp/proj_out.txt 2>&1); then
    echo "PASS  ${d#$root/}"
  else
    echo "FAIL  ${d#$root/}"; tail -n 15 /tmp/proj_out.txt | sed 's/^/      /'; fail=1
  fi
done < <(find "$root" -name project.py | sort)
exit $fail

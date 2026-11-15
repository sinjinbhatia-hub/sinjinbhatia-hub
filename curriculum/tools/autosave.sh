#!/usr/bin/env bash
# Commit and push curriculum work every 3 minutes while a build round runs, so a usage-limit
# kill or a reclaimed container loses at most a few minutes. Start from the repo root with:
#   setsid nohup bash curriculum/tools/autosave.sh > /tmp/autosave.log 2>&1 < /dev/null & disown
root="$(cd "$(dirname "$0")/../.." && pwd)"
cd "$root" || exit 1
branch="$(git branch --show-current)"
while true; do
  sleep 180
  if [ -n "$(git status --porcelain curriculum)" ]; then
    git add -A curriculum
    git commit -q -m "wip: autosave of module drafts in progress

Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_013cPnMkRA9yM69ycnWNbB9p" && git push -q -u origin "$branch" 2>/dev/null
    echo "$(date -u +%H:%M) autosaved"
  fi
done

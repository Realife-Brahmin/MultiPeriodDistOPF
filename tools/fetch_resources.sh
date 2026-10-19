#!/usr/bin/env bash
# Walks the repo for resources/MANIFEST.txt files and fetches any listed file
# that isn't already on disk. Run from anywhere; paths resolve relative to the
# repo root. Safe to re-run: present files are left untouched.
set -uo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

present=0
fetched=0
failed=0

while IFS= read -r manifest; do
  dir="$(dirname "$manifest")"
  while IFS='|' read -r filename url description; do
    filename="$(echo "$filename" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    [ -z "$filename" ] && continue
    case "$filename" in \#*) continue ;; esac
    url="$(echo "${url:-}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    description="$(echo "${description:-}" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')"
    target="$dir/$filename"

    if [ -f "$target" ]; then
      echo "[present] $target"
      present=$((present + 1))
    elif [ -z "$url" ]; then
      echo "[MISSING, no url given] $target"
      failed=$((failed + 1))
    else
      echo "[fetching] $target <- $url"
      if curl -fL -o "$target" "$url"; then
        echo "[fetched]  $target ($description)"
        fetched=$((fetched + 1))
      else
        echo "[FAILED]   $target <- $url"
        rm -f "$target"
        failed=$((failed + 1))
      fi
    fi
  done < "$manifest"
done < <(find . -path '*/external/*' -prune -o -path '*/resources/MANIFEST.txt' -print)

echo ""
echo "Summary: $present present, $fetched fetched, $failed missing/failed."
[ "$failed" -eq 0 ]

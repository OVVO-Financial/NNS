#!/usr/bin/env bash
# Stage the examples/ directory as a browsable GitHub Pages section.
#
# Copies the published artifact formats (html, pdf, images) verbatim, renders
# every Markdown study to styled standalone HTML with pandoc, and renders
# examples/index.md as the landing page. Relative *.md links are rewritten to
# *.html so the curated index works on the deployed site while staying a
# normal Markdown document in the repository.
#
# Usage: tools/stage_examples_site.sh <output-directory>

set -euo pipefail

OUT="${1:?usage: tools/stage_examples_site.sh <output-directory>}"
SRC="examples"
SITE_ROOT="https://ovvo-financial.github.io/NNS"

command -v pandoc >/dev/null || { echo "pandoc is required" >&2; exit 1; }

mkdir -p "$OUT"

rsync -a --delete \
  --include='*.html' --include='*.pdf' \
  --include='*.png' --include='*.jpg' --include='*.css' \
  --exclude='*' "$SRC/" "$OUT/"

# Shared navigation strip injected above every rendered Markdown page.
NAV_FILE="$(mktemp)"
cat > "$NAV_FILE" <<NAV
<nav class="nns-site-nav">
  <a href="$SITE_ROOT/">NNS home</a>
  <a href="$SITE_ROOT/examples/">Examples</a>
  <a href="$SITE_ROOT/book/">Book</a>
  <a href="https://github.com/OVVO-Financial/NNS">GitHub</a>
</nav>
NAV
trap 'rm -f "$NAV_FILE"' EXIT

render_markdown() {
  local source="$1" target="$2" title="$3"
  # Point relative Markdown links at their rendered pages. Links containing
  # ':' (absolute URLs) are left untouched.
  sed -E 's/\(([^):]*)\.md\)/(\1.html)/g' "$source" | pandoc \
    --from gfm \
    --to html5 \
    --standalone \
    --css=pages.css \
    --metadata pagetitle="$title" \
    --include-before-body="$NAV_FILE" \
    --output "$target"
  echo "rendered: $target"
}

while IFS= read -r -d '' md; do
  name="$(basename "$md" .md)"
  case "$name" in
    README|index) continue ;;
  esac
  render_markdown "$md" "$OUT/$name.html" "NNS examples: $name"
done < <(find "$SRC" -maxdepth 1 -name '*.md' -print0)

render_markdown "$SRC/index.md" "$OUT/index.html" "NNS Examples"

echo "staged $(find "$OUT" -type f | wc -l) files into $OUT"

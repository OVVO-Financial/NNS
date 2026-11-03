#!/usr/bin/env bash
# Stage the examples/ directory as a browsable GitHub Pages section.
#
# Copies published artifacts (html, pdf, images), marks historical HTML outputs,
# renders
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
command -v python3 >/dev/null || { echo "python3 is required" >&2; exit 1; }

mkdir -p "$OUT"

rsync -a --delete \
  --include='*/' --include='*.html' --include='*.pdf' \
  --include='*.R' --include='*.Rmd' \
  --include='*.png' --include='*.jpg' --include='*.css' \
  --exclude='*' "$SRC/" "$OUT/"

# Preserve the archived source files. Add a notice only to their staged copies
# until each affected numerical example has a current reproducible build.
python3 - "$OUT" "$SITE_ROOT" <<'PY'
import re
import sys
from pathlib import Path

output, site_root = Path(sys.argv[1]), sys.argv[2]
historical_pages = (
    "Causal_Inference_Amongst_Macroeconomic_Variables_Using_NNS.html",
    "NNSvignette_Overview.html",
    "Sunspots_example.html",
    "VAR_example.html",
    "prophet_NNS_comparison.html",
    "tides.html",
)
notice = (
    '<aside class="nns-historical-output" role="note" '
    'style="padding:1em;margin:1em;border:2px solid #9a6700;'
    'background:#fff8c5;color:#24292f">'
    '<strong>Historical output.</strong> This archived example contains '
    'results from earlier NNS implementations. Its numerical outputs have '
    'not been regenerated after the November 2026 corrections. See the '
    '<a href="' + site_root + '/examples/example-corrections.html'
    '#historical-html-output-regeneration">current corrections and '
    'regeneration checklist</a> and the current package vignettes before '
    'using these results.</aside>\n'
).encode("ascii")
for name in historical_pages:
    path = output / name
    if not path.is_file():
        continue
    content, count = re.subn(rb"(<body\b[^>]*>)",
                             lambda match: match[0] + b"\n" + notice,
                             path.read_bytes(), count=1, flags=re.IGNORECASE)
    if count != 1:
        raise SystemExit("Cannot place historical-output notice in " + name)
    path.write_bytes(content)
PY

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
    --mathjax=https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-mml-chtml.js \
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

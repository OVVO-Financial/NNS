# Bookdown Assembly Workspace

Yes — this repository can be assembled into a **bookdown** project.

This folder is a starter workspace that:

- Defines basic bookdown config files (`_bookdown.yml`, `_output.yml`, `index.Rmd`)
- Provides a chapter-generation script that converts existing chapter `.txt` files into chapter `.Rmd` files in build order

## Quick start

1. Install dependencies in R:

   ```r
   install.packages(c("bookdown", "rmarkdown"))
   ```

2. Generate chapter Rmd files from source text:

   ```bash
   bash scripts/generate_rmd_chapters.sh
   ```

3. Render the book from this folder:

   ```bash
   Rscript -e "bookdown::render_book('index.Rmd', 'bookdown::gitbook')"
   ```

## Notes

- Source chapters are read from `../Chapter <n>/NNS_chapter_<n>.txt`.
- Generated files are written in this folder as `chapter-XX.Rmd`.
- Re-run the generator whenever chapter text changes.

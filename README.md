# Direct download link:
[Course package](https://github.com/douglasadamoski/rcourse/releases/download/v20260921/RCourse_20260921.zip)

# Introduction to R

A one-day, hands-on introduction to R for people who work with biological data
and have never programmed before.

Taught by **Douglas Adamoski Meira**. Current edition: **21 September 2026**.

Everything used on the day lives in this repository: the slides, the script we
type together, and the example data.

---

## Who this is for

You if you analyze experiments in Excel or GraphPad and keep thinking there
must be a faster way. No prior programming of any kind is assumed. We start at
`1 + 1` and finish with volcano plots, heatmaps and a qPCR analysis.

By the end of the day you will be able to:

- read your own CSV and Excel files into R;
- build publication-quality figures with `ggplot2` and save them as vectors;
- write a loop that produces one figure per gene instead of clicking 30 times;
- run a delta-delta Ct qPCR analysis reproducibly;
- read someone else's R script well enough to adapt it.

---

## Before you arrive

Install both, in this order. They are free, and they are two separate programs:
R is the engine, RStudio is the dashboard you actually sit in front of.

1. **R** — <https://cran.r-project.org/>
2. **RStudio Desktop** — <https://posit.co/download/rstudio-desktop/>

Then download the course files:

1. Download **`RCourse_20260921.zip`** from the
   [latest release](https://github.com/douglasadamoski/rcourse/releases/latest).
2. Unzip it. It creates a folder with everything inside, already together.
3. Open `RCourse_20260921.R` in RStudio.
4. Tell R where that folder is: **Session → Set Working Directory → To Source
   File Location**. Skipping this step is the single most common reason the
   script cannot find the data files.

The script installs the R packages it needs in its first section. That takes a
few minutes, so if you can run those lines before the course starts, please do.

---

## What is in this repository

| File | What it is |
|---|---|
| `RCourse_20260921.R` | The script we work through, line by line. Heavily commented — it is meant to be readable months later. |
| `RCourse_20260921.pptx` | The slides. |
| `example_gene_table.csv` | A small RNA-seq-like result table: 30 genes, 3 controls, 3 treated, fold change and FDR. |
| `example_gene_table.xlsx` | The same table as an Excel file, to practice reading both formats. |
| `MyCtValues.xlsx` | Raw qPCR Ct values in the long format instruments export. |
| `archive/2025-09-29/` | The previous edition, kept exactly as taught. |

The one-file download for each edition is attached to its
[release](https://github.com/douglasadamoski/rcourse/releases), not committed here:
`RCourse_20260921.zip` and `RCourse_20250929.zip`.

---

## The day

| Part | Topic | Where in the script |
|---|---|---|
| 1 | Why learn R? Reproducibility, automation, figures | — |
| 2 | R and RStudio: the panes, scripts, projects, the working directory | *Running commands* |
| 3 | Basics: variables, types, vectors, indexing, functions | *Assigning variables* → *Calling positions* |
| 4 | Data in and out: CSV, Excel, data frames, `ggplot2` | *Libraries* → *Adding information to the table* |
| 5 | Automating: `for` loops, one figure per gene, qPCR delta-delta Ct | *What about plotting every gene?*, *qPCR example* |
| 6 | Omics figures: volcano plots and heatmaps | *Volcano plots*, *Heatmaps* |
| 7 | Sharing your work: Shiny | *Sharing your work: Shiny* |
| 8 | Wrap-up, where to go next, questions | *Where to go next* |

---

## Packages used

Installed from **[CRAN](https://cran.r-project.org/)**:

```r
install.packages("ggplot2")    # plotting
install.packages("openxlsx")   # read and write Excel files
install.packages("ggpubr")     # statistics annotations on ggplots
install.packages("pheatmap")   # heatmaps
install.packages("shiny")      # interactive web apps
```

Installed from **[Bioconductor](https://www.bioconductor.org/)**, the
repository for biology and omics packages:

```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")   # volcano plots
BiocManager::install("ddCt")              # qPCR delta-delta Ct
```

Package documentation:
[ggplot2](https://ggplot2.tidyverse.org/) ·
[openxlsx](https://ycphs.github.io/openxlsx/) ·
[ggpubr](https://rpkgs.datanovia.com/ggpubr/) ·
[pheatmap](https://cran.r-project.org/package=pheatmap) ·
[EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano/) ·
[ddCt](https://bioconductor.org/packages/ddCt/) ·
[shiny](https://shiny.posit.co/)

---

## Where to go next

**Learning R**
- [R for Data Science](https://r4ds.hadley.nz/) — the free book to read after this course
- [Posit cheatsheets](https://posit.co/resources/cheatsheets/) — one page per package, worth printing
- [Posit Primers](https://posit.cloud/learn/primers) — short interactive exercises
- [The R Manuals](https://cran.r-project.org/manuals.html) — the official reference

**Finding packages**
- [CRAN task views](https://cran.r-project.org/web/views/) — curated package lists by topic
- [Bioconductor workflows](https://bioconductor.org/packages/release/BiocViews.html#___Workflow) — worked end-to-end analyses
- [Bioconductor support forum](https://support.bioconductor.org/) — where to ask omics questions

**Going further**
- [Shiny gallery](https://shiny.posit.co/r/gallery/) — interactive apps, with source code
- [Quarto](https://quarto.org/) — reports and papers with the code and figures inside
- [The tidyverse](https://www.tidyverse.org/) — a different, very popular style of R

**When something breaks**

Copy the exact error message into a search engine. Someone has had it before.
If that fails, [Stack Overflow's r tag](https://stackoverflow.com/questions/tagged/r)
and [Posit Community](https://forum.posit.co/) both answer beginner questions
without fuss.

---

## Versions

| Edition | Tag | Notes |
|---|---|---|
| 21 September 2026 | [`v20260921`](https://github.com/douglasadamoski/rcourse/releases/tag/v20260921) | `RCourse_20260921.zip` · Added Shiny; rewrote the script comments; fixed the bugs that stopped the 2025 script running end to end. |
| 29 September 2025 | [`v20250929`](https://github.com/douglasadamoski/rcourse/releases/tag/v20250929) | `RCourse_20250929.zip` · First edition. Files kept under `archive/2025-09-29/`. |

---

## Reuse and contact

These materials are shared so that participants can keep them and reread them.
If you would like to reuse or adapt them for your own teaching, please get in
touch first — I am happy to say yes, I would just like to know.

Questions, corrections and typo reports are welcome as
[issues](https://github.com/douglasadamoski/rcourse/issues).

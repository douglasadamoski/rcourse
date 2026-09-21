# =============================================================================
# Introduction to R
# Douglas Adamoski Meira - 21.09.2026
# https://github.com/douglasadamoski/rcourse/
#
# How to use this script:
#   Run it line by line. Put the cursor on a line and press Ctrl+Enter
#   (Cmd+Enter on a Mac). Read the output in the Console before moving on.
#   Nothing here needs to be run all at once.
#
# What you need in your working directory:
#   example_gene_table.csv    a small RNA-seq-like result table
#   example_gene_table.xlsx   the same table as an Excel file
#   MyCtValues.xlsx           raw qPCR Ct values
#
# Everything is in RCourse_20260921.zip, on the Releases page of the repo
# linked above. Unzip it, then point R at the folder it creates
# (Session > Set Working Directory > To Source File Location).
# =============================================================================



# -----------------------------------------------------------------------------
# Running commands
# -----------------------------------------------------------------------------
# The Console is a calculator that never gets tired. Type an expression,
# press Enter, and R prints the answer straight back.

1 + 1

10 * 2

20 / 4



# -----------------------------------------------------------------------------
# Assigning variables
# -----------------------------------------------------------------------------
# "<-" stores a value under a name so you can reuse it later.
# Read it out loud as "banana gets 2".

banana <- 2

# Typing the name on its own asks R to print what is stored in it.
banana

# Once a name holds a number, it behaves exactly like that number.
2 * 3
banana * 3

# "=" also assigns. Both work, but "<-" is the R convention and it is what
# you will see in almost every script and textbook.
youcanchooseyourownname = 5

# You can print it anytime
youcanchooseyourownname

# Names are yours to choose. Choose ones your future self will understand.
2 * 5
banana * youcanchooseyourownname



# -----------------------------------------------------------------------------
# Variable types
# -----------------------------------------------------------------------------
# Every object in R has a type. class() tells you which one.

class(banana)

# Text goes inside quotes. Note that the quotes matter, not the word:
# "banana" here is a piece of text that happens to spell the same as the
# variable banana above. They are two completely different things.
pineapple <- "banana"

pineapple

class(pineapple)



# -----------------------------------------------------------------------------
# Vectors
# -----------------------------------------------------------------------------
# A vector is an ordered collection of values of the SAME type.
# c() means "combine".

myfirstvector <- c(1, 2, 3)
myfirstvector

mysecondvector <- c("banana", "pineapple")
mysecondvector

# Careful: this combines the VALUES stored in banana (2) and pineapple
# ("banana"). A vector cannot mix types, so R quietly converts the number
# into text. This silent conversion is a classic source of confusion.
mythirdvector <- c(banana, pineapple)
mythirdvector

class(mythirdvector)



# -----------------------------------------------------------------------------
# Calling positions
# -----------------------------------------------------------------------------
# Square brackets pick elements out of a vector.
# R counts the way people do: the first element is 1, not 0.

mysecondvector[1]

mysecondvector[2]

# There is no third element, so R returns NA ("not available") rather than
# an error. Getting NA where you expected a value usually means you asked
# for something that is not there.
mysecondvector[3]

# The position can itself be a variable. banana holds 2, so this is the
# same as mysecondvector[2].
mysecondvector[banana]


# Functions can work on a whole vector at once, with no loop needed.
# This is one of the things R is genuinely good at.

myfirstvector

mean(myfirstvector)

sum(myfirstvector)

length(myfirstvector)

class(myfirstvector)



# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
# Base R is deliberately small. Packages (also called libraries) add the
# rest: plotting, Excel files, statistics, genomics.
#
# CRAN is the general-purpose package repository:
# https://cran.r-project.org/

# Run these installs ONCE per computer. After that, skip straight to the
# library() calls further down.
# Pay attention to the quotes: you install a package by its NAME as text.
install.packages("ggplot2")
install.packages("openxlsx")
install.packages("ggpubr")
install.packages("pheatmap")


# Bioconductor is the repository for biology and omics packages.
# It has its own installer, BiocManager, which you install from CRAN first.
# https://www.bioconductor.org/
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("ddCt")


# Installing puts a package on your disk. library() loads it into the
# current session. You install once, but you load in every new session.
library(ggplot2)          # plotting
library(openxlsx)         # read and write Excel files
library(ggpubr)           # statistics annotations on ggplots
library(EnhancedVolcano)  # volcano plots
library(pheatmap)         # heatmaps
library(ddCt)             # qPCR delta-delta Ct



# -----------------------------------------------------------------------------
# Reading data in
# -----------------------------------------------------------------------------
# First, let's check where R currently is. Every file name you type without
# a full path is interpreted relative to this folder.
getwd()

# The working directory is the single most common source of
# "cannot open file" errors. You can give a full path
# ("C:/Users/me/course/example_gene_table.csv") or a relative one
# ("example_gene_table.csv"), but the relative one only works if R is
# already in the right folder.

myGeneTable_fromCSV <- read.csv(file = "example_gene_table.csv",
                                header = TRUE,
                                row.names = 1)
myGeneTable_fromCSV

# How do I get help on a function? Put a question mark in front of it.
# The help page opens in the Help pane and lists every argument.
?read.csv

# Now let's read it again, the "harder" way: without row.names, so the gene
# names stay as an ordinary column instead of becoming row names. We will
# fix that by hand below, which is a good way to see what row.names = 1
# was doing for us.
myGeneTable_fromCSV <- read.csv("example_gene_table.csv")
myGeneTable_fromCSV

# A data.frame has names along both edges.
rownames(myGeneTable_fromCSV)   # just 1, 2, 3, ... for now
colnames(myGeneTable_fromCSV)   # Gene, Control_1, ... FDR

# The $ operator pulls out one column by name.
myGeneTable_fromCSV$Gene

# Square brackets on a data.frame take two positions: [rows, columns].
# Leaving one side empty means "all of them", so this is "all rows,
# column Gene" - the same thing the $ did.
myGeneTable_fromCSV[, "Gene"]

# The column name can come from a variable, which is what makes loops and
# functions possible later on.
mycolumnnow <- "Gene"
mycolumnnow
myGeneTable_fromCSV[, mycolumnnow]

# Assign the gene names as row names.
rownames(myGeneTable_fromCSV) <- myGeneTable_fromCSV[, "Gene"]

# See the row names now.
rownames(myGeneTable_fromCSV)

# Check the whole table again.
myGeneTable_fromCSV

# And that is the point of all this: you can now ask for a row by position
# OR by gene name. The second one is much harder to get wrong.
myGeneTable_fromCSV[6, ]
myGeneTable_fromCSV["BRCA1", ]


# The "easier" way, if your data is already in Excel. openxlsx reads a
# named sheet directly, and can take the row names for you.
myGeneTable_fromXLSX <- read.xlsx(xlsxFile = "example_gene_table.xlsx",
                                  sheet = "MyNiceTable",
                                  startRow = 1,
                                  colNames = TRUE,
                                  rowNames = TRUE)

myGeneTable_fromXLSX



# -----------------------------------------------------------------------------
# Let's do some plots!
# -----------------------------------------------------------------------------
# Before plotting, we need something to plot. Let's summarize the three
# replicates of each condition into a single mean per gene.

# How do we grab only the control columns? By name, in a vector.
myGeneTable_fromXLSX[, c("Control_1", "Control_2", "Control_3")]

# The same thing, but only the first row.
myGeneTable_fromXLSX[1, c("Control_1", "Control_2", "Control_3")]

# rowMeans() averages across the columns, one value per row.
# On a single row, that is one number.
rowMeans(myGeneTable_fromXLSX[1, c("Control_1", "Control_2", "Control_3")])

# On the whole table, it is one number per gene.
rowMeans(myGeneTable_fromXLSX[, c("Control_1", "Control_2", "Control_3")])

# Assigning into a column name that does not exist yet CREATES it.
myGeneTable_fromXLSX[, "Control_mean"] <- rowMeans(myGeneTable_fromXLSX[, c("Control_1", "Control_2", "Control_3")])

# Now the same for the treated samples.
myGeneTable_fromXLSX[, "Treated_mean"] <- rowMeans(myGeneTable_fromXLSX[, c("Treated_1", "Treated_2", "Treated_3")])

# Base R can plot straight away. It is quick, it is ugly, and it is
# perfectly fine for a first look at your own data.
plot(x = myGeneTable_fromXLSX$Control_mean,
     y = myGeneTable_fromXLSX$Treated_mean,
     xlab = "Control (mean)",
     ylab = "Treated (mean)",
     main = "Per-gene means (base R)",
     pch = 19)



# -----------------------------------------------------------------------------
# Adding information to the table
# -----------------------------------------------------------------------------
# We want to color the points by significance, so we need a column that
# says whether each gene is significant. Let's build it step by step.

# A comparison on a whole column gives a TRUE/FALSE for every gene.
myGeneTable_fromXLSX$FDR < 0.05

# ifelse() turns that TRUE/FALSE vector into something more readable.
ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No")

# A factor is R's type for a category. It remembers which values are
# allowed and, importantly for plots, in which ORDER they should appear.
factor(ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No"))

# Setting levels explicitly is what puts "Yes" before "No" in the legend
# instead of the default alphabetical order.
myGeneTable_fromXLSX[, "Signif"] <- factor(ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No"),
                                           levels = c("Yes", "No"))


# A ggplot is built by adding layers with "+".
#   ggplot(data, aes(...))  says which table and which columns map to which
#                           part of the picture
#   geom_point()            says "draw these as points"
ggplot(myGeneTable_fromXLSX, aes(x = Control_mean, y = Treated_mean)) +
  geom_point(size = 2.5) +
  labs(title = "Per-gene means",
       x = "Control (mean)", y = "Treated (mean)") +
  theme_minimal()

# The theme controls everything that is not the data itself. Swapping one
# line changes the whole look.
ggplot(myGeneTable_fromXLSX, aes(x = Control_mean, y = Treated_mean)) +
  geom_point(size = 2.5) +
  labs(title = "Per-gene means",
       x = "Control (mean)", y = "Treated (mean)") +
  theme_classic()


# Now let's make it more informative, one layer at a time.
ggplot(myGeneTable_fromXLSX, aes(x = Control_mean, y = Treated_mean, color = Signif)) +
  geom_point(size = 2.5) +
  # Add a regression line (lm = linear model), without the confidence band
  geom_smooth(method = "lm", se = FALSE) +
  # Force the same scale on both axes, so "no change" really is the diagonal
  coord_equal() +
  # Everything the reader needs to know, in words
  labs(title = "Per-gene means with regression",
       x = "Control (mean)", y = "Treated (mean)", color = "FDR < 0.05") +
  theme_minimal() +
  # Pick the colors yourself, by category name. Colors can be English
  # names or hex codes.
  scale_color_manual(values = c("Yes" = "purple",
                                "No"  = "#de2d26"))


# A plot can also be stored in a variable instead of being drawn straight
# away. Nothing appears on screen when you run this - that is expected.
p <- ggplot(myGeneTable_fromXLSX, aes(x = Control_mean, y = Treated_mean, color = Signif)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  coord_equal() +
  labs(title = "Per-gene means with regression",
       x = "Control (mean)", y = "Treated (mean)", color = "FDR < 0.05") +
  theme_minimal() +
  scale_color_manual(values = c("Yes" = "purple",
                                "No"  = "#de2d26"))

# Typing the name draws it.
p

# print() does the same thing explicitly. Inside a loop or a function,
# typing the name is NOT enough - only print() actually draws. Get into
# the habit now and the loop further down will just work.
print(p)


# Saving to a file: open a device, print into it, close the device.
# Everything drawn between svg() and dev.off() goes into the file and not
# to the screen. Forgetting dev.off() leaves the file empty or locked.
# SVG is a vector format, so it stays sharp at any size and can still be
# edited in Inkscape or Illustrator afterwards. Use it for figures.
svg(filename = "MyFirstPlot.svg",
    width = 4,
    height = 4)
print(p)
dev.off()



# -----------------------------------------------------------------------------
# What about a barplot?
# -----------------------------------------------------------------------------
# A scatter plot shows all genes at once. Sometimes you want the six
# individual replicates of a single gene, the way a qPCR figure looks.

# Let's look at our table again.
myGeneTable_fromXLSX

# Pick one gene to work with.
gene_of_interest <- "TP53"   # <-- change as you wish


# One row, six columns: the six measurements for this gene.
TemporaryTable <- myGeneTable_fromXLSX[gene_of_interest,
                                       c("Control_1", "Control_2", "Control_3",
                                         "Treated_1", "Treated_2", "Treated_3")]
TemporaryTable

# ggplot wants one ROW per observation, not one column per observation.
# t() transposes the table, so the six samples become six rows.
TemporaryTable_transposed <- data.frame(t(TemporaryTable))
TemporaryTable_transposed

# After transposing, the sample names live in the row names. Copy them
# into a real column, because ggplot can only see columns.
TemporaryTable_transposed$Sample <- rownames(TemporaryTable_transposed)
TemporaryTable_transposed


# We also need to know which samples are controls and which are treated.
# The blunt way: type it out.
TemporaryTable_transposed$Groups <- c("Control", "Control", "Control",
                                      "Treated", "Treated", "Treated")

# Is there a smarter way? Yes - the information is already in the names.
# Typing it out by hand breaks silently the moment the column order changes.

# grepl() asks "does this text match this pattern?" and returns TRUE/FALSE.
# "^Control" means "starts with Control" (^ anchors to the beginning).
TemporaryTable_transposed$Sample
grepl("^Control", TemporaryTable_transposed$Sample)

# Feed that into ifelse() and the groups build themselves.
?ifelse
ifelse(test = grepl("^Control", TemporaryTable_transposed$Sample),
       yes = "Control",
       no  = "Treated")

# Overwrite the hand-typed version with the derived one.
TemporaryTable_transposed$Groups <- ifelse(test = grepl("^Control", TemporaryTable_transposed$Sample),
                                           yes = "Control",
                                           no  = "Treated")


# The value column is currently named after the gene, which changes every
# time we change gene_of_interest. Rename it to something stable so the
# plotting code below never has to change.
colnames(TemporaryTable_transposed)

colnames(TemporaryTable_transposed)[1]

colnames(TemporaryTable_transposed)[1] <- "Counts"

colnames(TemporaryTable_transposed)[1]

colnames(TemporaryTable_transposed)

# paste() glues text together, so titles can mention the current gene.
paste("qPCR-like summary for", gene_of_interest)

# Now the plot. Read it one layer at a time:
#   geom_bar(stat = "summary")   bar height = mean of the group
#   stat_summary(...errorbar)    standard error whiskers on top
#   geom_jitter()                the individual replicates, nudged sideways
#                                so they do not overlap - always show them
#   stat_compare_means()         a t-test, printed on the plot
ggplot(TemporaryTable_transposed, aes(x = Groups, y = Counts)) +
  geom_bar(stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.9) +
  ggpubr::stat_compare_means(method = "t.test", label = "p.format") +
  labs(title = paste("qPCR-like summary for", gene_of_interest),
       y = "Expression / Abundance") +
  theme_classic()



# -----------------------------------------------------------------------------
# What about plotting every gene?
# -----------------------------------------------------------------------------
# Thirty genes, thirty plots. Doing that by hand is where mistakes and
# lost afternoons come from. This is exactly what a computer is for.

# Somewhere to put the files. showWarnings = FALSE stops R complaining
# when the folder already exists, so you can re-run the script safely.
dir.create("MyPlots", showWarnings = FALSE)


# A for loop takes each element of a vector in turn, puts it in a variable,
# and runs the block in { } once per element. Run this first to see the
# shape of it: it prints and builds file names, but draws nothing yet.
for (gene_of_interest in c("TP53", "MYC", "BRCA1")) {

  print(gene_of_interest)

  print(paste("MyPlots/", "MyGeneIs", gene_of_interest, ".svg", sep = ""))
}


# Now the real thing: exactly the code we wrote above for one gene, with
# the "pick a gene" line replaced by the loop, and the screen replaced by
# a file. This is the usual way a script grows - get it right once, then
# wrap it in a loop.

# Every gene name in the table.
rownames(myGeneTable_fromXLSX)

for (gene_of_interest in rownames(myGeneTable_fromXLSX)) {

  # Get the temporary table for this one gene
  TemporaryTable <- myGeneTable_fromXLSX[gene_of_interest,
                                         c("Control_1", "Control_2", "Control_3",
                                           "Treated_1", "Treated_2", "Treated_3")]

  # Transpose, so replicates become rows
  TemporaryTable_transposed <- data.frame(t(TemporaryTable))
  TemporaryTable_transposed$Sample <- rownames(TemporaryTable_transposed)

  # Derive the groups from the sample names
  TemporaryTable_transposed$Groups <- ifelse(test = grepl("^Control", TemporaryTable_transposed$Sample),
                                             yes = "Control",
                                             no  = "Treated")

  # Stable name for the value column
  colnames(TemporaryTable_transposed)[1] <- "Counts"

  # Build the plot and keep it in p
  p <- ggplot(TemporaryTable_transposed, aes(x = Groups, y = Counts)) +
    geom_bar(stat = "summary", fun = "mean", width = 0.6, alpha = 0.7) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
    geom_jitter(width = 0.08, size = 2, alpha = 0.9) +
    ggpubr::stat_compare_means(method = "t.test", label = "p.format", label.x.npc = "right") +
    labs(title = paste("qPCR-like summary for", gene_of_interest),
         y = "Expression / Abundance") +
    theme_classic()

  # One file per gene. Note the print(p) - without it the file comes out
  # empty, because inside a loop R does not auto-print.
  svg(file = paste("MyPlots/", "MyGeneIs", gene_of_interest, ".svg", sep = ""),
      width = 4,
      height = 4)
  print(p)
  dev.off()

}

# Have a look in the MyPlots folder. Thirty figures, one command.
list.files("MyPlots")



# -----------------------------------------------------------------------------
# Volcano plots
# -----------------------------------------------------------------------------
# A volcano plot puts effect size on x and significance on y, so the genes
# that are both large and reliable end up in the top corners.

# Genes we want labeled by name.
label_genes <- c("TP53", "MYC")


# Volcano plots use log2 fold change, not raw fold change: it makes
# "twice as much" and "half as much" symmetric around zero (+1 and -1).
myGeneTable_fromXLSX$log2FC <- log2(myGeneTable_fromXLSX$FoldChange)

# One wrinkle in this example table: two genes have an FDR written as
# exactly 0, which is really "smaller than the rounding in the file".
# The y-axis is -log10(FDR), and -log10(0) is infinity, which no plot can
# draw. The usual fix is to floor those values at something defensible -
# here, half of the smallest FDR we actually observed - and to say so in
# the figure legend. Never silently drop the genes instead.
smallest_nonzero_FDR <- min(myGeneTable_fromXLSX$FDR[myGeneTable_fromXLSX$FDR > 0])
smallest_nonzero_FDR

myGeneTable_fromXLSX$FDR_forPlot <- ifelse(myGeneTable_fromXLSX$FDR == 0,
                                           smallest_nonzero_FDR / 2,
                                           myGeneTable_fromXLSX$FDR)

# Lets save the plot into a variable, as sometimes it fails due to the size
p <- EnhancedVolcano(myGeneTable_fromXLSX,
                lab = rownames(myGeneTable_fromXLSX),
                x = "log2FC",
                y = "FDR_forPlot",      # FDR, with exact zeros floored
                pCutoff = 0.05,         # horizontal threshold line
                FCcutoff = 1.0,         # vertical lines at 2-fold up and down
                selectLab = label_genes,
                xlab = "log2(Fold Change) (Treated / Control)",
                ylab = "-log10(FDR)",
                title = "Volcano Plot (log2FC vs FDR)",
                subtitle = "Zero FDR values floored for display",
                legendPosition = "right",
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                boxedLabels = FALSE
)

# Now we can save the plot to the disk using a distinct strategy!
ggplot2::ggsave("volcano.png", p, width = 10, height = 8, dpi = 300)

# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------
# A heatmap shows the whole matrix at once and clusters rows and columns by
# similarity, so patterns across samples become visible.

# pheatmap can draw a colored strip above the columns telling the reader
# which sample is which. It expects a small data.frame whose ROW NAMES
# match the column names of the matrix being plotted.
annotation_col <- data.frame(
  Group = c("Control", "Control", "Control",
            "Treated", "Treated", "Treated")
)
rownames(annotation_col) <- c("Control_1", "Control_2", "Control_3",
                              "Treated_1", "Treated_2", "Treated_3")

# Check it. If these names do not match exactly, the strip silently
# disappears.
annotation_col

# scale = "row" centers and scales each gene separately. Without it, the
# few highly expressed genes wash out everything else and the heatmap
# shows you nothing but which genes are abundant.
pheatmap(myGeneTable_fromXLSX[, c("Control_1", "Control_2", "Control_3",
                                  "Treated_1", "Treated_2", "Treated_3")],
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Gene Expression Heatmap")



# -----------------------------------------------------------------------------
# qPCR example: delta-delta Ct
# -----------------------------------------------------------------------------
# Everything so far started from a table someone else had already
# processed. This section starts from raw machine output instead.

# Read the Ct values. This file is in "long" format - one row per well,
# with columns Sample, Detector, Ct and Platename - which is what most
# qPCR instruments export.
myCtValues <- read.xlsx(xlsxFile = "MyCtValues.xlsx",
                        startRow = 1,
                        colNames = TRUE,
                        rowNames = FALSE)

# Show the table.
myCtValues


# The delta-delta Ct method needs two decisions from you, and they are
# biological decisions, not technical ones:
#   housekeepingGene   the gene assumed not to change between conditions
#   calibrationSample  the sample everything else is expressed relative to
result <- ddCtExpression(InputFrame(myCtValues),
                         calibrationSample = "Sample2",
                         housekeepingGene = "Gene2")

# Pull the results out as a plain table.
out <- elist(result)

out

# Write it to Excel, so it can go straight into a report or to a
# collaborator who does not use R.
write.xlsx(x = out,
           file = "My_qPCR.xlsx")


# And a quick look at the same numbers as a bar chart with error bars.
br <- errBarchart(result)
print(br)



# -----------------------------------------------------------------------------
# Sharing your work: Shiny
# -----------------------------------------------------------------------------
# Everything above produces a figure YOU can regenerate. Shiny produces
# something your colleagues can use without touching R at all: a small web
# page with dropdowns and sliders, where the plotting code runs in the
# background as they click.
#
# It is worth the effort when you find yourself re-running the same script
# with one value changed, over and over, for someone else.

install.packages("shiny")   # run once
library(shiny)

# Shiny ships with a set of small, complete example apps. They are the
# fastest possible way to see one running - nothing to download, no files
# to create. This one is roughly 20 lines of code in total: a slider that
# controls the number of bins in a histogram.
runExample("01_hello")

# Close the app (the red stop sign in the Console) before running the next
# line - a running app holds the R session.

# Calling it with no arguments lists everything available:
runExample()

# The full list, in increasing order of complexity:
#   01_hello       slider controlling a histogram
#   02_text        printing text output
#   03_reactivity  how outputs update when inputs change
#   04_mpg         plotting from a built-in dataset
#   05_sliders     the different kinds of slider
#   06_tabsets     organizing output into tabs
#   07_widgets     the full catalog of input controls
#   09_upload      letting the user upload their own file
#   10_download    letting the user download a result
#   11_timer       output that updates on a schedule

# Read the source of any of them - this prints the folder they live in:
system.file("examples", "01_hello", package = "shiny")

# Now go to this place and open the file in RStudio!

# How would OUR bar plot become an app? Every Shiny app has the same two
# halves, and our code already contains both:
#
#   ui      what the user sees. Here: one dropdown listing
#           rownames(myGeneTable_fromXLSX), and a space for a plot.
#
#   server  what R does about it. Here: the ggplot block from the barplot
#           section, wrapped in renderPlot({ ... }), with
#           gene_of_interest replaced by input$gene.
#
# That is the whole idea: the loop we wrote above runs the plot 30 times
# in advance; a Shiny app runs it once, on demand, for whichever gene the
# user picked. Same plotting code either way.
#
# Where to go from here:
#   Gallery with source code   https://shiny.posit.co/r/gallery/
#   Step-by-step tutorial      https://shiny.posit.co/r/getstarted/
#   Free hosting for small apps  https://www.shinyapps.io/



# -----------------------------------------------------------------------------
# Where to go next
# -----------------------------------------------------------------------------
#   R for Data Science (free book)   https://r4ds.hadley.nz/
#   ggplot2 reference                https://ggplot2.tidyverse.org/reference/
#   Bioconductor workflows           https://bioconductor.org/packages/release/BiocViews.html#___Workflow
#   CRAN task views (by topic)       https://cran.r-project.org/web/views/
#   Posit Cheatsheets                https://posit.co/resources/cheatsheets/
#
# And the most useful habit of all: when something breaks, copy the exact
# error message into a search engine. Someone has already had it.

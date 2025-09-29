



# Running commands
1 + 1

10*2

20/4



# Assigning variables
banana <- 2

banana

2*3
banana*3

youcanchooseyourownname = 5

2*5
banana*youcanchooseyourownname



# Variable types
class(banana)

pineapple <- "banana"

pineapple

class(pineapple)



# Vectors
myfirstvector <- c(1,2,3)
myfirstvector


mysecondvector <- c("banana", "pineapple")
mysecondvector

mythirdvector <- c(banana, pineapple)
mythirdvector


# Calling positions
# In R we count as in the real world
# 1,2,3...

mysecondvector[1]

mysecondvector[2]

mysecondvector[3]

mysecondvector[banana]

# You can operate over vectors

myfirstvector

mean(myfirstvector)

sum(myfirstvector)

length(myfirstvector)

class(myfirstvector)



# Libraries
# They can extend R functions
# https://cran.r-project.org/


# Pay attention to the quotes
install.packages("ggplot2")
install.packages("openxlsx")
install.packages("ggpubr")
install.packages("reshape2")
install.packages("pheatmap")


# https://www.bioconductor.org/
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("ddCt")



# Now ggplot2 is part of your 
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(reshape2)
library(EnhancedVolcano)
library(pheatmap)
library(ddCt)


# Lets check where we are.
getwd()

# Working dir location is very important
# you can provide the full location or relative one.

myGeneTable_fromCSV <- read.csv(file="example_gene_table.csv",
                                  header=TRUE,
                                  row.names = 1)
myGeneTable_fromCSV

# How can I get quick help?
?read.csv

# You can go "harder"
myGeneTable_fromCSV <- read.csv("example_gene_table.csv")
myGeneTable_fromCSV

# row and col names
rownames(myGeneTable_fromCSV)
colnames(myGeneTable_fromCSV)

# Calling parts of my table
myGeneTable_fromCSV$Gene

# Calling parts of my table
myGeneTable_fromCSV[,"Gene"]

# What about mixing?
mycolumnnow <- "Gene"
mycolumnnow
myGeneTable_fromCSV[,mycolumnnow]

# Assign the rownames
rownames(myGeneTable_fromCSV) <- myGeneTable_fromCSV[,"Gene"]

# See the row names
rownames(myGeneTable_fromCSV)

# Check the whole table
myGeneTable_fromCSV

# get one line
myGeneTable_fromCSV[6,]
myGeneTable_fromCSV["BRCA1",]


# You can go "easier"
myGeneTable_fromXLSX <- read.xlsx(xlsxFile="example_gene_table.xlsx",
                                  sheet="MyNiceTable",
                                  startRow = 1,
                                  colNames = TRUE,
                                  rowNames = TRUE)

myGeneTable_fromXLSX



# Lets do some plots!


# Precompute per-gene means

# How to get only the control?
myGeneTable_fromXLSX[, c("Control_1","Control_2","Control_3")]

myGeneTable_fromXLSX[1, c("Control_1","Control_2","Control_3")]

rowMeans(myGeneTable_fromXLSX[1, c("Control_1","Control_2","Control_3")])

rowMeans(myGeneTable_fromXLSX[, c("Control_1","Control_2","Control_3")])

myGeneTable_fromXLSX[,"Control_mean"] <- rowMeans(myGeneTable_fromXLSX[, c("Control_1","Control_2","Control_3")])

# Now for the treated sample
myGeneTable_fromXLSX[,"Treated_mean"] <- rowMeans(myGeneTable_fromXLSX[, c("Treated_1","Treated_2","Treated_3")])

# Using base R
plot(myGeneTable_fromXLSX$Control_mean, myGeneTable_fromXLSX$Treated_mean,
     xlab="Control (mean)",
     ylab="Treated (mean)",
     main="Per-gene means (base R)",
     pch=19)



# Lets add more information to our table

#
myGeneTable_fromXLSX$FDR < 0.05

#
ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No")

# Factors!
factor(ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No"))

#
myGeneTable_fromXLSX[,"Signif"] <- factor(ifelse(myGeneTable_fromXLSX$FDR < 0.05, "Yes", "No"),
                     levels=c("Yes", "No"))


# Lets generate a simple plot
ggplot(myGeneTable_fromXLSX, aes(x=Control_mean, y=Treated_mean)) +
  geom_point(size=2.5) +
  labs(title="Per-gene means with regression",
       x="Control (mean)", y="Treated (mean)") +
  theme_minimal()

# Lets change the theme
ggplot(myGeneTable_fromXLSX, aes(x=Control_mean, y=Treated_mean)) +
  geom_point(size=2.5) +
  labs(title="Per-gene means with regression",
       x="Control (mean)", y="Treated (mean)") +
  theme_classic()




# Lets make it more complex
ggplot(myGeneTable_fromXLSX, aes(x=Control_mean, y=Treated_mean, color=Signif)) +
  geom_point(size=2.5) +
  # Add a regression line
  geom_smooth(method="lm", se=FALSE) +
  # Make the axis proportionate
  coord_equal() +
  #Label information
  labs(title="Per-gene means with regression",
       x="Control (mean)", y="Treated (mean)", color="FDR") +
  theme_minimal() +
  scale_color_manual(values = c("Yes" = "purple",
                               "No" = "#de2d26"))


# lets save it
p <- ggplot(myGeneTable_fromXLSX, aes(x=Control_mean, y=Treated_mean, color=Signif)) +
  geom_point(size=2.5) +
  # Add a regression line
  geom_smooth(method="lm", se=FALSE) +
  # Make the axis proportionate
  coord_equal() +
  #Label information
  labs(title="Per-gene means with regression",
       x="Control (mean)", y="Treated (mean)", color="FDR") +
  theme_minimal() +
  scale_color_manual(values = c("Yes" = "purple",
                                "No" = "#de2d26"))


# close the window
dev.off()

# get me the plot
p
# This print will be important in a near future, try to always use it!
print(p)


# Save it (and a good lesson about SVGs)
svg(filename="MyFirstPlot.svg",
    width=4,
    height=4)
# 
print(p)
dev.off()



# What about a barplot?

# Lets look again into our table
myGeneTable_fromXLSX

#
gene_of_interest <- "TP53"   # <-- change as you wish


TemporaryTable <- myGeneTable_fromXLSX[gene_of_interest, 
          c("Control_1","Control_2","Control_3",
            "Treated_1","Treated_2","Treated_3")]

# Melt requires a data.frame
TemporaryTable_transposed <- data.frame(t(row))           # transpose so replicates are rows
TemporaryTable_transposed$Sample <- rownames(row_df)


# What about our groups?
TemporaryTable_transposed$Groups <- c("Control", "Control", "Control",
                                      "Treated", "Treated", "Treated")

# Any smarter way?

# Regex codes
TemporaryTable_transposed$Sample
grepl("^Control", TemporaryTable_transposed$Sample)

# use one ifelse
?ifelse
ifelse(test=grepl("^Control",TemporaryTable_transposed$Sample),
       yes="Control",
       no="Treated")

# Lets just change one column to "Counts"
colnames(TemporaryTable_transposed)

colnames(TemporaryTable_transposed)[1]

colnames(TemporaryTable_transposed)[1] <- "Counts"

colnames(TemporaryTable_transposed)[1]

colnames(TemporaryTable_transposed)

# What about creating a new name for the axis?
paste("qPCR-like summary for", gene_of_interest)

# Lets plot (and move point by point)
ggplot(TemporaryTable_transposed, aes(x=Groups, y=Counts)) +
  geom_bar(stat="summary", fun="mean", width=0.6, alpha=0.7) +
  stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
  geom_jitter(width=0.08, size=2, alpha=0.9) +
  ggpubr::stat_compare_means(method="t.test", label="p.format") +
  labs(title=paste("qPCR-like summary for", gene_of_interest),
       y="Expression / Abundance") +
  theme_classic()




# What about plotting everyone?

# Lets make a folder
dir.create("MyPlots")


# For loops

for(gene_of_interest in c("Gene1", "Gene2", "Gene3")){
  
  print(gene_of_interest)
  
  paste("MyPlots/", "MyGeneIs", gene_of_interest, ".svg", sep="")
}



# lets put all together
# Just to test
gene_of_interest <- "TP53"

#
rownames(myGeneTable_fromXLSX)

#
for(gene_of_interest in rownames(myGeneTable_fromXLSX)   ){
  
  # Get the temporaty table
  TemporaryTable <- myGeneTable_fromXLSX[gene_of_interest, 
                                         c("Control_1","Control_2","Control_3",
                                           "Treated_1","Treated_2","Treated_3")]
  
  # Melt requires a data.frame
  TemporaryTable_transposed <- data.frame(t(row))           # transpose so replicates are rows
  TemporaryTable_transposed$Sample <- rownames(row_df)
  
  
  # What about our groups?
  # use one ifelse
  TemporaryTable_transposed$Groups <- ifelse(test=grepl("^Control",TemporaryTable_transposed$Sample),
         yes="Control",
         no="Treated")
  
  colnames(TemporaryTable_transposed)[1] <- "Counts"

  # Lets plot (and move point by point)
  p <- ggplot(TemporaryTable_transposed, aes(x=Groups, y=Counts)) +
    geom_bar(stat="summary", fun="mean", width=0.6, alpha=0.7) +
    stat_summary(fun.data=mean_se, geom="errorbar", width=0.2) +
    geom_jitter(width=0.08, size=2, alpha=0.9) +
    ggpubr::stat_compare_means(method="t.test", label="p.format", label.x.npc="right") +
    labs(title=paste("qPCR-like summary for", gene_of_interest),
         y="Expression / Abundance") +
    theme_classic()
  
  #
  svg(file=  paste("MyPlots/", "MyGeneIs", gene_of_interest, ".svg", sep=""),
      width=4,
      height=4)
  print(p)
  dev.off()
  
}




# Cool, lets do some volcano plots!
label_genes <- c("TP53","MYC")


# Compute log2 fold-change from provided FoldChange
myGeneTable_fromXLSX$log2FC <- log2(myGeneTable_fromXLSX$FoldChange)

#
EnhancedVolcano(myGeneTable_fromXLSX,
                lab = rownames(myGeneTable_fromXLSX),
                x = "log2FC",
                y = "FDR",              # using FDR on y-axis
                pCutoff = 0.05,         # threshold line at FDR < 0.05
                FCcutoff = 1.0,         # |log2FC| >= 1 (i.e., >= 2-fold)
                selectLab = label_genes,
                xlab = "log2(Fold Change) (Treated / Control)",
                ylab = "FDR",
                title = "Volcano Plot (log2FC vs FDR)",
                subtitle = "Highlighted selected genes",
                legendPosition = "right",
                max.overlaps = Inf,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                boxedLabels = FALSE
)


# Lets do some heatmaps!

#Lets add the groups to some annotation
myGeneTable_fromXLSX <- ifelse(test=grepl("^Control",TemporaryTable_transposed$Sample),
                               yes="Control",
                               no="Treated")

# Add a annotation table for us
annotation_col <- data.frame(
  Group = c("Control", "Control", "Control",
            "Treated", "Treated", "Treated")
)
rownames(annotation_col) <- c("Control_1", "Control_2", "Control_3",
                              "Treated_1", "Treated_2", "Treated_3")

# check it!
annotation_col

# Pretty heatmap
pheatmap(myGeneTable_fromXLSX[,c("Control_1", "Control_2", "Control_3",
                                 "Treated_1", "Treated_2", "Treated_3")],
         scale = "row",                      # center & scale each gene
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Gene Expression Heatmap")






# qPCR example
# Lets read our Cts
myCtValues <- read.xlsx(xlsxFile="MyCtValues.xlsx",
                                  startRow = 1,
                                  colNames = TRUE,
                                  rowNames = FALSE)
# Show the table
myCtValues



# Add information about samples and housekeepgin genes
result <- ddCtExpression(InputFrame(myCtValues),
                         calibrationSample="Sample2",
                         housekeepingGene="Gene2")

# Get the result
out <- elist(result)

#
out

# Save it
write.xlsx(x=out,
           file="My_qPCR.xlsx")



# Lets quickly check the results
br <- errBarchart(result)
print(br)


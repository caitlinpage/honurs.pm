#Installations of Bioconductor packages - I think everything is installed already? 
#actually with BiocManager I think I can download them just by loading them?? seems to be
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")

BiocManager::install("recount3")
BiocManager::install("snapcount")
BiocManager::install("megadepth")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")
BiocManager::install("derfinder")
BiocManager::install("limma")
BiocManager::install("GenomicRanges")
BiocManager::install("regionReport")
BiocManager::install("cluserProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("rtracklayer", force = TRUE)
BiocManager::install("GenomicFeatures", force = TRUE)
BiocManager::install("bumphunter", force = TRUE)
BiocManager::install("derfinderPlot", force = TRUE)
BiocManager::install("devtools")
BiocManager::install("EBSEA")
install.packages("pheatmap")

#See installed packages - left out a lot so not sure if I've installed everything I need
installed.packages()

#Load packages (recount3 tools)
#Honestly don't really understand what I need, will probably add
library(BiocManager)
library(recount3)
library(snapcount)
library(megadepth)
library(edgeR)
library(DESeq2)
library(GenomicRanges)
library(limma)
library(derfinder)
library(IRanges)
library(EBSEA)
library(SummarizedExperiment)

#Load packages (normal R tools/visualisation)
#Honestly not sure how many I need, will probably add to this
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(readr)
library(cowplot)
library(RColorBrewer)
library(pheatmap)


#DESeq
#example
BiocManager::install("airway")
library("airway")
data("airway")
se <- airway
ddsSE_example <- DESeqDataSet(se, design = ~ cell + dex)
ddsSE_example
ddsSE_ovary.exon <- DESeqDataSet(count_data = assay, colDat)

dim(ovary_abcb1.exon)
class(ovary_abcb1.exon)
metadata(ovary_abcb1.exon)
#this gives the big list of everything - this is the column names of the ranged summ exp
colData(ovary_abcb1.exon)
colnames(ovary_abcb1.exon)
rownames(ovary_abcb1.exon)
#this gives the counts this is good
#now just need to figure out how to view to aggregate for exons
assay(ovary_abcb1.exon)
assays(ovary_abcb1.exon)
head(ovary_abcb1.exon)
head(assay(ovary_abcb1.exon, 3))
#this gives something that looks good and interesting
##dataframe with 65 rows, 12 columns
rowData(ovary_abcb1.exon)
#mostly same info as above but this version is probably better
##GRanges object with 65 ranges and 12 metadata columns
rowRanges(ovary_abcb1.exon)[, "exon_idx"]



#EBSEA
#example 
data("exonCounts")
head(exonCounts)

#SNAPCOUNT

#Query coverage for gene:ABCB1 using TCGA
sb <- QueryBuilder(compilation = "tcga", regions = "abcb1")
abcb1.gene <- query_gene(sb)
#Dimensions of object(r,c)
dim(abcb1.gene)
#details
head(abcb1.gene)
#subfilter to Ovary tissue type
ovary_abcb1.gene <- query_gene(sb_ovary)
dim(ovary_abcb1.gene)
class(ovary_abcb1.gene)
ovary_abcb1.gene
assay(ovary_abcb1.gene)
colData(ovary_abcb1.gene)
colnames(colData(ovary_abcb1.gene))
#OK I can do recount2 paper stuff with the snapcount data yay!
colData(ovary_abcb1.gene)[,c("read_count_as_reported_by_sra", "reads_downloaded")]
summary(colData(ovary_abcb1.gene)$proportion_of_reads_reported_by_sra_downloaded)
colData(ovary_abcb1.gene)[, c("mapped_read_count", "paired_end", "avg_read_length")]
colData(ovary_abcb1.gene)[, "geo_accession"]
colData(ovary_abcb1.gene)[, "characteristics"]


#Query for junctions
abcb1.jx.all <- query_jx(sb)
dim(abcb1.jx.all)
head(abcb1.jx.all)

##Subfilter junction by tissue type: ovary
sb_ovary <- set_column_filters(sb, all == "Ovary")
ovary_abcb1.jx.all <- query_jx(sb_ovary)
dim(ovary_abcb1.jx.all)
head(ovary_abcb1.jx.all)
###vignette wants me to filter again by fully annotated - finds splice sites,
###but would that filter out the fusion? because it wouldn't have exon 1???

#Query for exons
abcb1.exon <- query_exon(sb)
dim(abcb1.exon)
head(abcb1.exon)
##subfilter by tissue type
ovary_abcb1.exon <- query_exon(sb_ovary)
dim(ovary_abcb1.exon)
colnames(ovary_abcb1.exon)

#Percent Spliced In (PSI) - yeah I don't think I can do this
#I don't understand the regions part because it's chromosome ids
#Tells us how often an exon is included in surrounding transcript
#With a bit of luck maybe we can see how often exon 1 is included?
##Make a new query
##Left
lq <- QueryBuilder(compilation = "tcga", regions = )


#Following steps in recount paper
#Downloading TCGAovary data gene level - 430 samples
OVproject.gene <- recount3::create_rse_manual(
  project = "OV",
  project_home = "data_sources/tcga",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

dim(colData(OVproject.gene))
colnames(colData(OVproject.gene))

str(OVproject.gene)

rowRanges(OVproject.gene)
assays(OVproject.gene)
head(OVproject.gene)
colData(OVproject.gene)[,c("recount_qc.aligned_reads%.chrm", "recount_qc.aligned_reads%.chrm", "recount_qc.aligned_reads%.chrm")]
#drug treatment
colData(OVproject.gene)[, 225:237]
#metastasis
colData(OVproject.gene)[, 288:290]

class(OVproject.gene)
assays(OVproject.gene)
assays(OVproject.gene)$raw_counts
rowRanges(OVproject.gene)[, "gene_id"]$ENSG00000085563

colData(OVproject.gene)[, 177]
colData(OVproject.gene)[, 188]
colData(OVproject.gene)[, 176]
#Can see treatment for each sample - most have been treated, yay!
colData(OVproject.gene)["ENSG00000085563.1", 225]
cisOV.gene <-OVproject.gene[, OVproject.gene$"tcga.cgc_case_drug_therapy_drug_name" == "Cisplatin"]
cisOV.gene[, "rownames"]
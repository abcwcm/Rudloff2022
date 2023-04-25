library(magrittr)
library(data.table)
library(DESeq2)
library(DiffBind)
library(sva)

decoderFile <- "atac_decoder_diffBind_plusPublic_justn2.csv"
decoder.data <- read.table(decoderFile,header=T,stringsAsFactors=F, sep=",")
decoder.data$SampleID = make.names(decoder.data$SampleID)
decoder.data$Condition = factor(decoder.data$Condition)
decoder.data$Condition <- relevel(decoder.data$Condition, ref = "N")
decoder.data = subset(decoder.data, (Condition %in% c("N", "L10d" ,"P10d","L1d", "P1d", "L5d", "P5d", "M")))
decoder.data = subset(decoder.data, batch == "m2023" | Condition == "M")
decoder.data$Condition <- relevel(decoder.data$Condition, ref = "N")


# Load the DiffBind object and extract raw read count matrix 
# https://support.bioconductor.org/p/67307/
DB <- readRDS("DB.count.plusPublic_n2.Rds")
DB <- dba.count(DB, peaks=NULL, score=DBA_SCORE_READS)
consensus_peaks <- dba.peakset(DB, peaks=NULL, bRetrieve=TRUE)



# create dds
raw_counts = as.data.frame(mcols(consensus_peaks))
row.names(raw_counts) = paste0(as.data.frame(consensus_peaks)$seqnames, ":",as.data.frame(consensus_peaks)$start,"-",as.data.frame(consensus_peaks)$end)
raw_counts = raw_counts[,decoder.data$SampleID]
table(decoder.data$SampleID == colnames(raw_counts))

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = decoder.data,
                              design= ~ Condition)
dds <- estimateSizeFactors(dds)
dds = DESeq(dds)
save(dds,file="parking_dds.RData")
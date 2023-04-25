# Load  packages
library(magrittr)
library(data.table)
library(DESeq2)
library(DiffBind)
library(sva)

# Read metadata
decoderFile <- "atac_decoder_diffBind_plusPublic_justn2.csv"
decoder.data <- read.table(decoderFile,header=T,stringsAsFactors=F, sep=",")
decoder.data$Condition = factor(decoder.data$Condition)
decoder.data$Condition <- relevel(decoder.data$Condition, ref = "N")
decoder.data$SampleID = make.names(decoder.data$SampleID)



# Load the DiffBind object and extract raw read count matrix 
# https://support.bioconductor.org/p/67307/
DB <- readRDS("DB.count.plusPublic_n2.Rds")
DB <- dba.count(DB, peaks=NULL, score=DBA_SCORE_READS)
consensus_peaks <- dba.peakset(DB, peaks=NULL, bRetrieve=TRUE)

# Create DESeq2 object
raw_counts = as.data.frame(mcols(consensus_peaks))
row.names(raw_counts) = paste0(as.data.frame(consensus_peaks)$seqnames, ":",as.data.frame(consensus_peaks)$start,"-",as.data.frame(consensus_peaks)$end)
raw_counts = raw_counts[,decoder.data$SampleID]
table(decoder.data$SampleID == colnames(raw_counts))

dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = decoder.data,
                              design= ~ Condition)
dds <- estimateSizeFactors(dds)
vsd <- DESeq2::vst(dds, blind=TRUE)

# Estimate 3 surrogate variables using sva and add them to design matrix
mod1 <- model.matrix(~Condition, decoder.data)
mod0 <- model.matrix(~1, decoder.data)
dat <- assay(vsd)
svseq <- sva(dat, mod1, mod0, n.sv = 3)
dds$SV1 = svseq$sv[,1]
dds$SV2 = svseq$sv[,2]
dds$SV3 = svseq$sv[,3]
design(dds) = ~SV1+SV2+SV3+Condition

# Perform differential expression analysis
dds <- DESeq(dds, fitType = "parametric")

# Remove batch effects from the transformed counts and save the DESeqDataSet and transformed counts objects
assay(vsd) =  limma::removeBatchEffect(assay(vsd), covariates = svseq$sv, design = mod1)
save(dds, vsd, file="sva_dds_vsd.RData")
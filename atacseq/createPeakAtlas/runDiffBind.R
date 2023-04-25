# cd  /home/paz2005/archive/2023_01_mary_philip/all_atac_old_new_public
library(DiffBind)
library(magrittr)
library(DESeq2)

setwd("/home/paz2005/archive/2023_01_mary_philip/all_atac_old_new_public")

DB <- dba(sampleSheet = "atac_decoder_diffBind_plusPublic_justn2.csv", peakCaller = "macs", peakFormat = "narrow", config=data.frame(AnalysisMethod=DBA_DESEQ2))
DB <- dba.blacklist(DB, blacklist=DBA_BLACKLIST_MM10, greylist=FALSE)
DB_consensus<-dba.peakset(DB,consensus=DBA_CONDITION, minOverlap=0.75)
DB_consensus <- dba(DB_consensus,mask=DB_consensus$masks$Consensus,minOverlap=1)
consensus_peaks <- dba.peakset(DB_consensus, bRetrieve=TRUE)
DB <- dba.count(DB,  bUseSummarizeOverlaps=T, bParallel=TRUE, peaks=consensus_peaks)
saveRDS(DB, "DB.count.plusPublic_n2.Rds")


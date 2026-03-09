setwd("/data/wangyy/CUT_Tag/reanalysis/12_diff_peaks/H3K27me3")

library(GenomicRanges)
blacklist <- rtracklayer::import("/exec/wangyy/software/ChIP_seq/RefGenome/mouse/blacklist/ENCFF547MET_fixed.bed")

library(DiffBind)
samples <- read.csv('H3K27me3.csv')
H3K27me3 <- dba(sampleSheet=samples, config=data.frame(fragmentSize=FALSE, 
                                                       doBlacklist=TRUE, doGreylist=FALSE))
H3K27me3 <- dba.count(H3K27me3, summits = FALSE)
H3K27me3 <- dba.normalize(H3K27me3)
H3K27me3 <- dba.contrast(H3K27me3, minMembers = 2, design = "~ Condition", 
                         contrast = c("Condition","Silica_WT","PBS_WT"))
# H3K27me3 <- dba.contrast(H3K27me3, minMembers = 2, design = "~ Condition", 
#                          contrast = c("Condition","Silica_CKO","Silica_WT"))
H3K27me3 <- dba.blacklist(H3K27me3, blacklist = blacklist)
H3K27me3 <- dba.analyze(H3K27me3)
dba.show(H3K27me3, bContrasts=TRUE)
H3K27me3_result <- dba.report(H3K27me3)
H3K27me3_result 

cor <- plot(H3K27me3)
# cor <- dba.plotHeatmap(H3K27me3)
pca <- dba.plotPCA(H3K27me3, attributes=DBA_CONDITION, label=DBA_ID)



library(ChIPseeker)
options(ChIPseeker.ignore_1st_exon = T)
options(ChIPseeker.ignore_1st_intron = T)
options(ChIPseeker.ignore_downstream = T)
options(ChIPseeker.ignore_promoter_subcategory = T)

library(TxDb.Mmusculus.UCSC.mm39.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene

# consensus_peaks <- dba.peakset(H3K27me3, bRetrieve=TRUE)
# peak_anno <- annotatePeak(consensus_peaks, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb = "org.Mm.eg.db")
peak_anno <- annotatePeak(H3K27me3_result, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb = "org.Mm.eg.db")

plotAnnoPie(peak_anno)
# plotAnnoBar(peak_anno)
# plotDistToTSS(peak_anno)

write.table(peak_anno, file = 'H3K27me3_Silica_WT_vs_PBS_WT.txt',  sep = '\t', quote = FALSE, row.names = FALSE)
# write.table(peak_anno, file = 'H3K27me3_Silica_CKO_vs_Silica_WT.txt',  sep = '\t', quote = FALSE, row.names = FALSE)

library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(ReactomePA)
library(ggplot2)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

Pfiles <- list.files(path = "/nfs/turbo/umms-thahoang/sherine/GSE181251/macs", pattern = "\\.narrowPeak$", full.names = TRUE, recursive = TRUE)
#peaks = lapply(Pfiles[1:2], readPeakFile) #Ignore NFI_input

peaks=GenomicRanges::GRangesList(Nfi_1=readPeakFile(Pfiles[[1]]),
                                Nfi_2=readPeakFile(Pfiles[[2]]))


#In case you wanna read one by one 
#nfi_1 <- toGRanges("macs/Nfi_1_peaks.narrowPeak", format="narrowPeak", header =FALSE)
#nfi_1 = peaks[1]  
#seqlevels(nfi_1) <- paste0("chr", seqlevels(nfi_1))

peaks <- lapply(peaks, function(gr) {
    seqlevels(gr) <- paste0("chr", seqlevels(gr))
    return(gr)
  })

#specific chr, remove parameter for all chromosomes 
figure_name = paste("Nfi_replicates", "Coverage.pdf", sep="_")
pdf(file =figure_name)
#specify a specific peak 
p <- covplot(peaks, weightCol="V5", chrs=c("chr17"), xlim=c(4.5e7, 4.7e7)) 
#plot all Chromosomes
#p <- covplot(peaks)
col <- c(Nfi_1='red', Nfi_2='green')
p + facet_grid(chr ~ .id) + scale_color_manual(values=col) + scale_fill_manual(values=col)
dev.off() 

promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter) 
names(tagMatrixList) <- c("Nfi_1", "Nfi_2") 

figure_name = paste("Nfi_replicates", "AvgProfile.pdf", sep="_")
pdf(file =figure_name)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off() 


figure_name = paste("Nfi_replicates", "HM.pdf", sep="_")
pdf(file =figure_name)
tagHeatmap(tagMatrixList)
dev.off() 

peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
names(peakAnnoList) <- c("Nfi_1", "Nfi_2")

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) <- c("Nfi_1", "Nfi_2")
comp  <- compareCluster(
  geneCluster   = genes,
  fun           = "enrichKEGG",
  organism      = "mmu",
  pvalueCutoff  = 0.01,
  pAdjustMethod = "BH")

figure_name = paste("Nfi_replicates", "KEGGpathways.pdf", sep="_")
pdf(file =figure_name)
dotplot(comp , showCategory = 10, title = "KEGG Pathway Enrichment Analysis") + theme(axis.text.x = element_text(size = 8)) 
dev.off() 

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) <- c("Nfi_1", "Nfi_2")
figure_name = paste("Nfi_replicates", "vennplot.pdf", sep="_")
pdf(file =figure_name)
vennplot(genes) 
dev.off()



figure_name = paste("Nfi_replicates", "DistToTSS.pdf", sep="_")
pdf(file =figure_name)
plotDistToTSS(peakAnnoList)
dev.off()

figure_name = paste("Nfi_replicates", "Anno.pdf", sep="_")
pdf(file =figure_name)
plotAnnoBar(peakAnnoList)
dev.off() 

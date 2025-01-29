library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(ReactomePA)
library(ggplot2)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 


#This merged narrow peaks comming from: bedtools intersect -a macs/Nfi_1_peaks.narrowPeak -b macs/Nfi_2_peaks.narrowPeak -wo > nfi_peaks.narrowPeak

Pfiles <- list.files(path = "/nfs/turbo/umms-thahoang/sherine/GSE181251/", pattern = "\\.narrowPeak$", full.names = TRUE, recursive = TRUE)

peaks=GenomicRanges::GRangesList(Nfi=readPeakFile(Pfiles[[5]])) #To pull the merged version and ignore the rest 
seqlevels(peaks) <- paste0("chr", seqlevels(peaks))


peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
comp  <- compareCluster(
  geneCluster   = genes,
  fun           = "enrichKEGG",
  organism      = "mmu",
  pvalueCutoff  = 0.01,
  pAdjustMethod = "BH")


comp_df <- as.data.frame(comp)
write.csv(comp_df, "compareClusterResult.csv", row.names = FALSE)


#Plot selected pathways
figure_name = paste("Nfi", "KEGGSelectedpathways.pdf", sep="_")
pdf(file =figure_name)
enriched_results <- comp@compareClusterResult
selected_pathway_names <- c("PI3K-Akt signaling pathway - Mus musculus (house mouse)", "Ras signaling pathway - Mus musculus (house mouse)")
selected_pathways <- enriched_results[enriched_results$Description %in% selected_pathway_names, ]
dotplot(comp, showCategory = selected_pathway_names)
dev.off()





figure_name = paste("Nfi", "KEGGpathways.pdf", sep="_")
pdf(file =figure_name)
dotplot(comp , showCategory = 10, title = "KEGG Pathway Enrichment Analysis") + theme(axis.text.x = element_text(size = 8)) 
dev.off() 

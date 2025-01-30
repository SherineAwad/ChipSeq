library(ChIPseeker)
library(clusterProfiler)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(ReactomePA)
library(ggplot2)


txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene 

args <- commandArgs(trailingOnly = TRUE)
peakfile <- args[1]
peaks=GenomicRanges::GRangesList(Nfi=readPeakFile(peakfile)) #To pull the merged version and ignore the rest 
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

#Remove the  - Mus musculus (house mouse suffix 
comp@compareClusterResult$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", comp@compareClusterResult$Description)


#Plot selected pathways
enriched_results <- comp@compareClusterResult
selected_pathway_names <- c("Rap1 signaling pathway", "Hippo signaling pathway",
"Signaling pathways regulating pluripotency of stem cells","Wnt signaling pathway","MAPK signaling pathway","PI3K-Akt signaling pathway","Ras signaling pathway","Notch signaling pathway","Cell cycle","Polycomb repressive complex","Hippo signaling pathway - multiple species")
selected_pathways <- enriched_results[enriched_results$Description %in% selected_pathway_names, ]

figure_name = paste("Nfi", "KEGGSelectedpathways.pdf", sep="_")
pdf(file =figure_name)
#dotplot(comp, showCategory = selected_pathway_names)
dotplot(comp,  x = "Count", color = "p.adjust",showCategory = selected_pathway_names)
dev.off()

figure_name = paste("Nfi", "piechart.pdf", sep="_")
pdf(file =figure_name)
plotAnnoPie(peakAnnoList$Nfi)
dev.off()


figure_name = paste("Nfi", "KEGGpathways.pdf", sep="_")
pdf(file =figure_name)
dotplot(comp , showCategory = 10, title = "KEGG Pathway Enrichment Analysis") + theme(axis.text.x = element_text(size = 8)) 
dev.off() 

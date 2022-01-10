enrichment <- function(geneList_, fileName, outputPath = F){
 
  # # library(UniProt.ws)
  # up <- UniProt.ws(taxId=9606)
  # select(up, geneList.uniprot, "ENTREZ_GENE")
  convert2Symbol <- function(x, map){
    l<-unlist(strsplit(x, split="/"))
    map <- select(org.Hs.eg.db, l, "SYMBOL", "ENTREZID")
    return(paste(map$SYMBOL, collapse = "/"))
  }
  
   map <- select(org.Hs.eg.db, geneList_, "ENTREZID", "SYMBOL")
  geneList <- na.omit(map$ENTREZID)
  
  ego <- enrichGO(gene           = geneList,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05)
  
  tmp <- data.frame(ego)
  genes <- sapply(tmp$geneID,convert2Symbol)
  tmp$geneSymbols <- genes
  write.csv(tmp , file = paste0("./functional_analysis/",fileName,"_GO_BP.csv"), row.names=FALSE)
  
  ego <- enrichGO(gene           = geneList,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05)
  
  tmp <- data.frame(ego)
  genes <- sapply(tmp$geneID,convert2Symbol)
  tmp$geneSymbols <- genes
  write.csv(tmp , file = paste0("./functional_analysis/",fileName,"_GO_MF.csv"), row.names=FALSE)
  
  ego <- enrichGO(gene          = geneList,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "CC",
                  pAdjustMethod = "fdr",
                  pvalueCutoff  = 0.05)
  
 
  
  
  tmp <- data.frame(ego)
  genes <- sapply(tmp$geneID,convert2Symbol)
  tmp$geneSymbols <- genes
  write.csv(tmp , file = paste0("./functional_analysis/",fileName,"_GO_CC.csv"), row.names=FALSE)
  
  

# -------------------------------------------------------------------------
  ekegg <- enrichKEGG(gene         = geneList,
                      #keyType      = 'uniprot',
                      pAdjustMethod = "fdr",
                      pvalueCutoff  = 0.05)
 
  tmp <- data.frame(ekegg)
  genes <- sapply(tmp$geneID,convert2Symbol)
  tmp$geneSymbols <- genes
  write.csv(tmp , file = paste0("./functional_analysis/",fileName,"_Kegg.csv") , row.names=FALSE)

# PathView ----------------------------------------------------------------
# Not needed at this stage
  pathViewOut <- function(keggID, geneEntrez){
    gene.data <- rep(1, length(geneEntrez))
    names(gene.data) <- geneEntrez
    msg <- capture.output(pathview::pathview(gene.data  = geneEntrez,
                                             pathway.id = keggID,
                                             limit = list(gene = 2, cpd = 1), bins = list(gene = 20, cpd= 10),
                                             low = list(gene = "green", cpd = "blue"),
                                             species    = "hsa"), type = "message")
    
    msg <- grep("image file", msg, value = T)
    filename <- sapply(strsplit(msg, " "), function(x) x[length(x)])
    img <- png::readPNG(filename)
    png::writePNG(img, target = paste0("./functional_analysis/",fileName,"_", filename))
    invisible(file.remove(filename))
  }  
  # if(outputPath & nrow(tmp)>0){
  #   N = 5
  #   if(N>nrow(tmp))N=nrow(tmp)
  #   for (i in 1:nrow(tmp)) {
  #     print(tmp$ID[i])
  #     pathViewOut(tmp$ID[i], geneEntrez = geneList)
  #   }
  # }
}
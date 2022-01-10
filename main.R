library(stringr)
library(gdata)
library(venn)



serum <- read.csv("./input/Serum.csv", header = T, fileEncoding = "UTF-8-BOM")
s <- as.data.frame(apply(serum, 2, function(x) sapply(x,function(y) str_split(y, "\\|")[[1]][2])))

plasma <- read.csv("./input/Plasma.csv", header = T, fileEncoding = "UTF-8-BOM")
p <- as.data.frame(apply(plasma, 2, function(x) sapply(x,function(y) str_split(y, "\\|")[[1]][2])))

urine <- read.csv("./input/Urine.csv", header = T, fileEncoding = "UTF-8-BOM")
u <- as.data.frame(apply(urine, 2, function(x) sapply(x,function(y) str_trim(y, side = "both"))))
u[u==""]<-NA

dat <- as.data.frame(cbindX(s,p,u))
pheno <- read.csv("./input/pheno.csv", header = T, fileEncoding = "UTF-8-BOM", row.names = 1)


# all.prot <- unique(na.omit(union(unique(unlist(s)),union(unique(unlist(p)),unique(unlist(u))))))
# all.prot <- unique(gsub("-[0-9]+", "", all.prot))
# map <- select(org.Hs.eg.db, all.prot, c("SYMBOL", "ENTREZID"), "UNIPROT")
# write.csv(map, "./input/map.na.csv", quote = F, row.names = F)
# Venn Diagram ------------------------------------------------------------


par(mfrow=c(2,3))


s.prot <- na.omit(unique(unlist(s)))
p.prot <- na.omit(unique(unlist(p)))
u.prot <- na.omit(unique(unlist(u)))

keep <- which(pheno$Source == "serum" & pheno$Type == "PCa")
s.pca <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Type == "Control")
s.ctl <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.pca, s.ctl), ilab=TRUE, zcolor = "style", snames = c("Serum PCa", "Serum Control"))

s.ctl.only <- setdiff(s.ctl, s.pca)


keep <- which(pheno$Source == "urine" & pheno$Type == "PCa")
u.pca <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "urine" & pheno$Type == "Control")
u.ctl <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(u.pca, u.ctl), ilab=TRUE, zcolor = "style", snames = c("Urine PCa", "Urine Control"))


keep <- which(pheno$Source == "urine" & pheno$Type == "PCa")
u.pca <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Type == "PCa")
s.pca <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.pca, u.pca), ilab=TRUE, zcolor = "style", snames = c("Serum PCa", "Urine PCa"))

keep <- which(pheno$Source == "urine" & pheno$Type == "Control")
u.ctl <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Type == "Control")
s.ctl <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.ctl, u.ctl), ilab=TRUE, zcolor = "style", snames = c("Serum Control", "Urine Control"))

venn::venn(list(s.prot, u.prot), ilab=TRUE, zcolor = "style", snames = c("Serum", "Urine"))

venn::venn(list(s.prot, p.prot, u.prot), ilab=TRUE, zcolor = "style", snames = c("Serum", "Plasma", "Urine"))

keep <- which(pheno$Source == "serum" & pheno$Type == "PCa")
s.pca <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Type == "Control")
s.ctl <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "urine" & pheno$Type == "PCa")
u.pca <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "urine" & pheno$Type == "Control")
u.ctl <- na.omit(unique(unlist(dat[,keep])))

venn::venn(list(s.pca, s.ctl, u.pca, u.ctl), ilab=TRUE, zcolor = "style", 
           snames = c("Serum PCa", "Serum Control", "Urine PCa", "Urine Control"))


# Patient level Venn ------------------------------------------------------
#patients in common: "100" "101" "102" "103" "105" "109"
keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "100")
u.100 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "100")
s.100 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.100, u.100), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 100", "Urine Patient 100"))

keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "101")
u.101 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "101")
s.101 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.101, u.101), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 101", "Urine Patient 101"))

keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "102")
u.102 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "102")
s.102 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.102, u.102), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 102", "Urine Patient 102"))

keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "103")
u.103 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "103")
s.103 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.103, u.103), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 103", "Urine Patient 103"))

keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "105")
u.105 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "105")
s.105 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.105, u.105), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 105", "Urine Patient 105"))

keep <- which(pheno$Source == "urine"  & pheno$Patient.ID == "109")
u.109 <- na.omit(unique(unlist(dat[,keep])))
keep <- which(pheno$Source == "serum" & pheno$Patient.ID == "109")
s.109 <- na.omit(unique(unlist(dat[,keep])))
venn::venn(list(s.109, u.109), ilab=TRUE, zcolor = "style", snames = c("Serum Patient 109", "Urine Patient 109"))


# Patient level Heatmap ---------------------------------------------------
common <- c("100", "101", "102", "103", "105", "109")

# Serum only proteins (not in Urine)
s.only <- setdiff(na.omit(unique(unlist(dat[,which(pheno$Source == "serum" & pheno$Type == "PCa" & pheno$Patient.ID %in% common)]))), 
        na.omit(unique(unlist(dat[,which(pheno$Source == "urine" & pheno$Type == "PCa" & pheno$Patient.ID %in% common)]))))
x <- matrix(nrow = length(s.only), ncol = 6)
rownames(x) <- s.only
colnames(x) <- common

for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    prot <- na.omit(unique(unlist(dat[, which(pheno$Patient.ID == colnames(x)[j])])))
    if(rownames(x)[i] %in% prot)
      x[i,j] = 1
    else
      x[i,j] = 0
  }
}
# 74 proteins are Serum specific, 37 appear in at least 2 patients

x.sort <- x[order(rowSums(x), decreasing = T),]

#x.sort <- x.sort[which(rowSums(x.sort)>1),] # to select those appear in at least 2 patients

library(pheatmap)
cols = colorRampPalette(c("white", "purple"))(10)
pheatmap(x.sort,cluster_rows = F, cluster_cols = F, color = cols, legend = F, fontsize_row = 6, main = "Serum-Specific Proteins (not found in Urine) in Patients")


# Urine only proteins (not in Serume)
u.only <- setdiff(na.omit(unique(unlist(dat[,which(pheno$Source == "urine" & pheno$Type == "PCa" & pheno$Patient.ID %in% common)]))), 
                  na.omit(unique(unlist(dat[,which(pheno$Source == "serum" & pheno$Type == "PCa" & pheno$Patient.ID %in% common)]))))
x <- matrix(nrow = length(u.only), ncol = 6)
rownames(x) <- u.only
colnames(x) <- common

for (i in 1:nrow(x)) {
  for (j in 1:ncol(x)) {
    prot <- na.omit(unique(unlist(dat[, which(pheno$Patient.ID == colnames(x)[j])])))
    if(rownames(x)[i] %in% prot)
      x[i,j] = 1
    else
      x[i,j] = 0
  }
}
# 35 proteins are urine specific, 7 appear in at least 2 patients

x.sort <- x[order(rowSums(x), decreasing = T),]

#x.sort <- x.sort[which(rowSums(x.sort)>1),]

library(pheatmap)
cols = colorRampPalette(c("white", "skyblue"))(10)
pheatmap(x.sort,cluster_rows = F, cluster_cols = F, color = cols, legend = F, fontsize_row = 6, main = "Urine-Specific Proteins (not found in Urine) in Patients")

# Functional analysis -----------------------------------------------------
library(enrichplot) #http://yulab-smu.top/clusterProfiler-book/chapter12.html
library(ggnewscale)
library(clusterProfiler)
library(org.Hs.eg.db)
library("pathview")
require(grid)
require(png)

source("./enrichment.R")
l <- list()
map <- read.csv("./input/map.csv",header = T, fileEncoding = "UTF-8-BOM", row.names = NULL)

g <- s.only
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Serum-Specific-6Patients",outputPath = T)
l[[1]] <- g

g <- u.only
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Urine-Specific-6Patients",outputPath = T)
l[[2]] <- g


# Run lines 65 - 75 first
g <- setdiff(setdiff(setdiff(s.pca,s.ctl),u.pca),u.ctl)
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Serum-Specific-99Proteins",outputPath = T)
l[[3]] <- g

g <- setdiff(setdiff(setdiff(u.pca,u.ctl),s.pca),s.ctl)
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Urine-Specific-31Proteins",outputPath = T)
l[[4]] <- g

g <- intersect(setdiff(s.pca,s.ctl), setdiff(u.pca,u.ctl))
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Serum-and-Urine-6Proteins",outputPath = T)
l[[5]] <- g

g <- s.ctl.only
g <- unique(gsub("-[0-9]+", "", g))
g <- map$SYMBOL[which(map$UNIPROT %in% g)]
enrichment(geneList_ = g , fileName = "Serum-Ctl-Only-23Proteins",outputPath = T)

# rrvgo visualisation ----------------------------------------------------
#Read: http://www.bioconductor.org/packages/release/bioc/vignettes/rrvgo/inst/doc/rrvgo.html
library(rrvgo)
library(ggplot2)
goVis <- function(fileName){
  f.bp <- paste0("./functional_analysis/", fileName, "_GO_BP.csv")
  go.bp <- read.csv(f.bp, fileEncoding = "UTF-8-BOM")
  rownames(go.bp) <- go.bp$ID
  simMatrix <- calculateSimMatrix(go.bp$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  
  scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  scatterPlot(simMatrix, reducedTerms)
  ggsave(paste0("./functional_analysis/scatterPlots/", fileName, "_GO_BP.pdf"), width = 10, height = 6)

  pdf(file = paste0("./functional_analysis/treemapPlots/", fileName, "_GO_BP.pdf"), width = 8, height = 8)
  treemapPlot(reducedTerms)
  dev.off()
  
  # MF
  f.mf <- paste0("./functional_analysis/", fileName, "_GO_MF.csv")
  go.mf <- read.csv(f.mf, fileEncoding = "UTF-8-BOM")
  rownames(go.mf) <- go.mf$ID
  simMatrix <- calculateSimMatrix(go.mf$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="MF",
                                  method="Rel")
  
  scores <- setNames(-log10(go.mf$qvalue), go.mf$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  scatterPlot(simMatrix, reducedTerms)
  ggsave(paste0("./functional_analysis/scatterPlots/", fileName, "_GO_MF.pdf"), width = 10, height = 6)
  
  pdf(file = paste0("./functional_analysis/treemapPlots/", fileName, "_GO_MF.pdf"), width = 8, height = 8)
  treemapPlot(reducedTerms)
  dev.off()
  
  # CC
  f.cc <- paste0("./functional_analysis/", fileName, "_GO_CC.csv")
  go.cc <- read.csv(f.cc, fileEncoding = "UTF-8-BOM")
  rownames(go.cc) <- go.cc$ID
  simMatrix <- calculateSimMatrix(go.cc$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="CC",
                                  method="Rel")
  
  scores <- setNames(-log10(go.cc$qvalue), go.cc$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.7,
                                  orgdb="org.Hs.eg.db")
  
  scatterPlot(simMatrix, reducedTerms)
  ggsave(paste0("./functional_analysis/scatterPlots/", fileName, "_GO_CC.pdf"), width = 10, height = 6)
  
  pdf(file = paste0("./functional_analysis/treemapPlots/", fileName, "_GO_CC.pdf"), width = 8, height = 8)
  treemapPlot(reducedTerms)
  dev.off()
  
}

goVis(fileName = "Serum-Ctl-Only-23Proteins")
goVis(fileName = "Serum-Specific-6Patients")
goVis(fileName = "Urine-Specific-6Patients")
goVis(fileName = "Serum-Specific-99Proteins")
goVis(fileName = "Urine-Specific-31Proteins")
goVis(fileName = "Serum-and-Urine-6Proteins")


# bar plots ----------------------------------------------------------------
fileName = "Serum-Ctl-Only-23Proteins"
f.cc <- paste0("./functional_analysis/", fileName, "_GO_CC.csv")
go.cc <- read.csv(f.cc, fileEncoding = "UTF-8-BOM")

ggplot(data=go.cc, aes(x=Description, y=-log10(p.adjust))) + geom_bar(stat="identity", fill = "blue")+ 
  geom_hline(yintercept = -log10(0.05))+coord_flip() + theme_bw()

# go.bp <- read.csv("./functional_analysis/Serum-Specific-6Patients_GO_BP.csv", fileEncoding = "UTF-8-BOM")
# rownames(go.bp) <- go.bp$ID
# simMatrix <- calculateSimMatrix(go.bp$ID,
#                                 orgdb="org.Hs.eg.db",
#                                 ont="BP",
#                                 method="Rel")
# 
# scores <- setNames(-log10(go.bp$qvalue), go.bp$ID)
# reducedTerms <- reduceSimMatrix(simMatrix,
#                                 scores,
#                                 threshold=0.7,
#                                 orgdb="org.Hs.eg.db")
# heatmapPlot(simMatrix,
#             reducedTerms,
#             annotateParent=TRUE,
#             annotationLabel="parentTerm",
#             fontsize=6)
# scatterPlot(simMatrix, reducedTerms)
# treemapPlot(reducedTerms)

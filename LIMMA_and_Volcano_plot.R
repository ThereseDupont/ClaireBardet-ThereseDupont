library(GEOquery)
library(limma)
library(data.table)
library(EnhancedVolcano)

list.ds <- c("GSE110224","GSE35279","GSE23878","GSE9348")
list.gsms <- c("0101010101010101010101010101010101",
               "1111111111111111111111111111111111111111111111111111111111111111110000011111111",
               "11111111111111111111111111111111111000000000000000000000000",
               "1111111111111111111111111111111111111111111111111111111111111111111111000000000000"
               )

i <- 5
ds <- list.ds[i]
gsms <- list.gsms[i]

# -------------------------------------------------------------------------

gset <- getGEO(ds, GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

fvarLabels(gset) <- make.names(fvarLabels(gset))


sml <- strsplit(gsms, split="")[[1]]

ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

gs <- factor(sml)
groups <- make.names(c("Normal","Cancer"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

fit <- lmFit(gset, design)  # fit linear model

cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT$logFC <- - tT$logFC

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
write.table(tT, file=paste0("/home/claire/Documents/ULM/L3/open_science/projet/output_limma/",ds,".top.table.csv"), row.names=F, sep=",")



##################################################################################################
#                                        Volcano Plot                                            #
##################################################################################################



res1 <- fread(paste0("/home/claire/Documents/ULM/L3/open_science/projet/output_limma/",ds,".top.table.csv"), header=TRUE)

png(paste0('/home/claire/Documents/ULM/L3/open_science/projet/output_limma/',ds,'.png'),width = 900, height = 480, units = "px")


EnhancedVolcano(res1,
                     lab = res1$Gene.symbol,
                     x = 'logFC',
                     y = 'adj.P.Val',
                     selectLab = c('CXCL8','CEMIP','MMP7','CA4','ADH1C','GUCA2A',
                                   'GUCA2B','ZG16','CLCA4','MS4A12','CLDN1'),
                     xlab = bquote(~Log[2]~ 'fold change'),
                     pCutoff = .05,
                     FCcutoff = 1.0,
                     pointSize = 2.0,
                     labSize = 3.0,
                     labCol = 'black',
                     labFace = 'plain',
                     boxedLabels = TRUE,
                     colAlpha = 4/5,
                     legendPosition = 'right',
                     legendLabSize = 14,
                     legendIconSize = 4.0,
                     drawConnectors = TRUE,
                     widthConnectors = 1.0,
                     colConnectors = 'black',
                     title = ds)

dev.off() 


# RNA-Seq -----------------------------------------------------------------


res2 <- fread("/home/claire/Documents/ULM/L3/open_science/projet/metastase/GSE50760_cancer.top.table.tsv")

geneF<-fread('/home/claire/Documents/ULM/L3/open_science/projet/metastase/liver_genes.txt' ) # les genes suprexpirmmer dans le foie
genesFS <- geneF[,Gene_Symbol]
res2<-res2[!(Symbol %in% genesFS),]

png(paste0('/home/claire/Documents/ULM/L3/open_science/projet/metastase/genes_controle_cancer.png'),
    width = 900, height = 480, units = "px")

EnhancedVolcano(res2,
                lab = res2$Symbol,
                x = 'log2FoldChange',
                y = "padj",
                selectLab = c('CXCL8','CEMIP','MMP7','CA4','ADH1C','GUCA2A',
                              'GUCA2B','ZG16','CLCA4','MS4A12','CLDN1'),
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'plain',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
                title = "cancer vs controle")

dev.off() 

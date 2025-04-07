if (!require(devtools)) install.packages("devtools")
devtools::install_github("yanlinlin82/ggvenn")

library(ggvenn)
ggvenn(
  file1, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
install.packages("VennDiagram")

library(VennDiagram)
venn.diagram(x, filename = "venn-4-dimensions.png")

file1<- read.csv('/Users/therese/Documents/ENS/L3/S2/Open science/genes_controle_primaire.tsv')
file2<-read.csv('/Users/therese/Documents/ENS/L3/S2/Open science/GSE50760.top.table.csv')
x<-read.csv('/Users/therese/Documents/ENS/L3/S2/Open science/symbole_gene.csv')
x3<-list(x[1],x[2])


ggvenn(
  x3, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 2
)



# Extraire les deux premières colonnes (listes de gènes)
genes3 <- file3bis[,Symbol]
genes2 <- file2bis[,Symbols]
genes1 <- file1bis[,Symbol]
genes4 <- file4bis[,Symbol]
genes5 <- file4bis[,Symbol]
# Créer une liste avec les deux jeux de gènes
gene_lists <- list(" "= genes1, " " = genes2, " "=genes3)

# Créer le diagramme de Venn
venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,  # Pour afficher dans la fenêtre graphique
  fill = c("skyblue", "pink","yellow"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  #cat.pos = c(-20, 20),
  cat.dist = 0.05
)


# Afficher le diagramme
grid.newpage()
grid.draw(venn.plot)

######
install.packages('data.table')
library(data.table)
file2bis<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/GSE50760.top.table.csv' )
file2bis<-file2bis[abs(logFC)>1&P.Value<0.01,]
file1bis<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/genes_controle_primaire.tsv' )
file1bis<-file1bis[abs(log2FC)>1&Pvalue>=2,]
file3bis<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/gene_AllCancer_controle.tsv' )
file3bis<-file3bis[abs(log2FC)>1&Pvalue>=2,]
file4bis<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/genes_controle_metastase.tsv' )
file4bis<-file4bis[abs(log2FC)>1&Pvalue>=2,]
file5bis<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/genes_primaire_metastase.tsv' )
file5bis<-file5bis[abs(log2FC)>1&Pvalue>=2,]

# en filtrant sur pvalue=0.05 et abs(logFC)>1 on a plus que 2748 gene, avec pvalue à 0.01 2391







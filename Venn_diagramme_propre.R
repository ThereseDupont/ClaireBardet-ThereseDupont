### Ce code permet de faire un Venn Diagrammme entre 3 base de donnée qui contiennent le nom de gènes.


## Import des librairies
library(VennDiagram)
library(data.table)

## lire les fichiers

file1<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/genes_controle_primaire.tsv' )
file1bis<-file1[abs(log2FC)>1&Pvalue>=2,]
file2<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/GSE50760.top.table.csv' )
file2bis<-file2[abs(logFC)>1&P.Value<0.01,]
file3<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/gene_AllCancer_controle.tsv' )
file3bis<-file3[abs(log2FC)>1&Pvalue>=2,]

# Supprimer les gènes spécifique du foie

geneF<-fread('/Users/therese/Documents/ENS/L3/S2/Open science/liver_genes_TIGER.txt' ) # les genes suprexpirmmer dans le foie
genesFS <- geneF[,Gene_Symbol]
file3ter<-file3bis[!(Symbol %in% genesFS),]

# Extraire la bonne collones

genes1 <- file1bis[,Symbol]
genes2 <- file2bis[,Symbols]
genes3 <- file3ter[,Symbol]


# Créer une liste avec les différents jeux de gènes
gene_lists <- list("DESeq controle vs primaire "= genes1, " LIMMA controle vs primaire" = genes2, "DESeq cancer(primaire+métastase)vs controle "=genes3)


# Créer le diagramme de Venn
venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,  # Pour afficher dans la fenêtre graphique
  fill = c("skyblue", "pink","yellow"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  #cat.pos = c(-20, 20),
  cat.dist = 0.05)

# Afficher le diagramme
grid.newpage()
grid.draw(venn.plot)

# Fonction pour comparer deux fichiers
compare_files <- function(file1, file2) {
  # Lire les fichiers
  data1 <- read.csv('/Users/therese/Documents/ENS/L3/S2/Open science/top_table_GSE9348_OK', header = TRUE, stringsAsFactors = FALSE)
  data2 <- read.csv('/Users/therese/Documents/ENS/L3/S2/Open science/GSE9348_top_table.csv', header = TRUE, stringsAsFactors = FALSE)
  data1bis <- data1[1,54674]
 
  # Comparer les données
  differences <- which(rowSums(data1bis != data2) > 0)
  
  if (length(differences) == 0) {
    cat("✅ Les fichiers sont identiques !\n")
    return(TRUE)
  } else {
    cat("❌ Les fichiers sont différents !\n")
    cat("Nombre de lignes différentes :", length(differences), "\n")
    
    # Afficher les premières lignes différentes
    cat("\nDétails des différences :\n")
    for (index in head(differences, 5)) {
      cat("\nLigne", index, ":\n")
      cat("Fichier 1 :", paste(as.character(data1[index,]), collapse = " | "), "\n")
      cat("Fichier 2 :", paste(as.character(data2[index,]), collapse = " | "), "\n")
    }
    
    return(FALSE)
  }
}

# Exemple d'utilisation
compare_files( '/Users/therese/Downloads/top_table_GSE9348_OK' , '/Users/therese/Downloads/GSE9348_top_table.csv' )


# ---------------------------------
# 03_similarity_matrix.R
# Build similarity and distance matrices
# ---------------------------------

cat("Reading parsed BLAST file...\n")

blast <- read.table(
  "results/blast/parsed_blast.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# obtener lista de proteínas
proteins <- sort(unique(c(blast$query, blast$subject)))
n <- length(proteins)

cat("Number of proteins:", n, "\n")

# crear matriz de bitscores
bit_matrix <- matrix(
  0,
  nrow = n,
  ncol = n,
  dimnames = list(proteins, proteins)
)

# llenar matriz
for (i in 1:nrow(blast)) {
  q <- blast$query[i]
  s <- blast$subject[i]
  b <- blast$bitscore[i]

  bit_matrix[q, s] <- b
  bit_matrix[s, q] <- b
}

# máximo bitscore fuera de la diagonal
max_score <- max(bit_matrix[bit_matrix != 0])

cat("Max bitscore:", max_score, "\n")

# normalizar matriz
B <- bit_matrix / max_score

# asegurar diagonal = 1
diag(B) <- 1

# convertir a distancia
D <- 1 - B

# convertir a objeto dist
dist_matrix <- as.dist(D)

# crear carpeta si no existe
dir.create("results/matrix", showWarnings = FALSE, recursive = TRUE)

# guardar resultados
write.table(
  B,
  "results/matrix/similarity_matrix.tsv",
  sep = "\t",
  quote = FALSE
)

write.table(
  D,
  "results/matrix/distance_matrix.tsv",
  sep = "\t",
  quote = FALSE
)

saveRDS(
  dist_matrix,
  "results/matrix/distance_matrix.rds"
)

cat("Matrices saved in results/matrix/\n")

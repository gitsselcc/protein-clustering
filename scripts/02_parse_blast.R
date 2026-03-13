# -------------------------------
# 02_parse_blast.R
# Clean BLAST output (outfmt 7)
# -------------------------------

# rutas
input_file <- "results/blast/psbA_blast.tsv"
output_file <- "results/blast/parsed_blast.tsv"

cat("Reading BLAST output...\n")

# leer archivo completo
lines <- readLines(input_file)

# remover líneas de comentarios (#)
blast_lines <- lines[!grepl("^#", lines)]

# convertir a tabla
blast_table <- read.table(
  text = blast_lines,
  sep = "\t",
  stringsAsFactors = FALSE
)

# asignar nombres de columnas (formato BLAST estándar)
colnames(blast_table) <- c(
  "query",
  "subject",
  "identity",
  "alignment_length",
  "mismatch",
  "gapopen",
  "qstart",
  "qend",
  "sstart",
  "send",
  "evalue",
  "bitscore"
)

cat("Rows parsed:", nrow(blast_table), "\n")

# guardar tabla limpia
write.table(
  blast_table,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Parsed BLAST saved to:", output_file, "\n")

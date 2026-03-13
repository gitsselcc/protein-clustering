#!/bin/bash

set -e

# rutas
FASTA="data/raw/psbA_sequences.fasta"
OUTDIR="results/blast"
DB="$OUTDIR/psbA_db"

mkdir -p $OUTDIR

echo "Creating BLAST database..."

makeblastdb \
-in $FASTA \
-dbtype prot \
-out $DB

echo "Running BLAST all-vs-all..."

blastp \
-query $FASTA \
-db $DB \
-out $OUTDIR/psbA_blast.tsv \
-outfmt 6 \
-num_threads 4

echo "BLAST finished."
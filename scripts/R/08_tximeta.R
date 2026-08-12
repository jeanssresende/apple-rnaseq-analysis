# ==========================================================
# Apple RNA-seq Analysis Pipeline
# Module 08 - tximeta import with linked transcriptome
# Author: Jean Resende
# ==========================================================

suppressPackageStartupMessages({
  library(tximeta)
  library(jsonlite)
  library(SummarizedExperiment)
  library(rtracklayer)
})

# ----------------------------------------------------------
# Reference paths
# ----------------------------------------------------------

reference_root <- "../reference/Ensembl_Malus_domestica"

index_dir <- file.path(reference_root, "salmon_index")

fasta_file <- file.path(
  reference_root,
  "transcriptome",
  "Malus_domestica_golden.ASM211411v1.cdna.all.fa.gz"
)

gff3_file <- file.path(
  reference_root,
  "annotation",
  "Malus_domestica_golden.ASM211411v1.63.gff3.gz"
)

json_file <- file.path(reference_root, "linkedTxome.json")

cat("index_dir = ", index_dir, "\n")

# ----------------------------------------------------------
# Quantification directory
# ----------------------------------------------------------

quant_dir <- "../results/salmon_quant"

# ----------------------------------------------------------
# Output directory
# ----------------------------------------------------------

output_dir <- "results/tximeta"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------------------------------------
# Check required reference files
# ----------------------------------------------------------

required_files <- c(fasta_file, gff3_file)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing reference files:\n",
    paste(missing_files, collapse = "\n")
  )
}

if (!file.exists(file.path(index_dir, "info.json"))) {
  stop("Salmon index not found in: ", index_dir)
}

# ----------------------------------------------------------
# Register transcriptome if needed
# ----------------------------------------------------------

if (!file.exists(json_file)) {
  
  cat("====================================================\n")
  cat("Linked transcriptome not found.\n")
  cat("Creating linkedTxome...\n")
  cat("====================================================\n\n")
  
  info_file <- file.path(index_dir, "info.json")
  
  digest_info <- read_json(info_file)
  
  digest <- digest_info$seq_hash
  
  if (is.null(digest) || digest == "") {
    stop("Could not extract seq_hash from Salmon info.json")
  }
  
  cat("Digest: ", substr(digest, 1, 20), "...\n\n")
  
  makeLinkedTxome(
    digest = digest,
    indexName = "Malus_domestica_ASM211411v1_EnsemblPlants63",
    source = "Ensembl Plants",
    organism = "Malus domestica",
    release = "63",
    genome = "ASM211411v1",
    fasta = fasta_file,
    gtf = gff3_file,
    jsonFile = json_file
  )
  
  cat("\nlinkedTxome created successfully:\n")
  cat(json_file, "\n\n")
  
} else {
  
  cat("====================================================\n")
  cat("Loading existing linkedTxome...\n")
  cat("====================================================\n\n")
  
  loadLinkedTxome(json_file)
}

# ----------------------------------------------------------
# Locate quant.sf files
# ----------------------------------------------------------

files <- list.files(
  quant_dir,
  pattern = "quant.sf$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No quant.sf files found in: ", quant_dir)
}

sample_names <- basename(dirname(files))

# ----------------------------------------------------------
# Build coldata
# ----------------------------------------------------------

coldata <- data.frame(
  files = files,
  names = sample_names,
  stringsAsFactors = FALSE
)

cat("Samples detected: ", nrow(coldata), "\n\n")

# ----------------------------------------------------------
# Import with tximeta
# ----------------------------------------------------------

cat("Importing Salmon quantifications with tximeta...\n\n")

se <- tximeta(coldata)

# ----------------------------------------------------------
# Summarize transcripts to genes
# ----------------------------------------------------------

gse <- summarizeToGene(se)

# ----------------------------------------------------------
# Import annotation from GFF3
# ----------------------------------------------------------

cat("Importing GFF3 annotation...\n\n")

gff <- import(gff3_file)

cat("Feature types in GFF3:\n")
print(table(gff$type))

# Keep only genes and ncRNA genes
genes_gff <- gff[gff$type %in% c("gene", "ncRNA_gene")]

# Build annotation table
annotation_df <- data.frame(
  gene_id       = genes_gff$gene_id,
  gene_name     = genes_gff$Name,
  external_name = genes_gff$external_name,
  biotype       = genes_gff$biotype,
  description   = genes_gff$description,
  stringsAsFactors = FALSE
)

# Remove duplicated gene IDs
annotation_df <- annotation_df[
  !duplicated(annotation_df$gene_id),
]

# Match annotation to gse
idx <- match(rowData(gse)$gene_id, annotation_df$gene_id)

# Add annotation to gse
rowData(gse)$gene_name     <- annotation_df$gene_name[idx]
rowData(gse)$external_name <- annotation_df$external_name[idx]
rowData(gse)$biotype       <- annotation_df$biotype[idx]
rowData(gse)$description   <- annotation_df$description[idx]

cat("Annotation added to gse.\n\n")

# ----------------------------------------------------------
# Create protein-coding subset
# ----------------------------------------------------------

gse_pc <- gse[rowData(gse)$biotype == "protein_coding", ]

cat("Protein-coding genes: ", nrow(gse_pc), "\n\n")

# ----------------------------------------------------------
# Save objects
# ----------------------------------------------------------

se_file  <- file.path(output_dir, "se_transcript.rds")
gse_file <- file.path(output_dir, "gse_gene_annotated.rds")
pc_file  <- file.path(output_dir, "gse_protein_coding.rds")

saveRDS(se, se_file)
saveRDS(gse, gse_file)
saveRDS(gse_pc, pc_file)

# ----------------------------------------------------------
# Summary
# ----------------------------------------------------------

cat("\n====================================================\n")
cat("tximeta completed successfully!\n")
cat("====================================================\n")
cat("Transcript-level features : ", nrow(se), "\n")
cat("Gene-level features       : ", nrow(gse), "\n")
cat("Protein-coding genes      : ", nrow(gse_pc), "\n")
cat("Samples                   : ", ncol(gse), "\n")
cat("Output directory          : ", output_dir, "\n")
cat("LinkedTxome               : ", json_file, "\n")
cat("====================================================\n\n")

# ----------------------------------------------------------
# Preview
# ----------------------------------------------------------

cat("Preview of annotated genes:\n\n")

print(
  rowData(gse)[1:5,
               c("gene_id",
                 "gene_name",
                 "external_name",
                 "biotype")]
)

# ----------------------------------------------------------
# Validate output files
# ----------------------------------------------------------

cat("\nFiles created:\n")
cat("  se_transcript.rds       : ", file.exists(se_file), "\n")
cat("  gse_gene_annotated.rds : ", file.exists(gse_file), "\n")
cat("  gse_protein_coding.rds : ", file.exists(pc_file), "\n")
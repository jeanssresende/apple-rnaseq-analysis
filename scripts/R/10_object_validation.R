# ============================================================================
# MÓDULO 10: VALIDAÇÃO DOS OBJETOS IMPORTADOS
# ============================================================================
# Objetivo: Inspecionar dimensões, distribuição de contagens, valores ausentes
#           e estatísticas descritivas antes da análise exploratória.
#
# Entrada: gse_with_metadata.rds (SummarizedExperiment com metadata)
# Saídas:
#   - Relatório de validação (console + arquivo de log)
#   - Gráficos de QC
#   - Arquivo de estatísticas descritivas
# ============================================================================

# Limpar ambiente
rm(list = ls())

# ============================================================================
# 1. CARREGAR PACOTES
# ============================================================================

cat("\n=== MÓDULO 10: VALIDAÇÃO DOS OBJETOS IMPORTADOS ===\n")

library(SummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(reshape2)

# ============================================================================
# 2. DEFINIR CAMINHOS
# ============================================================================

base_dir <- here::here()
results_dir <- file.path(base_dir, "results", "tximeta")
qc_dir <- file.path(base_dir, "results", "qc")
log_dir <- file.path(base_dir, "logs")

# Criar diretórios se não existirem
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

# Arquivo de log
log_file <- file.path(log_dir, "module_10_validation.log")
sink(log_file, append = TRUE)

cat("\n--- Iniciando Módulo 10 ---\n")
cat("Data/Hora:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ============================================================================
# 3. CARREGAR OBJETO
# ============================================================================

cat("\n--- Carregando objeto gse ---\n")

gse_file <- file.path(results_dir, "gse_with_metadata.rds")

if (!file.exists(gse_file)) {
  cat("ERRO: Arquivo não encontrado:", gse_file, "\n")
  stop("Execute Módulo 09 primeiro.")
}

gse <- readRDS(gse_file)
cat("✓ Objeto carregado com sucesso\n")

# ============================================================================
# 4. VALIDAÇÃO BÁSICA DE DIMENSÕES
# ============================================================================

cat("\n--- DIMENSÕES DO OBJETO ---\n")

n_genes <- nrow(gse)
n_samples <- ncol(gse)
n_assays <- length(assays(gse))

cat("Genes:", n_genes, "\n")
cat("Amostras:", n_samples, "\n")
cat("Assays disponíveis:", n_assays, "\n")
cat("Nomes dos assays:", paste(names(assays(gse)), collapse = ", "), "\n")

# Validar que temos pelo menos 2 assays (abundance e counts)
if (n_assays < 2) {
  cat("Aviso: Esperado pelo menos 2 assays (abundance, counts)\n")
}

# ============================================================================
# 5. VALIDAÇÃO DE METADATA
# ============================================================================

cat("\n--- METADATA (colData) ---\n")

metadata <- colData(gse)
cat("Número de amostras:", nrow(metadata), "\n")
cat("Colunas de metadata:\n")
print(colnames(metadata))

cat("\nEstrutura de metadata:\n")
print(str(metadata))

cat("\nResumo de metadata:\n")
print(summary(metadata))

# ============================================================================
# 6. VALIDAÇÃO DE ANOTAÇÃO GÊNICA
# ============================================================================

cat("\n--- ANOTAÇÃO GÊNICA (rowData) ---\n")

gene_data <- rowData(gse)
cat("Colunas de anotação:\n")
print(colnames(gene_data))

# Verificar gene_id
if ("gene_id" %in% colnames(gene_data)) {
  cat("\n✓ gene_id presente\n")
  cat("  Únicos:", n_distinct(gene_data$gene_id), "\n")
  cat("  NA:", sum(is.na(gene_data$gene_id)), "\n")
} else {
  cat("gene_id NÃO encontrado\n")
}

# Verificar gene_name
if ("gene_name" %in% colnames(gene_data)) {
  n_with_name <- sum(!is.na(gene_data$gene_name))
  pct_with_name <- (n_with_name / n_genes) * 100
  cat("\n✓ gene_name presente\n")
  cat("  Com gene_name:", n_with_name, sprintf("(%.2f%%)", pct_with_name), "\n")
  cat("  Sem gene_name (NA):", sum(is.na(gene_data$gene_name)), "\n")
} else {
  cat("gene_name NÃO encontrado\n")
}

# Verificar biotype
if ("biotype" %in% colnames(gene_data)) {
  cat("\n✓ biotype presente\n")
  cat("Distribuição de biotypes:\n")
  print(table(gene_data$biotype))
} else {
  cat("biotype NÃO encontrado\n")
}

# ============================================================================
# 7. VALIDAÇÃO DE CONTAGENS
# ============================================================================

cat("\n--- DISTRIBUIÇÃO DE CONTAGENS ---\n")

# Extrair counts
counts_matrix <- assay(gse, "counts")

cat("Dimensões da matriz de contagens:", dim(counts_matrix), "\n")

# Verificar valores ausentes
n_na <- sum(is.na(counts_matrix))
cat("Valores NA na matriz:", n_na, "\n")

if (n_na > 0) {
  cat("Aviso: Encontrados", n_na, "valores NA\n")
}

# Verificar valores negativos (não deve haver)
n_negative <- sum(counts_matrix < 0, na.rm = TRUE)
cat("Valores negativos:", n_negative, "\n")

if (n_negative > 0) {
  cat("ERRO: Encontrados valores negativos em contagens\n")
}

# Verificar zeros
n_zeros <- sum(counts_matrix == 0, na.rm = TRUE)
pct_zeros <- (n_zeros / (n_genes * n_samples)) * 100
cat("Zeros na matriz:", n_zeros, sprintf("(%.2f%%)", pct_zeros), "\n")

# ============================================================================
# 8. ESTATÍSTICAS POR AMOSTRA
# ============================================================================

cat("\n--- ESTATÍSTICAS POR AMOSTRA ---\n")

sample_stats <- tibble(
  sample = colnames(gse),
  total_counts = colSums(counts_matrix),
  mean_counts = colMeans(counts_matrix),
  median_counts = apply(counts_matrix, 2, median),
  sd_counts = apply(counts_matrix, 2, sd),
  min_counts = apply(counts_matrix, 2, min),
  max_counts = apply(counts_matrix, 2, max),
  n_expressed = colSums(counts_matrix > 0),
  pct_expressed = (colSums(counts_matrix > 0) / n_genes) * 100
)

cat("\nResumo de estatísticas por amostra:\n")
print(summary(sample_stats[, -1]))

cat("\nTabela completa:\n")
print(sample_stats)

# Salvar estatísticas
stats_file <- file.path(qc_dir, "sample_statistics.csv")
write_csv(sample_stats, stats_file)
cat("\n✓ Estatísticas salvas em:", stats_file, "\n")

# ============================================================================
# 9. ESTATÍSTICAS POR GENE
# ============================================================================

cat("\n--- ESTATÍSTICAS POR GENE ---\n")

gene_stats <- tibble(
  gene_id = rownames(gse),
  mean_counts = rowMeans(counts_matrix),
  median_counts = apply(counts_matrix, 1, median),
  sd_counts = apply(counts_matrix, 1, sd),
  min_counts = apply(counts_matrix, 1, min),
  max_counts = apply(counts_matrix, 1, max),
  n_samples_expressed = rowSums(counts_matrix > 0),
  pct_samples_expressed = (rowSums(counts_matrix > 0) / n_samples) * 100
)

cat("\nResumo de estatísticas por gene:\n")
print(summary(gene_stats[, -1]))

# Genes não expressos (zero em todas as amostras)
n_not_expressed <- sum(gene_stats$n_samples_expressed == 0)
cat("\nGenes não expressos (zero em todas amostras):", n_not_expressed, "\n")

# Genes expressos em pelo menos uma amostra
n_expressed_any <- sum(gene_stats$n_samples_expressed > 0)
cat("Genes expressos (>0 em ≥1 amostra):", n_expressed_any, "\n")

# Genes expressos em todas as amostras
n_expressed_all <- sum(gene_stats$n_samples_expressed == n_samples)
cat("Genes expressos em TODAS as amostras:", n_expressed_all, "\n")

# Salvar estatísticas
gene_stats_file <- file.path(qc_dir, "gene_statistics.csv")
write_csv(gene_stats, gene_stats_file)
cat("\n✓ Estatísticas salvas em:", gene_stats_file, "\n")

# ============================================================================
# 10. IDENTIFICAÇÃO DE POSSÍVEIS OUTLIERS
# ============================================================================

cat("\n--- IDENTIFICAÇÃO DE POSSÍVEIS OUTLIERS ---\n")

# Amostras com total de contagens muito baixo/alto
q1_total <- quantile(sample_stats$total_counts, 0.25)
q3_total <- quantile(sample_stats$total_counts, 0.75)
iqr_total <- q3_total - q1_total
lower_bound <- q1_total - 1.5 * iqr_total
upper_bound <- q3_total + 1.5 * iqr_total

outliers_total <- sample_stats %>%
  filter(total_counts < lower_bound | total_counts > upper_bound)

if (nrow(outliers_total) > 0) {
  cat("Possíveis outliers por total de contagens:\n")
  print(outliers_total)
} else {
  cat("✓ Nenhum outlier detectado por total de contagens\n")
}

# ============================================================================
# 11. CORRELAÇÃO ENTRE REPLICATAS
# ============================================================================

cat("\n--- CORRELAÇÃO ENTRE REPLICATAS ---\n")

# Extrair informações de replicata e condição
if ("status" %in% colnames(metadata)) {
  cat("status disponível em metadata\n")
  
  # Agrupar por status
  status_groups <- unique(metadata$status)
  cat("Grupos de status:", paste(status_groups, collapse = ", "), "\n")
  
  for (status in status_groups) {
    samples_status <- which(metadata$status == status)
    if (length(samples_status) > 1) {
      # Calcular correlação entre replicatas
      counts_subset <- counts_matrix[, samples_status]
      # Log2 transform para melhor visualização (adicionar pseudocount)
      counts_log2 <- log2(counts_subset + 1)
      cor_matrix <- cor(counts_log2)
      
      cat("\n  status:", status, "\n")
      cat("  Número de replicatas:", length(samples_status), "\n")
      cat("  Correlação média entre replicatas:\n")
      
      # Extrair correlações fora da diagonal
      diag(cor_matrix) <- NA
      mean_cor <- mean(cor_matrix, na.rm = TRUE)
      min_cor <- min(cor_matrix, na.rm = TRUE)
      max_cor <- max(cor_matrix, na.rm = TRUE)
      
      cat("    Média:", sprintf("%.4f", mean_cor), "\n")
      cat("    Mín:", sprintf("%.4f", min_cor), "\n")
      cat("    Máx:", sprintf("%.4f", max_cor), "\n")
    }
  }
} else {
  cat("status não encontrado em metadata\n")
}

# ============================================================================
# 12. GRÁFICOS DE QC
# ============================================================================

cat("\n--- GERANDO GRÁFICOS DE QC ---\n")

# Gráfico 1: Distribuição de total de contagens por amostra
p1 <- ggplot(sample_stats, aes(x = reorder(sample, total_counts), y = total_counts)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Total de contagens por amostra",
    x = "Amostra",
    y = "Total de contagens"
  )

ggsave(file.path(qc_dir, "01_total_counts_per_sample.png"), p1, width = 12, height = 6)
cat("✓ Salvo: 01_total_counts_per_sample.png\n")

# Gráfico 2: Distribuição de genes expressos por amostra
p2 <- ggplot(sample_stats, aes(x = reorder(sample, n_expressed), y = n_expressed)) +
  geom_col(fill = "darkgreen", alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Número de genes expressos (>0) por amostra",
    x = "Amostra",
    y = "Genes expressos"
  )

ggsave(file.path(qc_dir, "02_expressed_genes_per_sample.png"), p2, width = 12, height = 6)
cat("✓ Salvo: 02_expressed_genes_per_sample.png\n")

# Gráfico 3: Boxplot de contagens por amostra (log scale)
counts_log2 <- log2(counts_matrix + 1)
counts_long <- melt(counts_log2)
colnames(counts_long) <- c("gene", "sample", "log2_counts")

p3 <- ggplot(counts_long, aes(x = sample, y = log2_counts)) +
  geom_boxplot(fill = "lightblue", alpha = 0.7) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Distribuição de contagens por amostra (log2 scale)",
    x = "Amostra",
    y = "log2(counts + 1)"
  )

ggsave(file.path(qc_dir, "03_counts_distribution_boxplot.png"), p3, width = 12, height = 6)
cat("✓ Salvo: 03_counts_distribution_boxplot.png\n")

# Gráfico 4: Distribuição de contagens por gene (histograma)
p4 <- ggplot(gene_stats, aes(x = mean_counts)) +
  geom_histogram(bins = 50, fill = "coral", alpha = 0.7) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "Distribuição de expressão média por gene",
    x = "Mean counts (log10 scale)",
    y = "Frequência"
  )

ggsave(file.path(qc_dir, "04_gene_expression_distribution.png"), p4, width = 10, height = 6)
cat("✓ Salvo: 04_gene_expression_distribution.png\n")

# Gráfico 5: Genes vs amostras expressos
p5 <- ggplot(gene_stats, aes(x = pct_samples_expressed)) +
  geom_histogram(bins = 50, fill = "purple", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Distribuição de genes por porcentagem de amostras com expressão",
    x = "% de amostras com expressão",
    y = "Número de genes"
  )

ggsave(file.path(qc_dir, "05_gene_detection_across_samples.png"), p5, width = 10, height = 6)
cat("✓ Salvo: 05_gene_detection_across_samples.png\n")

# ============================================================================
# 13. RESUMO FINAL
# ============================================================================

cat("\n=== RESUMO FINAL DE VALIDAÇÃO ===\n")

cat("\n✓ VALIDAÇÕES CONCLUÍDAS COM SUCESSO\n")
cat("\nArquivos gerados:\n")
cat("  - results/qc/sample_statistics.csv\n")
cat("  - results/qc/gene_statistics.csv\n")
cat("  - results/qc/01_total_counts_per_sample.png\n")
cat("  - results/qc/02_expressed_genes_per_sample.png\n")
cat("  - results/qc/03_counts_distribution_boxplot.png\n")
cat("  - results/qc/04_gene_expression_distribution.png\n")
cat("  - results/qc/05_gene_detection_across_samples.png\n")

cat("\n--- Próximo passo: Módulo 11 (PCA e Análise Exploratória) ---\n")

# Fechar log
sink()

cat("\n✓ Log salvo em:", log_file, "\n")

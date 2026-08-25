# ============================================================================
# MÓDULO 11: ANÁLISE EXPLORATÓRIA - PCA E CLUSTERING
# ============================================================================
# Objetivo: Explorar estrutura dos dados, identificar padrões, outliers
#           e validar desenho experimental antes de DESeq2.
#
# Entrada: gse_with_metadata.rds
# Saídas:
#   - PCA plots (colorido por Status, Variety, etc)
#   - Heatmap de correlação entre amostras
#   - Gráficos de distância entre amostras
#   - Identificação de outliers
# ============================================================================

rm(list = ls())

# ============================================================================
# 1. CARREGAR PACOTES
# ============================================================================

cat("\n=== MÓDULO 11: ANÁLISE EXPLORATÓRIA - PCA ===\n")

library(SummarizedExperiment)
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(here)

# ============================================================================
# 2. DEFINIR CAMINHOS
# ============================================================================

base_dir <- here::here()
results_dir <- file.path(base_dir, "results", "tximeta")
pca_dir <- file.path(base_dir, "results", "pca_eda")
log_dir <- file.path(base_dir, "logs")

dir.create(pca_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(log_dir, "module_11_pca_eda.log")
sink(log_file, append = TRUE)

cat("\n--- Iniciando Módulo 11 ---\n")
cat("Data/Hora:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# ============================================================================
# 3. CONSTRUIR DESeqDataSet
# ============================================================================

cat("\n--- Construindo DESeqDataSet ---\n")

# Verificar colData
cat("Colunas do colData:\n")
print(colnames(colData(gse)))

# Construir DDS com design correto
dds <- DESeqDataSet(se = gse, design = ~ status)

cat("✓ DESeqDataSet construído\n")
cat("Dimensões do dds:\n")
print(dim(dds))

# ============================================================================
# 4. CRIAR OBJETO DESeqDataSet
# ============================================================================

cat("\n--- Criando DESeqDataSet ---\n")

# Extrair design experimental
# Vamos usar ~ Status como design principal
dds <- DESeqDataSet(gse, design = ~ status)

cat("✓ DESeqDataSet criado\n")
cat("Design:", as.character(design(dds)), "\n")

# ============================================================================
# 5. PRÉ-FILTRAGEM: REMOVER GENES COM EXPRESSÃO MUITO BAIXA
# ============================================================================

cat("\n--- Pré-filtragem de genes ---\n")

# Remover genes com < 10 contagens em < 3 amostras
# Justificativa: Genes com muito poucos reads são ruído
keep <- rowSums(counts(dds) >= 10) >= 3

cat("Genes antes da filtragem:", nrow(dds), "\n")
cat("Genes após filtragem:", sum(keep), "\n")
cat("Genes removidos:", nrow(dds) - sum(keep), "\n")

dds <- dds[keep, ]

cat("✓ Filtragem concluída\n")

# ============================================================================
# 6. NORMALIZAÇÃO COM rlog (REGULARIZED LOG)
# ============================================================================

cat("\n--- Normalizando com rlog ---\n")

# rlog é mais apropriado para PCA/clustering que vst
# quando n_samples < 30 (nosso caso)
rld <- rlog(dds, blind = TRUE)

cat("✓ rlog concluído\n")
cat("Dimensões após normalização:", dim(rld), "\n")

# ============================================================================
# 7. PCA - MANUAL (PARA MAIS CONTROLE)
# ============================================================================

cat("\n--- Executando PCA ---\n")

# Extrair dados normalizados
normalized_counts <- assay(rld)

# Centralizar e escalar
data_centered <- scale(t(normalized_counts), center = TRUE, scale = TRUE)

# PCA
pca_result <- prcomp(data_centered)

# Variância explicada
variance_explained <- (pca_result$sdev^2 / sum(pca_result$sdev^2)) * 100

cat("Variância explicada pelos primeiros 5 PCs:\n")
print(variance_explained[1:5])

# ============================================================================
# 8. PREPARAR DADOS PARA VISUALIZAÇÃO
# ============================================================================

cat("\n--- Preparando dados para visualização ---\n")

# Extrair scores (PCA coordinates)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$sample <- rownames(pca_scores)

# Extrair metadata diretamente de colData(rld)
# Usar rownames para garantir correspondência correta
metadata <- as.data.frame(colData(rld))
metadata$sample <- rownames(metadata)

# Fazer merge pelo nome da amostra (que já é a chave correta)
pca_scores <- pca_scores %>%
  left_join(metadata %>% select(sample, sample_id, status, variety), 
            by = "sample")

# Validar se houve correspondência
if (any(is.na(pca_scores$status))) {
  cat("ERRO: Alguns metadados não foram encontrados\n")
  print(pca_scores %>% filter(is.na(status)))
  stop("Falha no merge de metadata")
}

cat("✓ Dados preparados para visualização\n")
cat("Primeiras 5 linhas:\n")
print(head(pca_scores))

# ============================================================================
# 9. GRÁFICO 1: PCA - PC1 vs PC2 (colorido por Status)
# ============================================================================

cat("\n--- Gerando PCA plots ---\n")

p1 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = status, shape = variety)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 10) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "PCA: PC1 vs PC2",
    x = sprintf("PC1 (%.2f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.2f%%)", variance_explained[2]),
    color = "status",
    shape = "variety"
  )

ggsave(file.path(pca_dir, "01_pca_pc1_pc2.png"), p1, width = 12, height = 8)
cat("✓ Salvo: 01_pca_pc1_pc2.png\n")

# ============================================================================
# 10. GRÁFICO 2: PCA - PC1 vs PC3
# ============================================================================

p2 <- ggplot(pca_scores, aes(x = PC1, y = PC3, color = status, shape = variety)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(aes(label = sample), size = 3, max.overlaps = 10) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(
    title = "PCA: PC1 vs PC3",
    x = sprintf("PC1 (%.2f%%)", variance_explained[1]),
    y = sprintf("PC3 (%.2f%%)", variance_explained[3]),
    color = "status",
    shape = "variety"
  )

ggsave(file.path(pca_dir, "02_pca_pc1_pc3.png"), p2, width = 12, height = 8)
cat("✓ Salvo: 02_pca_pc1_pc3.png\n")

# ============================================================================
# 11. GRÁFICO 3: Variância explicada (Scree plot)
# ============================================================================

scree_data <- tibble(
  PC = 1:length(variance_explained),
  Variance = variance_explained,
  Cumulative = cumsum(variance_explained)
)

p3 <- ggplot(scree_data[1:10, ], aes(x = factor(PC), y = Variance)) +
  geom_col(fill = "steelblue", alpha = 0.7) +
  geom_line(aes(x = PC, y = Cumulative), color = "red", size = 1, group = 1) +
  geom_point(aes(x = PC, y = Cumulative), color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = "Scree Plot - Variância Explicada",
    x = "Componente Principal",
    y = "Variância Explicada (%)",
    subtitle = "Linha vermelha = variância cumulativa"
  )

ggsave(file.path(pca_dir, "03_scree_plot.png"), p3, width = 10, height = 6)
cat("✓ Salvo: 03_scree_plot.png\n")

# ============================================================================
# 12. GRÁFICO 4: Heatmap de correlação entre amostras
# ============================================================================

cat("\n--- Gerando heatmap de correlação ---\n")

# Extrair metadata de colData e converter para data.frame
metadata <- as.data.frame(colData(gse))

# Calcular matriz de distância Euclidiana
sample_dist <- dist(data_centered)
sample_dist_matrix <- as.matrix(sample_dist)

# Converter para correlação (1 - distância normalizada)
sample_cor <- 1 - (sample_dist_matrix / max(sample_dist_matrix))

# Preparar anotações - extrair apenas status e variety
annotation_col <- metadata[, c("status", "variety"), drop = FALSE]

# Validar correspondência
stopifnot(rownames(annotation_col) == rownames(sample_cor))

# Heatmap
png(file.path(pca_dir, "04_sample_correlation_heatmap.png"), width = 10, height = 8, units = "in", res = 300)
pheatmap(
  sample_cor,
  annotation_col = annotation_col,
  annotation_row = annotation_col,
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Correlação entre amostras (rlog-normalized)",
  color = colorRampPalette(c("white", "yellow", "red"))(50),
  breaks = seq(0, 1, length.out = 51)
)
dev.off()
cat("✓ Salvo: 04_sample_correlation_heatmap.png\n")

# ============================================================================
# 13. DETECÇÃO DE OUTLIERS (BASEADO EM DISTÂNCIA)
# ============================================================================

cat("\n--- Detectando possíveis outliers ---\n")

# Calcular distância média de cada amostra para as outras
mean_dist <- colMeans(sample_dist_matrix)

# Threshold: média + 2*sd
threshold <- mean(mean_dist) + 2 * sd(mean_dist)

outliers <- names(mean_dist[mean_dist > threshold])

if (length(outliers) > 0) {
  cat("Possíveis outliers detectados:\n")
  for (out in outliers) {
    cat("  -", out, "(distância média:", sprintf("%.2f", mean_dist[out]), ")\n")
  }
} else {
  cat("✓ Nenhum outlier detectado\n")
}

# ============================================================================
# 14. GRÁFICO 5: Distância média entre amostras
# ============================================================================

dist_data <- tibble(
  sample = names(mean_dist),
  mean_distance = mean_dist
) %>%
  left_join(metadata, by = c("sample" = "sample_id"))

p5 <- ggplot(dist_data, aes(x = reorder(sample, mean_distance), y = mean_distance, fill = status)) +
  geom_col(alpha = 0.7) +
  geom_hline(yintercept = threshold, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(
    title = "Distância média de cada amostra para as outras",
    x = "Amostra",
    y = "Distância média",
    subtitle = "Linha vermelha = threshold para outliers"
  )

ggsave(file.path(pca_dir, "05_sample_distances.png"), p5, width = 12, height = 6)
cat("✓ Salvo: 05_sample_distances.png\n")

# ============================================================================
# 15. GRÁFICO 6: PCA colorido apenas por Status (simplificado)
# ============================================================================

p6 <- ggplot(pca_scores, aes(x = PC1, y = PC2, color = status)) +
  geom_point(size = 5, alpha = 0.7) +
  stat_ellipse(level = 0.95, linetype = "dashed") +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    text = element_text(size = 12)
  ) +
  labs(
    title = "PCA: Clean vs Infected",
    x = sprintf("PC1 (%.2f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.2f%%)", variance_explained[2]),
    color = "Status"
  ) +
  scale_color_manual(values = c("clean" = "green", "infected" = "red"))

ggsave(file.path(pca_dir, "06_pca_status_ellipse.png"), p6, width = 10, height = 8)
cat("✓ Salvo: 06_pca_status_ellipse.png\n")

# ============================================================================
# 16. SALVAR OBJETOS PARA PRÓXIMOS MÓDULOS
# ============================================================================

cat("\n--- Salvando objetos para próximos módulos ---\n")

# Salvar DDS (será usado em DESeq2)
dds_file <- file.path(results_dir, "dds_filtered.rds")
saveRDS(dds, dds_file)
cat("✓ Salvo: dds_filtered.rds\n")

# Salvar rld (será usado em visualizações posteriores)
rld_file <- file.path(results_dir, "rld_normalized.rds")
saveRDS(rld, rld_file)
cat("✓ Salvo: rld_normalized.rds\n")

# ============================================================================
# 17. RESUMO FINAL
# ============================================================================

cat("\n=== RESUMO DA ANÁLISE EXPLORATÓRIA ===\n")

cat("\nVariância explicada pelos primeiros 5 PCs:\n")
print(data.frame(
  PC = 1:5,
  Variance_pct = variance_explained[1:5],
  Cumulative_pct = cumsum(variance_explained)[1:5]
))

cat("\nSeparação de grupos (Status):\n")
status_groups <- pca_scores %>%
  group_by(status) %>%
  summarise(
    mean_PC1 = mean(PC1),
    mean_PC2 = mean(PC2),
    sd_PC1 = sd(PC1),
    sd_PC2 = sd(PC2),
    .groups = "drop"
  )
print(status_groups)

cat("\n✓ ANÁLISE EXPLORATÓRIA CONCLUÍDA\n")

cat("\nArquivos gerados:\n")
cat("  - results/pca_eda/01_pca_pc1_pc2.png\n")
cat("  - results/pca_eda/02_pca_pc1_pc3.png\n")
cat("  - results/pca_eda/03_scree_plot.png\n")
cat("  - results/pca_eda/04_sample_correlation_heatmap.png\n")
cat("  - results/pca_eda/05_sample_distances.png\n")
cat("  - results/pca_eda/06_pca_status_ellipse.png\n")

cat("\nObjetos salvos:\n")
cat("  - results/tximeta/dds_filtered.rds\n")
cat("  - results/tximeta/rld_normalized.rds\n")

cat("\n--- Próximo passo: Módulo 12 (DESeq2 - Expressão Diferencial) ---\n")

sink()

cat("\n✓ Log salvo em:", log_file, "\n")

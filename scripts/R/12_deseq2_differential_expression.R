# ============================================================================
# MÓDULO 12: ANÁLISE DE EXPRESSÃO DIFERENCIAL COM DESeq2
# ============================================================================

library(DESeq2)
library(tidyverse)

cat("\n=== MÓDULO 12: DESeq2 - Análise de Expressão Diferencial ===\n")

# ============================================================================
# 1. CARREGAR OBJETO
# ============================================================================

cat("\n--- Carregando objeto gene-level com anotação ---\n")

gse <- readRDS("results/tximeta/gse_with_metadata.rds")

cat(sprintf("✓ Objeto carregado: %d genes × %d amostras\n", nrow(gse), ncol(gse)))
print(gse)

# ============================================================================
# 2. CONSTRUIR DESeqDataSet
# ============================================================================

cat("\n--- Construindo DESeqDataSet ---\n")

# Design: ~variety + status
# Isso controla para a variação entre cultivares
dds <- DESeqDataSet(
  se = gse,
  design = ~ variety + status
)

cat("✓ DESeqDataSet construído\n")
print(dds)

# ============================================================================
# 3. CONVERTER PARA FATORES E VALIDAR
# ============================================================================

cat("\n--- Convertendo para fatores ---\n")

# Converter para fatores se ainda não estiverem
if (!is.factor(dds$status)) {
  dds$status <- factor(dds$status)
}
if (!is.factor(dds$variety)) {
  dds$variety <- factor(dds$variety)
}

cat("Status:\n")
print(table(dds$status))

cat("\nVariety:\n")
print(table(dds$variety))

# Definir referência (baseline)
# Clean como referência (baseline)
dds$status <- relevel(dds$status, ref = "clean")

cat("✓ Referência definida: status = 'clean'\n")

# ============================================================================
# 4. EXECUTAR DESeq2
# ============================================================================

cat("\n--- Executando DESeq2 (estimativa de size factors, dispersão e teste) ---\n")

dds <- DESeq(dds, parallel = TRUE)

cat("✓ DESeq2 concluído\n")

# ============================================================================
# 5. EXTRAIR RESULTADOS
# ============================================================================

cat("\n--- Extraindo resultados: Infected vs Clean ---\n")

# Contraste: Infected vs Clean (referência)
res <- results(
  dds,
  contrast = c("status", "infected", "clean"),
  alpha = 0.05
)

cat(sprintf("✓ Resultados extraídos\n"))
cat(sprintf("  - Total de genes testados: %d\n", nrow(res)))
cat(sprintf("  - Genes significativos (padj < 0.05): %d\n", sum(res$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  - Upregulados (log2FC > 0): %d\n", sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)))
cat(sprintf("  - Downregulados (log2FC < 0): %d\n", sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)))

print(head(res))

# ============================================================================
# 6. SALVAR RESULTADOS
# ============================================================================

cat("\n--- Salvando resultados ---\n")

# Criar diretório
dir.create("results/deseq2", showWarnings = FALSE)

# Salvar objetos RDS
saveRDS(dds, "results/deseq2/dds.rds")
saveRDS(res, "results/deseq2/results_infected_vs_clean.rds")

# Converter resultados para data.frame
res_df <- as.data.frame(res)

# Verificar se gene_id já existe como coluna
if ("gene_id" %in% colnames(res_df)) {
  # Se já existe, apenas adicionar rownames como coluna com outro nome
  res_df$gene_id_rowname <- rownames(res_df)
  res_df <- res_df %>%
    select(gene_id_rowname, everything()) %>%
    rename(gene_id = gene_id_rowname)
} else {
  # Se não existe, converter rownames para coluna
  res_df <- res_df %>%
    rownames_to_column("gene_id")
}

# Ordenar por p-value ajustado
res_df <- res_df %>% arrange(padj)

# Salvar CSV
write.csv(res_df, "results/deseq2/results_infected_vs_clean.csv", row.names = FALSE)

cat("✓ Resultados salvos em results/deseq2/\n")
cat("Arquivo: results/deseq2/results_infected_vs_clean.csv\n")

# Visualizar top resultados
cat("\nTop 20 genes diferencialmente expressos:\n")
print(head(res_df, 20))

cat("\n=== FIM DO MÓDULO 12 ===\n")
# ============================================================================
# 7. VISUALIZAÇÕES BÁSICAS
# ============================================================================

cat("\n--- Gerando visualizações ---\n")

# MA plot
png(file.path(de_dir, "01_ma_plot.png"), width = 10, height = 8, units = "in", res = 300)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: Infected vs Clean")
dev.off()
cat("✓ Salvo: 01_ma_plot.png\n")

# Volcano plot
png(file.path(de_dir, "02_volcano_plot.png"), width = 10, height = 8, units = "in", res = 300)
plot(
  res$log2FoldChange, -log10(res$padj),
  main = "Volcano Plot: Infected vs Clean",
  xlab = "log2(Fold Change)",
  ylab = "-log10(padj)",
  pch = 16,
  col = ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "red", "gray")
)
abline(v = c(-1, 1), col = "blue", lty = 2)
abline(h = -log10(0.05), col = "blue", lty = 2)
legend("topright", legend = "padj < 0.05, |log2FC| > 1", col = "red", pch = 16)
dev.off()
cat("✓ Salvo: 02_volcano_plot.png\n")

# Distribuição de log2FC
png(file.path(de_dir, "03_log2fc_distribution.png"), width = 10, height = 8, units = "in", res = 300)
hist(
  res$log2FoldChange,
  breaks = 50,
  main = "Distribuição de log2(Fold Change)",
  xlab = "log2(Fold Change)",
  col = "steelblue"
)
dev.off()
cat("✓ Salvo: 03_log2fc_distribution.png\n")

cat("\n=== MÓDULO 12 CONCLUÍDO ===\n")

# ============================================================================
# 7b. VOLCANO PLOT (versão ggplot2)
# ============================================================================

cat("\n--- Gerando Volcano Plot (ggplot2) ---\n")

res_df_plot <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Sig"),
    direction = ifelse(log2FoldChange > 0, "Up", "Down")
  )

p_volcano <- ggplot(res_df_plot, aes(x = log2FoldChange, y = -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Sig" = "gray")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", alpha = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue", alpha = 0.5) +
  theme_minimal() +
  labs(
    title = "Volcano Plot: Infected vs Clean",
    x = "log2(Fold Change)",
    y = "-log10(adjusted p-value)",
    color = "Status"
  )

ggsave(file.path(de_dir, "02_volcano_plot.png"), p_volcano, width = 10, height = 8, dpi = 300)
cat("✓ Salvo: 02_volcano_plot.png (versão ggplot2)\n")

# ============================================================================
# 8. ORDENAR E VISUALIZAR
# ============================================================================

cat("\n--- Top genes diferencialmente expressos ---\n")

# Converter para data.frame e remover NAs
res_df <- as.data.frame(res) %>%
  rownames_to_column("gene_id") %>%
  filter(!is.na(padj))

# Ordenar por p-value ajustado
res_ordered <- res_df %>% arrange(padj)

# Top 20 genes significativos (padj < 0.05)
top_genes <- res_ordered %>%
  filter(padj < 0.05) %>%
  head(20)

cat("\nTop 20 genes significativos (padj < 0.05):\n")
print(top_genes)

# Top upregulated em infected
top_up <- res_ordered %>%
  filter(log2FoldChange > 0, padj < 0.05) %>%
  head(10)

cat("\nTop 10 upregulated em Infected:\n")
print(top_up)

# Top downregulated em infected
top_down <- res_ordered %>%
  filter(log2FoldChange < 0, padj < 0.05) %>%
  head(10)

cat("\nTop 10 downregulated em Infected:\n")
print(top_down)

# ============================================================================
# 9. HEATMAP DOS TOP GENES
# ============================================================================

# ============================================================================
# 10. HEATMAP DOS TOP GENES
# ============================================================================

cat("\n--- Gerando heatmap dos top genes ---\n")

# Selecionar top 20 genes significativos
top_20_genes <- res_ordered %>%
  filter(padj < 0.05) %>%
  head(20) %>%
  pull(gene_id)

# Extrair contagens normalizadas (rlog) para esses genes
rld <- rlog(dds, blind = FALSE)
top_genes_counts <- assay(rld)[top_20_genes, ]

# Criar anotação de colunas a partir de colData(dds)
# IMPORTANTE: usar os nomes EXATOS de colnames(dds)
metadata_dds <- as.data.frame(colData(dds))
metadata_dds$sample_name <- colnames(dds)

annotation_col <- metadata_dds %>%
  select(status, variety)

# Garantir correspondência
stopifnot(rownames(annotation_col) == colnames(top_genes_counts))

# Criar diretório para gráficos
dir.create("results/deseq2/plots", showWarnings = FALSE)

# Gerar heatmap
png("results/deseq2/plots/01_heatmap_top20_genes.png", width = 10, height = 8, units = "in", res = 300)
pheatmap(
  top_genes_counts,
  annotation_col = annotation_col,
  scale = "row",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Heatmap: Top 20 genes diferencialmente expressos",
  fontsize = 10,
  color = colorRampPalette(c("blue", "white", "red"))(50)
)
dev.off()

cat("✓ Salvo: 01_heatmap_top20_genes.png\n")

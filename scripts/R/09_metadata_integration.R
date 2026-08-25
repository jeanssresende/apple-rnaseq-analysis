################################################################################
# Módulo 09 — Integração de Metadata Experimental
################################################################################
#
# Objetivo:
#   Integrar o arquivo Metadata.csv ao colData do objeto gse (gene-level),
#   fazendo parsing dos nomes de arquivo e validando correspondência 100%.
#
# Entrada:
#   - results/tximeta/gse_gene_annotated.rds (objeto DESeq2-ready)
#   - metadata/Metadata.csv (metadados experimentais estruturados)
#
# Saída:
#   - results/tximeta/gse_with_metadata.rds (gse com colData atualizado)
#   - metadata/sample_metadata_processed.csv (metadata processado e validado)
#
# Autor: [Jean]
# Data: 2026-08-24
#
################################################################################

# ============================================================================
# 1. SETUP E CARREGAMENTO
# ============================================================================

library(SummarizedExperiment)
library(tidyverse)

# Definir diretórios
project_root <- here::here()
results_dir <- file.path(project_root, "results", "tximeta")
metadata_dir <- file.path(project_root, "metadata")

# Verificar arquivos de entrada
gse_file <- file.path(results_dir, "gse_gene_annotated.rds")
metadata_file <- file.path(metadata_dir, "Metadata.csv")

if (!file.exists(gse_file)) {
  stop("Arquivo não encontrado: ", gse_file)
}

if (!file.exists(metadata_file)) {
  stop("Arquivo não encontrado: ", metadata_file)
}

cat("Arquivos de entrada validados\n")

# ============================================================================
# 2. CARREGAR OBJETOS
# ============================================================================

cat("\n--- Carregando objetos ---\n")

# Carregar gse
gse <- readRDS(gse_file)
cat("✓ gse carregado:", dim(gse), "genes x", ncol(gse), "amostras\n")

# Carregar metadata
metadata <- read.csv(metadata_file, stringsAsFactors = FALSE)
cat("✓ Metadata carregado:", nrow(metadata), "amostras\n")

# ============================================================================
# 3. INSPECIONAR ESTRUTURA ATUAL
# ============================================================================

cat("\n--- Estrutura atual do gse ---\n")
cat("Nomes de amostra (primeiros 5):\n")
print(colnames(gse)[1:5])

cat("\ncolData atual:\n")
print(head(colData(gse), 3))

cat("\nEstrutura do metadata:\n")
print(head(metadata, 3))
print(paste("Colunas:", paste(colnames(metadata), collapse = ", ")))

# ============================================================================
# 4. PARSING DOS NOMES DE ARQUIVO
# ============================================================================

cat("\n--- Parsing dos nomes de amostra ---\n")

# Extrair componentes dos nomes de arquivo
# Formato: {Sample_ID}_{Status}_{Lane}
# Exemplo: 2H50475B2_rep1_clean_S21
# Exemplo: H50475B1_3_rep1_infected_S26 (com underscore extra no Sample_ID)

sample_names <- colnames(gse)

# Usar regex para extrair componentes
# Estratégia: extrair tudo antes de _(clean|infected)_S
parsed <- tibble(
  full_name = sample_names,
  # Extrair Sample_ID (tudo antes de _clean ou _infected)
  # Funciona para H50475B1_3_rep1 também
  sample_id = str_remove(full_name, "_(clean|infected)_S\\d+$"),
  # Extrair Status (clean ou infected)
  status = str_extract(full_name, "(?<=_)(clean|infected)(?=_)"),
  # Extrair Lane (S + número)
  lane = str_extract(full_name, "S\\d+$")
)

cat("Parsing resultado (primeiras 10):\n")
print(head(parsed, 10))

# Validar parsing
if (any(is.na(parsed$sample_id))) {
  cat("ERRO: Alguns sample_id não foram extraídos\n")
  print(parsed %>% filter(is.na(sample_id)))
  stop("Falha no parsing. Verifique os nomes de arquivo.")
}

if (any(is.na(parsed$status))) {
  cat("Aviso: Alguns status não foram extraídos\n")
  print(parsed %>% filter(is.na(status)))
}

cat("✓ Parsing concluído com sucesso\n")

# ============================================================================
# 5. MATCH COM METADATA
# ============================================================================

cat("\n--- Fazendo match com metadata ---\n")

# Preparar metadata para join
metadata_clean <- metadata %>%
  select(Sample_ID, Variety, Status, Virus_combo, Group) %>%
  rename(sample_id = Sample_ID,
         variety = Variety,
         status_meta = Status,
         virus_combo = Virus_combo,
         group = Group)

# Converter status para lowercase para match
metadata_clean$status_lower <- tolower(metadata_clean$status_meta)

# Join
merged <- parsed %>%
  left_join(
    metadata_clean %>% select(sample_id, variety, virus_combo, group, status_lower),
    by = c("sample_id", "status" = "status_lower")
  )

cat("Merged (primeiras 5):\n")
print(head(merged, 5))

# ============================================================================
# 6. VALIDAÇÃO DE CORRESPONDÊNCIA
# ============================================================================

cat("\n--- Validação de correspondência ---\n")

# Verificar se há NAs após o join
na_count <- sum(is.na(merged$variety))
if (na_count > 0) {
  cat("ERRO: ", na_count, " amostras não fizeram match com metadata\n")
  cat("Amostras sem match:\n")
  print(merged %>% filter(is.na(variety)))
  stop("Falha na correspondência. Verifique os nomes de arquivo vs metadata.")
}

cat("✓ 100% de correspondência validada\n")
cat("✓ Todas as", nrow(merged), "amostras foram mapeadas corretamente\n")

# ============================================================================
# 7. CONSTRUIR colData ESTRUTURADO
# ============================================================================

cat("\n--- Construindo colData estruturado ---\n")

# Criar data frame com colData
new_coldata <- merged %>%
  select(full_name, sample_id, status, lane, variety, virus_combo, group) %>%
  column_to_rownames("full_name") %>%
  as.data.frame()

# Converter colunas para fatores apropriados
new_coldata$status <- factor(new_coldata$status, levels = c("clean", "infected"))
new_coldata$variety <- factor(new_coldata$variety, levels = c("V1", "V2", "V3"))
new_coldata$group <- factor(new_coldata$group)

cat("colData novo (primeiras 5):\n")
print(head(new_coldata, 5))

cat("\nResumo de fatores:\n")
cat("Status:\n")
print(table(new_coldata$status))
cat("\nVariety:\n")
print(table(new_coldata$variety))
cat("\nVirus_combo:\n")
print(table(new_coldata$virus_combo))

# ============================================================================
# 8. ATUALIZAR colData DO GSE
# ============================================================================

cat("\n--- Atualizando colData do gse ---\n")

# Garantir que a ordem está correta
stopifnot(rownames(new_coldata) == colnames(gse))

# Converter data.frame para DataFrame (classe Bioconductor S4)
library(S4Vectors)
colData(gse) <- S4Vectors::DataFrame(new_coldata)

cat("✓ colData atualizado\n")
cat("Estrutura final do gse:\n")
print(gse)

# ============================================================================
# 9. SALVAR METADATA PROCESSADO (VERSIONÁVEL)
# ============================================================================

cat("\n--- Salvando metadata processado ---\n")

# Converter new_coldata para CSV (mantendo rownames como coluna)
metadata_processed <- new_coldata %>%
  rownames_to_column("sample_name")

output_metadata_file <- file.path(metadata_dir, "sample_metadata_processed.csv")
write.csv(metadata_processed, output_metadata_file, row.names = FALSE)
cat("✓ Metadata processado salvo em:", output_metadata_file, "\n")

# ============================================================================
# 10. SALVAR GSE ATUALIZADO
# ============================================================================

cat("\n--- Salvando gse atualizado ---\n")

output_gse_file <- file.path(results_dir, "gse_with_metadata.rds")
saveRDS(gse, output_gse_file)
cat("✓ gse com metadata salvo em:", output_gse_file, "\n")

# ============================================================================
# 11. RESUMO FINAL
# ============================================================================

cat("\n", strrep("=", 80), "\n")
cat("✓ MÓDULO 09 CONCLUÍDO COM SUCESSO\n")
cat(strrep("=", 80), "\n\n")

cat("RESUMO:\n")
cat("- Amostras processadas:", ncol(gse), "\n")
cat("- Genes:", nrow(gse), "\n")
cat("- Cultivares:", nlevels(gse$variety), "(" , paste(levels(gse$variety), collapse = ", "), ")\n")
cat("- Condições:", nlevels(gse$status), "(", paste(levels(gse$status), collapse = ", "), ")\n")
cat("- Replicatas por grupo:\n")
print(table(gse$variety, gse$status))

cat("\nARQUIVOS GERADOS:\n")
cat("- Metadata processado: metadata/sample_metadata_processed.csv\n")
cat("- GSE atualizado: results/tximeta/gse_with_metadata.rds\n")

cat("\nPRÓXIMO MÓDULO:\n")
cat("Módulo 10 — Validação dos Objetos Importados\n")

################################################################################

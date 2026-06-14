# =============================================================================
# Random Forest classification of successional stages in subtropical
# Atlantic Forest (Santa Catarina, Brazil)
#
# Pipeline: VIF (collinearity removal) -> Boruta (variable selection) ->
#           PCoA (floristic gradients) -> train/test split ->
#           ntree optimisation -> Random Forest -> evaluation ->
#           comparison against CONAMA 04/94 thresholds
# =============================================================================

# -----------------------------------------------------------------------------
# Libraries
# -----------------------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(vegan)
library(ape)
library(randomForest)
library(caret)
library(ggplot2)
library(Boruta)
library(pROC)
library(corrplot)
library(car)
library(janitor)
library(pdp)
library(gridExtra)
library(cowplot)

# Reproducibility
set.seed(42)


# =============================================================================
# 1. DATA IMPORT
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 1: DATA IMPORT\n")
cat(rep("=", 70), "\n\n")

florestasc_arvores     <- read_excel("data/florestasc_arvores.xlsx")
florestasc_metadados   <- read_excel("data/florestasc_metadados.xlsx") %>% clean_names()
florestasc_parametros  <- read_excel("data/florestasc_parametros.xlsx") %>% clean_names()

# Standardise the plot identifier across tables
names(florestasc_parametros)[1] <- "plot"
names(florestasc_metadados)[1]   <- "plot"

df_complete <- left_join(florestasc_parametros, florestasc_metadados, by = "plot")

variaveis_conama <- df_complete %>%
  select(plot, ab_ha, ht_med, dap_med, estagio_sucessional)

# Response variable as an ordered factor (early -> medium -> advanced)
df_complete <- df_complete %>%
  mutate(estagio_sucessional = factor(estagio_sucessional,
                                      levels = c("Inicial", "Médio", "Avançado")))

cat("Plots imported:", nrow(df_complete), "\n")
cat("Total variables:", ncol(df_complete) - 1, "(+ response)\n\n")

# Drop categorical/geographic fields not used as predictors
df_complete <- df_complete[, !(names(df_complete) %in%
  c("ab", "s_reg", "classe", "x_longitude", "y_latitude",
    "municipio", "bacia_hidrografica", "regiao_fitoecologica"))]

# Preserve plot IDs before any row-altering transformation
plot_ids <- data.frame(
  row_id = seq_len(nrow(df_complete)),
  plot   = df_complete$plot
)

cat("Plot IDs preserved:", nrow(plot_ids), "records\n\n")


# =============================================================================
# 2. PRE-PROCESSING
# =============================================================================
cat(rep("=", 70), "\n")
cat("STEP 2: PRE-PROCESSING\n")
cat(rep("=", 70), "\n\n")

# Inspect missing values before imputation
na_count          <- colSums(is.na(df_complete))
variaveis_com_na  <- na_count[na_count > 0]
total_nas         <- sum(na_count)

if (total_nas > 0) {
  cat("Missing values detected:\n")
  cat("  Total NAs:", total_nas, "\n")
  cat("  Affected variables:", length(variaveis_com_na), "\n\n")
  cat("Breakdown by variable:\n")
  for (var in names(variaveis_com_na)) {
    percentual <- round((variaveis_com_na[var] / nrow(df_complete)) * 100, 2)
    cat(sprintf("  %s: %d NAs (%.2f%%)\n", var, variaveis_com_na[var], percentual))
  }
  cat("\n")
} else {
  cat("No missing values detected\n\n")
}

# Median imputation for numeric predictors
df_preprocessado <- df_complete %>%
  mutate(across(where(is.numeric), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))

cat("Missing values imputed with the median\n")
cat("Dataset ready:", ncol(df_preprocessado) - 1, "variables\n\n")

# Drop the plot identifier before VIF/Boruta
df_preprocessado_sem_plot <- df_preprocessado %>% select(-plot)


# =============================================================================
# 3. COLLINEARITY REMOVAL (VIF)
# =============================================================================
cat(rep("=", 70), "\n")
cat("STEP 3: VIF ANALYSIS AND COLLINEARITY REMOVAL\n")
cat(rep("=", 70), "\n\n")

# Compute the Variance Inflation Factor for every numeric predictor
calcular_vif <- function(dados, variavel_resposta) {

  dados_num <- dados %>%
    select(-all_of(variavel_resposta)) %>%
    select(where(is.numeric))

  # Keep only variables with non-zero, defined variance
  vars_validas <- sapply(dados_num, function(x) {
    var_x <- var(x, na.rm = TRUE)
    !is.na(var_x) && var_x > 0
  })
  dados_num <- dados_num[, vars_validas]

  if (ncol(dados_num) < 2) {
    cat("Fewer than two numeric variables available for VIF\n")
    return(NULL)
  }

  tryCatch({
    dados_completos <- dados_num %>% na.omit()

    if (nrow(dados_completos) < ncol(dados_completos) + 1) {
      cat("Insufficient data to compute VIF\n")
      return(NULL)
    }

    vif_valores <- numeric(ncol(dados_completos))
    names(vif_valores) <- names(dados_completos)

    for (i in seq_len(ncol(dados_completos))) {
      formula_temp <- as.formula(paste(names(dados_completos)[i], "~ ."))
      modelo_temp  <- lm(formula_temp, data = dados_completos)
      r_squared    <- summary(modelo_temp)$r.squared
      vif_valores[i] <- 1 / (1 - r_squared)
    }

    data.frame(
      Variavel = names(vif_valores),
      VIF      = as.numeric(vif_valores),
      stringsAsFactors = FALSE
    ) %>%
      mutate(
        Interpretacao = case_when(
          VIF < 5            ~ "Low",
          VIF >= 5 & VIF < 10 ~ "Moderate",
          VIF >= 10          ~ "High"
        )
      ) %>%
      arrange(desc(VIF))

  }, error = function(e) {
    cat("Error computing VIF:", e$message, "\n")
    return(NULL)
  })
}

# Iterative, stepwise collinearity pruning:
# remove the highest-VIF predictor first; once all VIF < threshold,
# remove one member of any pair exceeding the correlation threshold.
remover_colinearidade <- function(dados, variavel_resposta,
                                  vif_threshold = 10,
                                  cor_threshold = 0.8,
                                  max_iter = 20) {

  cat("Starting automatic collinearity removal\n")
  cat("  Maximum allowed VIF:", vif_threshold, "\n")
  cat("  Maximum allowed correlation:", cor_threshold, "\n\n")

  dados_limpos   <- dados
  vars_removidas <- character(0)
  iter <- 0

  repeat {
    iter <- iter + 1
    if (iter > max_iter) {
      cat("Maximum of", max_iter, "iterations reached\n")
      break
    }

    vif_atual <- calcular_vif(dados_limpos, variavel_resposta)
    if (is.null(vif_atual)) break

    vif_alto <- vif_atual %>% filter(VIF > vif_threshold)

    if (nrow(vif_alto) == 0) {
      dados_num <- dados_limpos %>%
        select(-all_of(variavel_resposta)) %>%
        select(where(is.numeric))

      if (ncol(dados_num) > 1) {
        matriz_cor    <- cor(dados_num, use = "complete.obs")
        matriz_cor_df <- as.data.frame(as.table(matriz_cor))
        names(matriz_cor_df) <- c("Var1", "Var2", "Correlacao")

        cor_altas <- matriz_cor_df %>%
          filter(abs(Correlacao) > cor_threshold & abs(Correlacao) < 1) %>%
          filter(as.character(Var1) < as.character(Var2))

        if (nrow(cor_altas) > 0) {
          par_max     <- cor_altas %>% arrange(desc(abs(Correlacao))) %>% slice(1)
          var_remover <- as.character(par_max$Var2)
          cat("  Iter", iter, "-> removing", var_remover,
              "(cor =", round(par_max$Correlacao, 3), ")\n")
          dados_limpos   <- dados_limpos %>% select(-all_of(var_remover))
          vars_removidas <- c(vars_removidas, var_remover)
        } else {
          cat("\nCleaning complete\n")
          break
        }
      } else {
        break
      }
    } else {
      var_remover <- vif_alto %>% arrange(desc(VIF)) %>% slice(1) %>% pull(Variavel)
      vif_valor   <- vif_alto %>% filter(Variavel == var_remover) %>% pull(VIF)
      cat("  Iter", iter, "-> removing", var_remover,
          "(VIF =", round(vif_valor, 2), ")\n")
      dados_limpos   <- dados_limpos %>% select(-all_of(var_remover))
      vars_removidas <- c(vars_removidas, var_remover)
    }
  }

  list(dados = dados_limpos, removidas = vars_removidas)
}

# Run VIF cleaning (without the plot identifier)
resultado_vif <- remover_colinearidade(df_preprocessado_sem_plot, "estagio_sucessional",
                                       vif_threshold = 10,
                                       cor_threshold = 0.8)

df_pos_vif <- resultado_vif$dados

cat("\nVIF summary:\n")
cat("  Variables before:", ncol(df_preprocessado_sem_plot) - 1, "\n")
cat("  Variables after:", ncol(df_pos_vif) - 1, "\n")
cat("  Variables removed:", length(resultado_vif$removidas), "\n")

if (length(resultado_vif$removidas) > 0) {
  cat("\nVariables removed by VIF:\n")
  cat(paste("  -", resultado_vif$removidas, collapse = "\n"), "\n")
}

# Final VIF diagnostics
cat("\nFinal VIF diagnostics:\n")
cat("  Numeric variables in df_pos_vif:",
    sum(sapply(df_pos_vif %>% select(-estagio_sucessional), is.numeric)), "\n")

vars_numericas <- df_pos_vif %>%
  select(-estagio_sucessional) %>%
  select(where(is.numeric))

cat("  Zero-variance variables:\n")
vars_zero <- sapply(vars_numericas,
                    function(x) var(x, na.rm = TRUE) == 0 | is.na(var(x, na.rm = TRUE)))
if (any(vars_zero)) {
  cat("   ", names(vars_zero[vars_zero]), "\n")
} else {
  cat("    None\n")
}

vif_final <- calcular_vif(df_pos_vif, "estagio_sucessional")

if (is.null(vif_final)) {
  cat("\nWARNING: final VIF could not be computed\n\n")

} else if (any(is.na(vif_final$VIF)) || any(is.nan(vif_final$VIF)) ||
           any(is.infinite(vif_final$VIF))) {
  cat("\nWARNING: final VIF contains invalid values (NA/NaN/Inf)\n")
  print(vif_final[is.na(vif_final$VIF) | is.nan(vif_final$VIF) |
                  is.infinite(vif_final$VIF), ])

  vif_final_valido <- vif_final %>%
    filter(!is.na(VIF) & !is.nan(VIF) & !is.infinite(VIF))

  if (nrow(vif_final_valido) > 0) {
    cat("\nMaximum VIF (valid values):", round(max(vif_final_valido$VIF), 2), "\n\n")

    p_vif <- ggplot(vif_final_valido,
                    aes(x = reorder(Variavel, VIF), y = VIF, fill = Interpretacao)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 5,  linetype = "dashed", color = "orange", linewidth = 1) +
      geom_hline(yintercept = 10, linetype = "dashed", color = "red",    linewidth = 1) +
      coord_flip() +
      scale_fill_manual(values = c("Baixa" = "green3",
                                   "Moderada" = "orange",
                                   "Alta" = "red2")) +
      labs(title = "Final VIF (valid values)", x = "Variables", y = "VIF") +
      theme_minimal() +
      theme(legend.position = "bottom", axis.text.y = element_text(size = 8))

    print(p_vif)
  }

} else {
  cat("\nMaximum final VIF:", round(max(vif_final$VIF), 2), "\n\n")

  p_vif <- ggplot(vif_final,
                  aes(x = reorder(Variavel, VIF), y = VIF, fill = Interpretacao)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.5, width = 0.7) +
    geom_hline(yintercept = 5,  linetype = "dashed", color = "gray20", linewidth = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "gray20", linewidth = 0.7) +
    coord_flip() +
    scale_fill_manual(values = c("Low" = "gray70", "Moderate" = "gray40", "High" = "black")) +
    labs(x = NULL, y = "VIF", fill = "Status") +
    theme_classic(base_size = 11) +
    theme(
      legend.position    = "top",
      legend.background  = element_rect(fill = "white", color = "black", linewidth = 0.3),
      legend.title       = element_text(size = 10, face = "bold"),
      legend.text        = element_text(size = 9),
      axis.title         = element_text(size = 11, face = "bold"),
      axis.text          = element_text(size = 10, color = "black"),
      axis.line          = element_line(color = "black", linewidth = 0.5),
      axis.ticks         = element_line(color = "black", linewidth = 0.5),
      panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.major.y = element_blank()
    )

  print(p_vif)

  ggsave("figura_vif_final.tiff", plot = p_vif, width = 14, height = 10,
         units = "cm", dpi = 300, compression = "lzw")
  cat("figura_vif_final.tiff exported (14x10 cm, 300 dpi)\n")
}


# =============================================================================
# 4. VARIABLE SELECTION (BORUTA)
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 4: VARIABLE SELECTION (BORUTA)\n")
cat(rep("=", 70), "\n\n")

cat("Running Boruta (this may take a while)...\n\n")
boruta_result <- Boruta(estagio_sucessional ~ ., data = df_pos_vif, doTrace = 2)
boruta_result <- TentativeRoughFix(boruta_result)
stats_boruta  <- attStats(boruta_result)
cat("\nBoruta complete\n\n")
print(boruta_result)

# Importance plot (long predictor names rotated for readability)
par(mar = c(12, 5, 4, 2) + 0.1)
plot(boruta_result,
     las = 2, cex.axis = 0.7,
     xlab = "", ylab = "Importance", main = "")
par(mar = c(5, 4, 4, 2) + 0.1)

# Importance history across iterations
par(mar = c(12, 5, 4, 2) + 0.1)
plotImpHistory(boruta_result, ylab = "Importance", xlab = "")
legend("topright",
       legend = c("Important", "Not Important", "Undecided", "Shadow"),
       col = c("green", "red", "yellow", "black"),
       lty = 1, lwd = 2, cex = 0.8)
par(mar = c(5, 4, 4, 2) + 0.1)

rejected_vars  <- rownames(stats_boruta[stats_boruta$decision == "Rejected", ])
confirmed_vars <- rownames(stats_boruta[stats_boruta$decision == "Confirmed", ])
df_pos_boruta  <- df_pos_vif[, !(names(df_pos_vif) %in% rejected_vars)]

cat("\nBoruta summary:\n")
cat("  Confirmed variables:", length(confirmed_vars), "\n")
cat("  Rejected variables:", length(rejected_vars), "\n")
cat("  Final variables:", ncol(df_pos_boruta) - 1, "\n")

if (length(rejected_vars) > 0) {
  cat("\nVariables rejected by Boruta:\n")
  cat(paste("  -", rejected_vars, collapse = "\n"), "\n")
}


# =============================================================================
# 4.5 FLORISTIC GRADIENTS (PCoA)
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 4.5: PCoA INTEGRATION\n")
cat(rep("=", 70), "\n\n")

cat("Computing PCoA...\n\n")

# Species-by-plot abundance matrix
matriz_especies <- florestasc_arvores %>%
  group_by(across(c(1, 7))) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = 2, values_from = n, values_fill = 0) %>%
  as.data.frame()

rownames(matriz_especies) <- matriz_especies[, 1]
matriz_especies <- matriz_especies[, -1]

cat("Species matrix:", nrow(matriz_especies), "plots x",
    ncol(matriz_especies), "species\n")

# Bray-Curtis dissimilarity and PCoA ordination (first 10 axes)
dist_matrix <- vegdist(matriz_especies, method = "bray")
pcoa_result <- cmdscale(dist_matrix, k = 10, eig = TRUE)

eigenvalues <- pcoa_result$eig[pcoa_result$eig > 0]
var_exp     <- eigenvalues / sum(eigenvalues) * 100

cat("\nVariance explained by the 10 PCoA axes:\n")
for (i in 1:10) {
  cat(sprintf("  PCoA%d: %.2f%%\n", i, var_exp[i]))
}
cat(sprintf("\n  Cumulative: %.2f%%\n", sum(var_exp[1:10])))

pcoa_scores <- as.data.frame(pcoa_result$points[, 1:10])
colnames(pcoa_scores) <- paste0("PCoA", 1:10)
pcoa_scores$plot <- as.numeric(rownames(pcoa_scores))

# Merge PCoA axes back through the preserved plot IDs
df_pos_boruta_pcoa <- df_pos_boruta %>%
  mutate(row_id = row_number()) %>%
  left_join(plot_ids, by = "row_id") %>%
  left_join(pcoa_scores, by = "plot") %>%
  select(-row_id, -plot)

cat("\nPCoA integrated into the dataset\n")
cat("  Final variables:", ncol(df_pos_boruta_pcoa) - 1,
    "(", ncol(df_pos_boruta) - 1, "+ 10 PCoA)\n")

pcoa_cols <- paste0("PCoA", 1:10)
na_count  <- sum(is.na(df_pos_boruta_pcoa[, pcoa_cols]))

if (na_count > 0) {
  cat("\nImputing", na_count, "NAs in the PCoA axes...\n")
  df_pos_boruta_pcoa <- df_pos_boruta_pcoa %>%
    mutate(across(all_of(pcoa_cols), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))
}

# PCoA ordination plot (first two axes)
pcoa_viz <- as.data.frame(pcoa_result$points[, 1:2])
colnames(pcoa_viz) <- c("PCoA1", "PCoA2")

p_pcoa <- ggplot(pcoa_viz, aes(x = PCoA1, y = PCoA2)) +
  geom_point(size = 4, alpha = 0.7, color = "steelblue") +
  geom_text(aes(label = rownames(pcoa_viz)), vjust = -1, size = 3) +
  labs(
    title = "PCoA - species composition",
    x = paste0("PCoA1 (", round(var_exp[1], 1), "%)"),
    y = paste0("PCoA2 (", round(var_exp[2], 1), "%)")
  ) +
  theme_minimal()

print(p_pcoa)

ggsave("figura_pcoa_composicao.tiff", plot = p_pcoa, width = 12, height = 10,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_pcoa_composicao.tiff exported (12x10 cm, 300 dpi)\n")


# =============================================================================
# 5. TRAIN/TEST SPLIT
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 5: TRAIN/TEST SPLIT\n")
cat(rep("=", 70), "\n\n")

split  <- createDataPartition(df_pos_boruta_pcoa$estagio_sucessional, p = 0.7, list = FALSE)
treino <- df_pos_boruta_pcoa[split, ]
teste  <- df_pos_boruta_pcoa[-split, ]

cat("Training set:", nrow(treino), "plots (70%)\n")
cat("Test set:", nrow(teste), "plots (30%)\n")
cat("Variables (including PCoA):", ncol(treino) - 1, "\n\n")


# =============================================================================
# 5.5 NTREE OPTIMISATION VIA OOB ERROR
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 5.5: NTREE OPTIMISATION\n")
cat(rep("=", 70), "\n\n")

cat("Testing different ntree values (OOB error on the training set)\n\n")

ntree_seq  <- seq(50, 500, by = 10)
oob_errors <- numeric(length(ntree_seq))

cat("Progress:\n")
pb <- txtProgressBar(min = 0, max = length(ntree_seq), style = 3)

for (i in seq_along(ntree_seq)) {
  rf_temp <- randomForest(
    estagio_sucessional ~ .,
    data       = treino,
    ntree      = ntree_seq[i],
    mtry       = 5,
    importance = FALSE,
    proximity  = FALSE
  )
  oob_errors[i] <- rf_temp$err.rate[ntree_seq[i], "OOB"]
  setTxtProgressBar(pb, i)
}
close(pb)

oob_df <- data.frame(ntree = ntree_seq, oob_error = oob_errors)

# Minimum-error and stabilisation criteria
ntree_otimo  <- oob_df$ntree[which.min(oob_df$oob_error)]
erro_minimo  <- min(oob_df$oob_error)

ultimos_10pct <- tail(oob_df, n = ceiling(nrow(oob_df) * 0.1))
threshold     <- mean(ultimos_10pct$oob_error)
ntree_estavel <- oob_df$ntree[min(which(oob_df$oob_error <= threshold))]

cat("\n\nOOB analysis results:\n")
cat(sprintf("  Minimum OOB error: %.4f (ntree = %d)\n", erro_minimo, ntree_otimo))
cat(sprintf("  Stabilisation: ntree = %d (error <= %.4f)\n", ntree_estavel, threshold))
cat(sprintf("  Error reduction: %.2f%% (50 -> %d trees)\n",
            (oob_df$oob_error[1] - erro_minimo) / oob_df$oob_error[1] * 100,
            ntree_otimo))

ntree_recomendado <- ntree_estavel
cat(sprintf("\nRecommended ntree: %d\n", ntree_recomendado))

p_oob <- ggplot(oob_df, aes(x = ntree, y = oob_error)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "steelblue", size = 2, alpha = 0.6) +
  geom_vline(xintercept = ntree_otimo,   linetype = "dashed", color = "red",       linewidth = 1) +
  geom_vline(xintercept = ntree_estavel, linetype = "dashed", color = "darkgreen", linewidth = 1) +
  geom_hline(yintercept = threshold,     linetype = "dotted", color = "gray50",    linewidth = 0.8) +
  annotate("text", x = ntree_otimo, y = max(oob_df$oob_error) * 0.95,
           label = paste("Minimum\n(ntree =", ntree_otimo, ")"),
           hjust = -0.1, color = "red", size = 3.5, fontface = "bold") +
  annotate("text", x = ntree_estavel, y = max(oob_df$oob_error) * 0.85,
           label = paste("Stabilisation\n(ntree =", ntree_estavel, ")"),
           hjust = -0.1, color = "darkgreen", size = 3.5, fontface = "bold") +
  labs(
    title    = "OOB error vs number of trees",
    subtitle = paste("Based on", nrow(treino), "training plots"),
    x = "Number of trees (ntree)",
    y = "OOB error"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank()
  )

print(p_oob)

ggsave("figura_oob_ntree.tiff", plot = p_oob, width = 14, height = 9,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_oob_ntree.tiff exported (14x9 cm, 300 dpi)\n")

# Marginal improvement between consecutive ntree values
oob_df <- oob_df %>%
  mutate(
    melhoria     = c(NA, -diff(oob_error)),
    melhoria_pct = (melhoria / lag(oob_error)) * 100
  )

cat("\nLast 10 iterations - improvement rate:\n")
print(
  tail(oob_df[, c("ntree", "oob_error", "melhoria", "melhoria_pct")], 10) %>%
    mutate(melhoria = round(melhoria, 6), melhoria_pct = round(melhoria_pct, 4))
)

p_melhoria <- ggplot(oob_df[-1, ], aes(x = ntree, y = melhoria_pct)) +
  geom_line(color = "coral", linewidth = 1) +
  geom_point(color = "coral", size = 1.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(
    title    = "OOB error improvement rate",
    subtitle = "Percentage change between consecutive iterations",
    x = "Number of trees (ntree)",
    y = "Improvement (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.minor = element_blank()
  )

print(p_melhoria)

ggsave("figura_oob_melhoria.tiff", plot = p_melhoria, width = 14, height = 9,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_oob_melhoria.tiff exported (14x9 cm, 300 dpi)\n")

cat(sprintf("\nUse ntree = %d in the final model\n\n", ntree_recomendado))

ntree_final <- ntree_recomendado


# =============================================================================
# 6. RANDOM FOREST MODEL
# =============================================================================
cat(rep("=", 70), "\n")
cat("STEP 6: RANDOM FOREST TRAINING\n")
cat(rep("=", 70), "\n\n")

modelo_rf <- randomForest(estagio_sucessional ~ .,
                          data       = treino,
                          ntree      = ntree_final,
                          mtry       = 5,
                          importance = TRUE,
                          proximity  = TRUE)

print(modelo_rf)


# =============================================================================
# 7. MODEL EVALUATION
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 7: EVALUATION - TRAINING AND TEST\n")
cat(rep("=", 70), "\n\n")

pred_treino <- predict(modelo_rf, newdata = treino)
conf_treino <- confusionMatrix(pred_treino, treino$estagio_sucessional)

cat("TRAINING SET:\n")
print(conf_treino)

pred_teste <- predict(modelo_rf, newdata = teste)
conf_teste <- confusionMatrix(pred_teste, teste$estagio_sucessional)

cat("\nTEST SET:\n")
print(conf_teste)

# Confusion-matrix tile plot
plot_confusion_matrix <- function(confusion_matrix, title) {

  cm_df <- as.data.frame(confusion_matrix$table) %>%
    group_by(Reference) %>%
    mutate(
      Total      = sum(Freq),
      Percentage = round((Freq / Total) * 100, 1),
      Label      = paste0(Freq, "\n(", Percentage, "%)")
    ) %>%
    ungroup()

  ggplot(cm_df, aes(x = Reference, y = Prediction)) +
    geom_tile(aes(fill = Freq), color = "white", linewidth = 1) +
    geom_text(aes(label = Label), size = 5, fontface = "bold") +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Frequency") +
    labs(title = title, x = "True class", y = "Predicted class") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank()
    ) +
    coord_fixed()
}

p_treino <- plot_confusion_matrix(conf_treino, "Confusion matrix - training")
p_teste  <- plot_confusion_matrix(conf_teste,  "Confusion matrix - test")

print(p_treino)
ggsave("figura_conf_treino.tiff", plot = p_treino, width = 10, height = 10,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_conf_treino.tiff exported (10x10 cm, 300 dpi)\n")

print(p_teste)
ggsave("figura_conf_teste.tiff", plot = p_teste, width = 10, height = 10,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_conf_teste.tiff exported (10x10 cm, 300 dpi)\n")


# =============================================================================
# 8. ROC CURVES AND AUC
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 8: ROC CURVES AND AUC\n")
cat(rep("=", 70), "\n\n")

probs    <- predict(modelo_rf, newdata = teste, type = "prob")
auc_list <- multiclass.roc(teste$estagio_sucessional, probs)
print(auc_list)

roc_list <- list()
classes  <- levels(teste$estagio_sucessional)
for (cls in classes) {
  roc_list[[cls]] <- roc(response = ifelse(teste$estagio_sucessional == cls, 1, 0),
                         predictor = probs[, cls])
}

# Class labels for display (data in PT, figures in EN)
labels_classe <- c(
  "Inicial"  = "Early",
  "Médio"    = "Intermediate",
  "Avançado" = "Advanced"
)

roc_df <- do.call(rbind, lapply(seq_along(classes), function(i) {
  cls <- classes[i]
  r   <- roc_list[[cls]]
  data.frame(
    Especificidade = 1 - r$specificities,
    Sensibilidade  = r$sensitivities,
    Classe         = paste0(labels_classe[cls], " (AUC = ", round(auc(r), 3), ")"),
    stringsAsFactors = FALSE
  )
}))

p_roc <- ggplot(roc_df, aes(x = Especificidade, y = Sensibilidade,
                            color = Classe, linetype = Classe)) +
  geom_line(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("#E69F00", "#0072B2", "#009E73"), name = NULL) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash"), name = NULL) +
  guides(
    color    = guide_legend(nrow = 3, byrow = TRUE),
    linetype = guide_legend(nrow = 3, byrow = TRUE)
  ) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic(base_size = 10) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 8, family = "serif"),
    legend.key.width = unit(1.2, "cm"),
    axis.title       = element_text(size = 9, face = "bold", family = "serif"),
    axis.text        = element_text(size = 8, color = "black", family = "serif"),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black", linewidth = 0.4),
    panel.grid.major = element_line(color = "gray88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin      = margin(t = 4, r = 6, b = 2, l = 4, unit = "mm")
  )

print(p_roc)

ggsave("figura_roc_auc.tiff", plot = p_roc, width = 10, height = 11.5,
       units = "cm", dpi = 300, compression = "lzw")
cat("figura_roc_auc.tiff exported (10x11.5 cm, 300 dpi)\n")


# =============================================================================
# 9. PARTIAL DEPENDENCE PLOTS
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 9: PARTIAL DEPENDENCE PLOTS\n")
cat(rep("=", 70), "\n\n")

# Okabe & Ito (2008) colourblind-safe palette
cores_estagios <- c(
  "Inicial"  = "#E69F00",
  "Médio"    = "#0072B2",
  "Avançado" = "#009E73"
)

# Top 8 predictors by Mean Decrease Gini
importancia <- importance(modelo_rf)
importancia_df <- data.frame(
  Variavel         = rownames(importancia),
  MeanDecreaseGini = importancia[, "MeanDecreaseGini"]
) %>% arrange(desc(MeanDecreaseGini))

top8_vars <- head(importancia_df$Variavel, 8)

cat("Top 8 most important variables:\n")
for (i in 1:8) {
  cat(sprintf("  %2d. %s (%.2f)\n", i, top8_vars[i], importancia_df$MeanDecreaseGini[i]))
}

# Theme for individual panels (no legend)
tema_pdp_sem_legenda <- theme_classic(base_size = 9) +
  theme(
    legend.position    = "none",
    axis.title         = element_text(size = 8.5, face = "bold", family = "serif"),
    axis.text          = element_text(size = 8,   color = "black", family = "serif"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "gray88", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(t = 4, r = 6, b = 2, l = 4, unit = "mm")
  )

# Build one PDP panel for a given predictor (three class probabilities)
criar_pdp <- function(modelo, variavel) {

  pd_inicial  <- partial(modelo, pred.var = variavel, which.class = "Inicial",  prob = TRUE)
  pd_medio    <- partial(modelo, pred.var = variavel, which.class = "Médio",    prob = TRUE)
  pd_avancado <- partial(modelo, pred.var = variavel, which.class = "Avançado", prob = TRUE)

  pd_inicial$classe  <- "Inicial"
  pd_medio$classe    <- "Médio"
  pd_avancado$classe <- "Avançado"

  pd_all <- bind_rows(pd_inicial, pd_medio, pd_avancado) %>%
    rename(x_var = 1) %>%
    mutate(classe = factor(classe, levels = c("Inicial", "Médio", "Avançado")))

  ggplot(pd_all, aes(x = x_var, y = yhat, color = classe)) +
    geom_line(linewidth = 0.55) +
    scale_color_manual(values = cores_estagios, name = "Successional stage") +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.25),
      labels = scales::percent_format(accuracy = 1)
    ) +
    labs(x = variavel, y = "Predicted probability") +
    tema_pdp_sem_legenda
}

# Shared legend extracted from an auxiliary plot
tema_legenda <- theme_classic(base_size = 9) +
  theme(
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 9, face = "bold", family = "serif"),
    legend.text       = element_text(size = 8.5, family = "serif"),
    legend.key.size   = unit(0.4, "cm"),
    legend.key        = element_rect(fill = "white", color = NA),
    legend.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    legend.margin     = margin(t = 3, r = 8, b = 3, l = 8)
  )

p_aux_legenda <- ggplot(
  data.frame(x = 1, y = 1,
             classe = factor(c("Inicial", "Médio", "Avançado"),
                             levels = c("Inicial", "Médio", "Avançado"))),
  aes(x = x, y = y, color = classe)
) +
  geom_line(linewidth = 0.55) +
  scale_color_manual(
    values = cores_estagios,
    name   = "Successional stage",
    labels = c("Inicial" = "Early", "Médio" = "Medium", "Avançado" = "Advanced")
  ) +
  tema_legenda

legenda_compartilhada <- cowplot::get_legend(p_aux_legenda)

cat("\nGenerating PDPs (8 variables)...\n")

pdp_list <- lapply(seq_along(top8_vars), function(i) {
  cat("  Variable", i, ":", top8_vars[i], "\n")
  criar_pdp(modelo_rf, top8_vars[i])
})

cat("\nExporting PDP panels...\n\n")

# Four panels of two PDPs each, with a shared legend below
for (bloco in 1:4) {
  idx1 <- (bloco - 1) * 2 + 1
  idx2 <- idx1 + 1

  linha_graficos <- cowplot::plot_grid(
    pdp_list[[idx1]], pdp_list[[idx2]],
    ncol = 2, align = "hv", rel_widths = c(1, 1)
  )

  p_painel <- cowplot::plot_grid(
    linha_graficos, legenda_compartilhada,
    ncol = 1, rel_heights = c(1, 0.12)
  )

  nome_arq <- sprintf("figura_pdp_bloco%d.tiff", bloco)
  ggsave(nome_arq, plot = p_painel, width = 14, height = 8,
         units = "cm", dpi = 300, compression = "lzw")
  cat(sprintf("  %s saved (14x8 cm, 300 dpi)\n", nome_arq))
}

# Full panel: 4 rows x 2 columns plus a single legend
linha1 <- cowplot::plot_grid(pdp_list[[1]], pdp_list[[2]], ncol = 2, align = "hv")
linha2 <- cowplot::plot_grid(pdp_list[[3]], pdp_list[[4]], ncol = 2, align = "hv")
linha3 <- cowplot::plot_grid(pdp_list[[5]], pdp_list[[6]], ncol = 2, align = "hv")
linha4 <- cowplot::plot_grid(pdp_list[[7]], pdp_list[[8]], ncol = 2, align = "hv")

p_pdp_completo <- cowplot::plot_grid(
  linha1, linha2, linha3, linha4,
  legenda_compartilhada,
  ncol = 1, rel_heights = c(1, 1, 1, 1, 0.08)
)

ggsave("figura_pdp_completo.tiff", plot = p_pdp_completo, width = 14, height = 30,
       units = "cm", dpi = 300, compression = "lzw")
cat("  figura_pdp_completo.tiff saved (14x30 cm, 300 dpi)\n\n")

# =============================================================================
# 10. COMPARISON WITH CONAMA 04/94 THRESHOLDS
# Descriptive comparison performed on the full dataset (483 plots),
# contrasting observed and predicted structure against the legal thresholds.
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("STEP 10: COMPARISON WITH CONAMA 04/94\n")
cat(rep("=", 70), "\n\n")

# Predictions for the full dataset (training + test)
pred_completo <- predict(modelo_rf, newdata = df_pos_boruta_pcoa)

df_conama <- plot_ids %>%
  mutate(
    estagio_real    = df_complete$estagio_sucessional,
    estagio_predito = pred_completo
  ) %>%
  left_join(
    df_complete %>% select(plot, ab_ha, ht_med, dap_med),
    by = "plot"
  ) %>%
  pivot_longer(
    cols      = c(estagio_real, estagio_predito),
    names_to  = "origem",
    values_to = "estagio"
  ) %>%
  mutate(
    origem  = case_match(origem,
                         "estagio_real"    ~ "Real (Campo)",
                         "estagio_predito" ~ "Predito (Modelo)"),
    estagio = factor(estagio, levels = c("Inicial", "Médio", "Avançado"))
  )

cat("Dataset for the CONAMA comparison assembled\n")
cat("  Observations:", nrow(df_conama) / 2, "(field + predictions)\n\n")

# Box colours (Okabe & Ito)
cores_origem_box <- c(
  "Real (Campo)"     = "#0072B2",
  "Predito (Modelo)" = "#56B4E9"
)

# Source labels in English (data in PT, figures in EN)
labels_origem <- c(
  "Real (Campo)"     = "Observed (Field)",
  "Predito (Modelo)" = "Predicted (Model)"
)

labels_estagio <- c(
  "Inicial"  = "Early",
  "Médio"    = "Medium",
  "Avançado" = "Advanced"
)

# CONAMA threshold-line colours
cores_conama <- c(
  "Inicial" = "#D55E00",
  "Médio"   = "#CC79A7"
)

# CONAMA 04/94 upper thresholds for early and medium stages
conama_ref <- data.frame(
  variavel       = c("ab_ha",  "ab_ha",
                     "ht_med", "ht_med",
                     "dap_med","dap_med"),
  limite         = c(8,  15,
                     4,  12,
                     8,  15),
  rotulo_estagio = c("Inicial", "Médio",
                     "Inicial", "Médio",
                     "Inicial", "Médio"),
  rotulo_legenda = c("Early threshold (CONAMA)",  "Medium threshold (CONAMA)",
                     "Early threshold (CONAMA)",  "Medium threshold (CONAMA)",
                     "Early threshold (CONAMA)",  "Medium threshold (CONAMA)")
)

tema_conama <- theme_classic(base_size = 10) +
  theme(
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.box.spacing = unit(0.3, "cm"),
    legend.key.size    = unit(0.45, "cm"),
    legend.key.width   = unit(0.75, "cm"),
    legend.title       = element_text(size = 9, face = "bold",
                                      family = "serif", margin = margin(b = 4)),
    legend.text        = element_text(size = 8.5, family = "serif"),
    legend.spacing.y   = unit(0.25, "cm"),
    legend.background  = element_blank(),
    legend.margin      = margin(t = 7, r = 10, b = 7, l = 10),
    axis.title         = element_text(size = 9, face = "bold", family = "serif"),
    axis.text          = element_text(size = 8.5, color = "black", family = "serif"),
    axis.line          = element_line(color = "black", linewidth = 0.4),
    axis.ticks         = element_line(color = "black", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "gray88", linewidth = 0.25),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(t = 4, r = 6, b = 2, l = 4, unit = "mm")
  )

# Boxplots of observed vs predicted structure against the CONAMA limits
plot_conama_colorido <- function(var, titulo_y) {

  df_plot <- df_conama %>%
    select(estagio, origem, valor = all_of(var)) %>%
    filter(!is.na(valor)) %>%
    mutate(
      estagio = factor(estagio, levels = c("Inicial", "Médio", "Avançado")),
      origem  = factor(origem,  levels = c("Real (Campo)", "Predito (Modelo)"))
    )

  linhas <- conama_ref %>% filter(variavel == var)

  cores_linhas <- setNames(cores_conama[linhas$rotulo_estagio], linhas$rotulo_legenda)
  tipos_linhas <- setNames(c("dashed", "dashed"), linhas$rotulo_legenda)

  ggplot(df_plot, aes(x = estagio, y = valor, fill = origem)) +
    geom_boxplot(
      color          = "gray25",
      outlier.shape  = 21,
      outlier.size   = 1.3,
      outlier.stroke = 0.3,
      outlier.colour = "gray25",
      width          = 0.55,
      position       = position_dodge(width = 0.68),
      linewidth      = 0.35
    ) +
    geom_hline(
      data        = linhas,
      aes(yintercept = limite, color = rotulo_legenda, linetype = rotulo_legenda),
      linewidth   = 0.85,
      inherit.aes = FALSE
    ) +
    scale_x_discrete(labels = labels_estagio) +
    scale_fill_manual(
      values = cores_origem_box,
      name   = "Source",
      labels = labels_origem,
      guide  = guide_legend(
        order = 1, nrow = 1,
        override.aes = list(color = "gray25", linewidth = 0.35, linetype = 0)
      )
    ) +
    scale_color_manual(
      values = cores_linhas,
      name   = "CONAMA 04/94",
      guide  = guide_legend(order = 2, nrow = 2, override.aes = list(linewidth = 0.85))
    ) +
    scale_linetype_manual(
      values = tipos_linhas,
      name   = "CONAMA 04/94",
      guide  = guide_legend(
        order = 2, nrow = 2,
        override.aes = list(color = unname(cores_linhas), linewidth = 0.85)
      )
    ) +
    labs(x = "Successional stage", y = titulo_y) +
    tema_conama
}

p_ab  <- plot_conama_colorido("ab_ha",
                              expression(bold("Basal area (m"^2*" ha"^-1*")")))
p_ht  <- plot_conama_colorido("ht_med",  "Mean height (m)")
p_dap <- plot_conama_colorido("dap_med", "Mean DBH (cm)")

print(p_ab)
print(p_ht)
print(p_dap)

# Export individual panels (10x10 cm, 300 dpi, TIFF-LZW)
ggsave("figura_conama_ab_color.tiff",  plot = p_ab,  width = 10, height = 10,
       units = "cm", dpi = 300, compression = "lzw")
ggsave("figura_conama_ht_color.tiff",  plot = p_ht,  width = 10, height = 10,
       units = "cm", dpi = 300, compression = "lzw")
ggsave("figura_conama_dap_color.tiff", plot = p_dap, width = 10, height = 10,
       units = "cm", dpi = 300, compression = "lzw")

# Composite panel (30x10 cm)
p_conama_painel <- grid.arrange(p_ab, p_ht, p_dap, ncol = 3)

ggsave("figura_conama_painel_color.tiff", plot = p_conama_painel,
       width = 30, height = 10, units = "cm", dpi = 300, compression = "lzw")

cat("CONAMA figures exported\n\n")

# Summary table: observed/predicted medians vs legal thresholds
cat("Summary table - medians by stage vs CONAMA 04/94:\n\n")

tabela_resumo <- df_conama %>%
  group_by(estagio, origem) %>%
  summarise(
    ab_ha_mediana   = round(median(ab_ha,   na.rm = TRUE), 2),
    ht_med_mediana  = round(median(ht_med,  na.rm = TRUE), 2),
    dap_med_mediana = round(median(dap_med, na.rm = TRUE), 2),
    n               = n() / 3,
    .groups = "drop"
  ) %>%
  arrange(estagio, origem) %>%
  mutate(
    conama_ab  = case_when(estagio == "Inicial" ~ "<= 8",
                           estagio == "Médio"    ~ "<= 15",
                           estagio == "Avançado" ~ "> 15"),
    conama_ht  = case_when(estagio == "Inicial" ~ "<= 4",
                           estagio == "Médio"    ~ "<= 12",
                           estagio == "Avançado" ~ "> 12"),
    conama_dap = case_when(estagio == "Inicial" ~ "<= 8",
                           estagio == "Médio"    ~ "<= 15",
                           estagio == "Avançado" ~ "> 15")
  )

print(tabela_resumo, n = Inf)

# Compliance: percentage of plots within the legal thresholds (field data)
cat("\n\nCompliance with CONAMA 04/94 (% of plots within thresholds):\n\n")

conformidade <- df_conama %>%
  filter(origem == "Real (Campo)") %>%
  distinct(plot, estagio, ab_ha, ht_med, dap_med) %>%
  mutate(
    conf_ab = case_when(estagio == "Inicial"  ~ ab_ha   <= 8,
                        estagio == "Médio"    ~ ab_ha   <= 15,
                        estagio == "Avançado" ~ ab_ha   >  15),
    conf_ht = case_when(estagio == "Inicial"  ~ ht_med  <= 4,
                        estagio == "Médio"    ~ ht_med  <= 12,
                        estagio == "Avançado" ~ ht_med  >  12),
    conf_dap = case_when(estagio == "Inicial"  ~ dap_med <= 8,
                         estagio == "Médio"    ~ dap_med <= 15,
                         estagio == "Avançado" ~ dap_med >  15),
    conf_total = conf_ab & conf_ht & conf_dap
  ) %>%
  group_by(estagio) %>%
  summarise(
    n         = n(),
    pct_ab    = round(mean(conf_ab,    na.rm = TRUE) * 100, 1),
    pct_ht    = round(mean(conf_ht,    na.rm = TRUE) * 100, 1),
    pct_dap   = round(mean(conf_dap,   na.rm = TRUE) * 100, 1),
    pct_todos = round(mean(conf_total, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  rename(
    "Estágio"        = estagio,
    "n"              = n,
    "% conf. AB"     = pct_ab,
    "% conf. Height" = pct_ht,
    "% conf. DBH"    = pct_dap,
    "% conf. Total"  = pct_todos
  )

print(conformidade)

cat("\nStep 10 complete\n\n")


# =============================================================================
# FINAL SUMMARY
# =============================================================================
cat("\n", rep("=", 70), "\n")
cat("COMPLETE SUMMARY\n")
cat(rep("=", 70), "\n\n")

cat("Variable flow:\n")
cat("  1. Initial:", ncol(df_preprocessado_sem_plot) - 1, "\n")
cat("  2. After VIF:", ncol(df_pos_vif) - 1, "\n")
cat("  3. After Boruta:", ncol(df_pos_boruta) - 1, "\n")
cat("  4. With PCoA:", ncol(df_pos_boruta_pcoa) - 1, "\n\n")

cat("Performance:\n")
cat("  TRAINING - Accuracy:", round(conf_treino$overall["Accuracy"], 4), "\n")
cat("  TEST     - Accuracy:", round(conf_teste$overall["Accuracy"], 4), "\n\n")

cat("Pipeline complete\n")

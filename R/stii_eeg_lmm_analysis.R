# =============================================================================
# STII Study EEG ACL — Complete Statistical Analysis Script
# Probabilistic Go/No-Go Task: Behavioral, Oscillatory, and Brain-Behavior LMMs
# =============================================================================
# Requirements:
#   install.packages(c("readxl", "dplyr", "tidyr", "lme4", "lmerTest",
#                      "emmeans", "effectsize", "ggplot2"))
# =============================================================================

library(readxl)
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(effectsize)

# ── 0. LOAD DATA ──────────────────────────────────────────────────────────────
# Participant-level data are not included in this repository (see README).
# Place the data file in data/ and run this script from the repository root.

df_raw <- read_excel(file.path("data", "ST_EEG_ACL.xlsx"))

# ── 0.1 Recode variables ──────────────────────────────────────────────────────

df_raw <- df_raw %>%
  mutate(
    Group = factor(Group,
                   levels = c(1, 4),
                   labels = c("Non-ACLR", "ACLR")),
    Congruence = case_when(
      grepl("INCONG", Condition) ~ "Incongruent",
      grepl("CONG",   Condition) ~ "Congruent",
      TRUE                       ~ "Total"
    ),
    Probability = case_when(
      grepl("FACIL",   Condition) ~ "High",
      grepl("MEDIO",   Condition) ~ "Medium",
      grepl("DIFICIL", Condition) ~ "Low",
      TRUE                        ~ "Total"
    ),
    Probability = factor(Probability, levels = c("High", "Medium", "Low")),
    Congruence  = factor(Congruence,  levels = c("Congruent", "Incongruent")),
    Subject     = factor(Subject)
  )

cat("Total rows:", nrow(df_raw), "\n")
cat("Subjects:", nlevels(df_raw$Subject), "\n")
cat("Groups:", table(df_raw$Group[!duplicated(df_raw$Subject)]), "\n")


# =============================================================================
# 1. BEHAVIORAL ANALYSIS
# =============================================================================
# Design: DV ~ Group * Probability + (1 | Subject)
# One observation per Subject × Probability (behavioral vars do not vary by Congruence)
# =============================================================================

# Build behavioral dataset: one row per Subject × Probability
beh <- df_raw %>%
  filter(
    Congruence  == "Congruent",
    Probability != "Total"
  ) %>%
  group_by(Subject, Group, Probability) %>%
  summarise(
    RTm      = first(RTm),
    PART     = first(PART),
    IFRT     = first(IFRT),
    Accuracy = first(Accuracy),
    CER      = first(CER),
    PIF      = first(PIF),
    .groups  = "drop"
  )

cat("\nBehavioral dataset:", nrow(beh), "rows,",
    n_distinct(beh$Subject), "subjects\n")

# ── Helper: run one behavioral LMM and extract results ────────────────────────
run_beh_lmm <- function(dv, data = beh) {
  
  formula_str <- paste0(dv, " ~ Group * Probability + (1 | Subject)")
  m <- lmer(as.formula(formula_str), data = data, REML = TRUE)
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("BEHAVIORAL LMM:", dv, "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  # Type III F-tests
  at <- anova(m, type = 3)
  print(at)
  
  # Effect sizes (eta-squared partial)
  cat("\nEffect sizes (eta_p^2):\n")
  es <- eta_squared(m, partial = TRUE, ci = 0.95)
  print(es)
  
  # Marginal means and post-hoc: Group within each Probability level
  # (only if interaction is significant or marginal)
  emm <- emmeans(m, ~ Group | Probability)
  cat("\nPost-hoc: Group within each Probability (Bonferroni):\n")
  ph <- contrast(emm, method = "revpairwise", adjust = "bonferroni")
  ph_df <- as.data.frame(summary(ph, infer = TRUE))
  print(ph_df)
  
  # Cohen's d for each post-hoc contrast
  cat("\nCohen's d (using emmeans + effectsize):\n")
  eff <- eff_size(emm, sigma = sigma(m), edf = df.residual(m))
  print(eff)
  
  invisible(list(model = m, anova = at, eta = es, posthoc = ph_df))
}

beh_results <- list()
for (dv in c("RTm", "PART", "IFRT", "Accuracy", "CER", "PIF")) {
  beh_results[[dv]] <- run_beh_lmm(dv)
}


# =============================================================================
# 2. OSCILLATORY ANALYSIS
# =============================================================================
# Design: Power ~ Group * Probability * Congruence + (1 | Subject)
# Separate models per Marker × Lock
# Restricted to correct trials (S1_CORRECTO)
# =============================================================================

osc <- df_raw %>%
  filter(
    OutcomeLabel == "S1_CORRECTO",
    Congruence   %in% c("Congruent", "Incongruent"),
    Probability  != "Total"
  ) %>%
  group_by(Subject, Group, Marker, Lock, Probability, Congruence) %>%
  summarise(
    Power = mean(Power_dB_average, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nOscillatory dataset:", nrow(osc), "rows\n")
cat("Marker × Lock combinations:\n")
print(osc %>% count(Marker, Lock))

# ── Helper: run one oscillatory LMM ──────────────────────────────────────────
run_osc_lmm <- function(marker, lock, data = osc) {
  
  sub <- data %>%
    filter(Marker == marker, Lock == lock) %>%
    drop_na()
  
  cat("\n", rep("=", 70), "\n", sep = "")
  cat("OSCILLATORY LMM:", marker, "|", lock, "\n")
  cat("n rows:", nrow(sub), "| n subjects:", n_distinct(sub$Subject), "\n")
  cat(rep("=", 70), "\n", sep = "")
  
  m <- lmer(
    Power ~ Group * Probability * Congruence + (1 | Subject),
    data = sub, REML = TRUE
  )
  
  at <- anova(m, type = 3)
  print(at)
  
  cat("\nEffect sizes (eta_p^2):\n")
  es <- eta_squared(m, partial = TRUE, ci = 0.95)
  print(es)
  
  # ── Efecto principal Grupo (overall) ─────────────────────────────────────────
  emm_g <- emmeans(m, ~ Group)
  cat("\nMain effect Group:\n")
  print(summary(emm_g))
  ph_g <- contrast(emm_g, method = "revpairwise")
  print(summary(ph_g, infer = TRUE))
  eff_g <- eff_size(emm_g, sigma = sigma(m), edf = df.residual(m))
  print(eff_g)
  
  # ── Efecto principal Congruencia (Go vs. No-Go) ───────────────────────────────
  emm_c <- emmeans(m, ~ Congruence)
  cat("\nMain effect Congruence (Go vs. No-Go):\n")
  print(summary(emm_c))
  ph_c <- contrast(emm_c, method = "revpairwise")
  print(summary(ph_c, infer = TRUE))
  eff_c <- eff_size(emm_c, sigma = sigma(m), edf = df.residual(m))
  print(eff_c)
  
  # ── Post-hoc interacción Group × Probability (Bonferroni × 3) ────────────────
  emm_gp <- emmeans(m, ~ Group | Probability)
  cat("\nGroup × Probability post-hoc (Bonferroni × 3):\n")
  ph_gp <- contrast(emm_gp, method = "revpairwise", adjust = "bonferroni")
  print(summary(ph_gp, infer = TRUE))
  eff_gp <- eff_size(emm_gp, sigma = sigma(m), edf = df.residual(m))
  print(eff_gp)
  
  # ── Post-hoc interacción triple Group × Probability × Congruence (Bonferroni × 6) ──
  emm_3 <- emmeans(m, ~ Group | Probability | Congruence)
  cat("\nGroup × Probability × Congruence post-hoc (Bonferroni × 6):\n")
  ph_3 <- contrast(emm_3, method = "revpairwise", adjust = "bonferroni")
  print(summary(ph_3, infer = TRUE))
  eff_3 <- eff_size(emm_3, sigma = sigma(m), edf = df.residual(m))
  print(eff_3)
  
  invisible(list(model = m, anova = at, eta = es))
}

# Run all 6 Marker × Lock combinations
osc_results <- list()
for (marker in c("MidfrontalTheta", "LeftCentralBeta", "PosteriorAlpha")) {
  for (lock in c("stim", "resp")) {
    key <- paste(marker, lock, sep = "_")
    osc_results[[key]] <- run_osc_lmm(marker, lock)
  }
}


# =============================================================================
# 3. BRAIN-BEHAVIOR RELATIONSHIPS
# =============================================================================
# For each Marker × Lock: Pearson correlations between subject-level mean power
# and each behavioral variable, overall and within each group.
# Group difference in r tested via Fisher's z-transformation.
# =============================================================================

# Subject-level mean power (across all conditions)
pow_subj <- osc %>%
  group_by(Subject, Group, Marker, Lock) %>%
  summarise(Power = mean(Power, na.rm = TRUE), .groups = "drop")

# Subject-level mean behavioral (across all probability levels)
beh_subj <- beh %>%
  group_by(Subject, Group) %>%
  summarise(across(c(RTm, PART, IFRT, Accuracy, CER, PIF), mean, na.rm = TRUE),
            .groups = "drop")

bb_data <- pow_subj %>%
  left_join(beh_subj, by = c("Subject", "Group"))

cat("\n", rep("=", 70), "\n", sep = "")
cat("BRAIN-BEHAVIOR CORRELATIONS\n")
cat(rep("=", 70), "\n\n", sep = "")

fisher_z_test <- function(r1, r2, n1, n2) {
  z1   <- atanh(r1)
  z2   <- atanh(r2)
  se   <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  z    <- (z1 - z2) / se
  p    <- 2 * pnorm(-abs(z))
  list(z = z, p = p)
}

r_ci <- function(r, n, conf = 0.95) {
  z    <- atanh(r)
  se   <- 1 / sqrt(n - 3)
  zc   <- qnorm(1 - (1 - conf) / 2)
  c(lo = tanh(z - zc * se), hi = tanh(z + zc * se))
}

beh_vars <- c("RTm", "PART", "IFRT", "Accuracy", "CER", "PIF")

bb_results <- data.frame()

for (marker in c("MidfrontalTheta", "LeftCentralBeta", "PosteriorAlpha")) {
  for (lock in c("stim", "resp")) {
    
    sub <- bb_data %>% filter(Marker == marker, Lock == lock)
    g0  <- sub %>% filter(Group == "Non-ACLR")
    g1  <- sub %>% filter(Group == "ACLR")
    n_all <- nrow(sub); n0 <- nrow(g0); n1 <- nrow(g1)
    
    cat("\n---", marker, "|", lock, "---\n")
    cat(sprintf("  n = %d (Non-ACLR: %d, ACLR: %d)\n", n_all, n0, n1))
    
    for (bv in beh_vars) {
      
      ct_all <- cor.test(sub$Power, sub[[bv]])
      r_all  <- ct_all$estimate
      p_all  <- ct_all$p.value
      ci_all <- r_ci(r_all, n_all)
      
      ct0 <- cor.test(g0$Power, g0[[bv]])
      r0  <- ct0$estimate; p0 <- ct0$p.value; ci0 <- r_ci(r0, n0)
      
      ct1 <- cor.test(g1$Power, g1[[bv]])
      r1  <- ct1$estimate; p1 <- ct1$p.value; ci1 <- r_ci(r1, n1)
      
      fz  <- fisher_z_test(r0, r1, n0, n1)
      
      cat(sprintf(
        "  %-9s  Overall: r=%+.3f [%+.3f,%+.3f] P=%s  |  Non-ACLR: r=%+.3f P=%s  |  ACLR: r=%+.3f P=%s  |  GroupDiff P=%s\n",
        bv,
        r_all, ci_all["lo"], ci_all["hi"], format.pval(p_all, digits = 3, eps = 0.001),
        r0, format.pval(p0, digits = 3, eps = 0.001),
        r1, format.pval(p1, digits = 3, eps = 0.001),
        format.pval(fz$p, digits = 3, eps = 0.001)
      ))
      
      bb_results <- rbind(bb_results, data.frame(
        Marker        = marker,
        Lock          = lock,
        BehVar        = bv,
        r_overall     = round(r_all, 3),
        CI_lo         = round(unname(ci_all["lo"]), 3),
        CI_hi         = round(unname(ci_all["hi"]), 3),
        P_overall     = round(p_all, 4),
        r_NonACLR     = round(r0, 3),
        P_NonACLR     = round(p0, 4),
        r_ACLR        = round(r1, 3),
        P_ACLR        = round(p1, 4),
        FisherZ       = round(fz$z, 3),
        P_groupdiff   = round(fz$p, 4),
        stringsAsFactors = FALSE
      ))
    }
  }
}

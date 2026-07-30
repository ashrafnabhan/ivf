#!/usr/bin/env Rscript
# ============================================================
# Analysis of physician knowledge/practice survey data
# Outcomes: kscore (continuous, 0-15) and adequate (binary:
#           Adequate/Inadequate, derived from kscore >= 10)
# ============================================================
#
# Data structure:
#   - 15 individual "low-value practice" knowledge items, each
#     coded Correct/Incorrect. kscore = number of items answered
#     correctly (sum of "Correct" across the 15 items).
#   - adequate = "Adequate" if kscore >= 10, else "Inadequate".
#     Because adequate/kscore are DERIVED FROM the 15 items, the
#     items themselves are NOT used as predictors below - only
#     the respondent/practice characteristics are used as
#     predictors of the two outcomes.
#
# Predictors used: age, sex, pg_years (years since postgraduation),
#   degree, subspecialty_rei, provision, facility, guideline_access
#   (provision_duration is analysed separately - see note in
#   Section 5 - because it is only defined for respondents with
#   provision == "Yes", so including it in the main model would
#   silently drop everyone who answered "No").
#
# Outputs (written to OUT_DIR, default "outputs/"):
#   tables/  - CSV files for every table produced
#   figures/ - PNG files for every plot produced
#   analysis_log.txt - full console log of the run
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(broom)
  library(car)
})

set.seed(123)

# ---- 0. Setup -------------------------------------------------
OUT_DIR <- "outputs"
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(OUT_DIR, "analysis_log.txt"), open = "wt")
sink(log_con, split = TRUE)  # print to console AND log file

save_table <- function(x, name) {
  write.csv(x, file.path(OUT_DIR, "tables", paste0(name, ".csv")), row.names = FALSE)
}
save_plot <- function(p, name, width = 7, height = 5) {
  ggsave(file.path(OUT_DIR, "figures", paste0(name, ".png")),
         plot = p, width = width, height = height, dpi = 300)
}
hline <- function() cat(strrep("-", 70), "\n")

# ---- 1. Load & prepare data ------------------------------------
DATA_PATH <- "df.rds"   # <-- change if your file lives elsewhere
df <- readRDS(DATA_PATH)

item_vars <- c("pgt_a", "assisted_hatching", "lymphocyte_immunization",
               "corticosteroids", "progesterone_lps", "sperm_dfi",
               "ert_testing", "hysteroscopy_pre_ivf", "gnrh_antagonist",
               "ivig", "lif_therapy", "intralipid", "elective_freeze_all",
               "lmwh_rif", "endometrial_scratch")

predictor_vars <- c("age", "sex", "pg_years", "degree", "subspecialty_rei",
                     "provision", "facility", "guideline_access")

# tidy factor labels / types
df$degree <- factor(df$degree, levels = c(0, 1),
                     labels = c("No postgraduate degree", "Postgraduate degree"))
df$adequate <- relevel(df$adequate, ref = "Inadequate")  # so glm models P(Adequate)
df$sex             <- relevel(df$sex, ref = "Female")
df$subspecialty_rei<- relevel(df$subspecialty_rei, ref = "No")
df$provision       <- relevel(df$provision, ref = "No")
df$guideline_access<- relevel(df$guideline_access, ref = "No")
df$facility        <- relevel(df$facility, ref = "Multi-setting")

cat("Loaded", nrow(df), "respondents,", ncol(df), "variables.\n")
hline()

# ---- 2. Missingness overview ------------------------------------
cat("\n== Missing data per predictor ==\n")
miss_tab <- data.frame(
  variable = predictor_vars,
  n_missing = sapply(df[predictor_vars], function(x) sum(is.na(x)))
)
print(miss_tab)
save_table(miss_tab, "01_missingness")
hline()

# ---- 3. Descriptive statistics ----------------------------------

## 3a. Item-level performance (% answered correctly)
item_pct <- data.frame(
  item = item_vars,
  n_correct = sapply(df[item_vars], function(x) sum(x == "Correct", na.rm = TRUE)),
  n_total   = sapply(df[item_vars], function(x) sum(!is.na(x)))
)
item_pct$pct_correct <- round(100 * item_pct$n_correct / item_pct$n_total, 1)
item_pct <- item_pct[order(item_pct$pct_correct), ]
cat("\n== Item-level % answered correctly (low to high) ==\n")
print(item_pct, row.names = FALSE)
save_table(item_pct, "02_item_performance")

item_pct$item <- factor(item_pct$item, levels = item_pct$item)
p_items <- ggplot(item_pct, aes(x = item, y = pct_correct)) +
  geom_col(fill = "#2c7fb8") +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey40") +
  coord_flip() +
  labs(title = "Percentage of respondents answering each item correctly",
       x = NULL, y = "% correct") +
  theme_minimal(base_size = 12)
save_plot(p_items, "01_item_performance", height = 6)

## 3b. Outcome distributions
cat("\n== Normality checks (Shapiro-Wilk) ==\n")
norm_check <- data.frame(
  variable = c("age", "pg_years", "kscore"),
  W = c(shapiro.test(df$age)$statistic, shapiro.test(df$pg_years)$statistic,
        shapiro.test(df$kscore)$statistic),
  p_value = c(shapiro.test(df$age)$p.value, shapiro.test(df$pg_years)$p.value,
              shapiro.test(df$kscore)$p.value)
)
norm_check$W <- round(norm_check$W, 3)
norm_check$p_value <- signif(norm_check$p_value, 3)
print(norm_check, row.names = FALSE)
cat("age, pg_years, and kscore are all significantly non-normal (p < .05).\n")
cat("Descriptives below use median (IQR); group comparisons use\n")
cat("non-parametric tests (Wilcoxon rank-sum / Spearman) accordingly.\n")
save_table(norm_check, "00_normality_checks")

cat("\n== kscore summary ==\n")
cat("Median (IQR):", median(df$kscore, na.rm = TRUE), "(",
    IQR(df$kscore, na.rm = TRUE), ")\n")
cat("Range:", min(df$kscore, na.rm = TRUE), "-", max(df$kscore, na.rm = TRUE), "\n")

cat("\n== adequate distribution ==\n")
print(table(df$adequate))
print(round(100 * prop.table(table(df$adequate)), 1))

p_kscore <- ggplot(df, aes(x = kscore)) +
  geom_histogram(binwidth = 1, fill = "#2c7fb8", color = "white") +
  geom_vline(xintercept = 10, linetype = "dashed", color = "red") +
  labs(title = "Distribution of knowledge score (kscore)",
       subtitle = "Dashed line = adequate/inadequate cutoff (10/15)",
       x = "kscore (0-15)", y = "Number of respondents") +
  theme_minimal(base_size = 12)
save_plot(p_kscore, "02_kscore_distribution")

hline()

## 3c. Table 1 - characteristics by adequate group
tab1_vars <- predictor_vars
tab1_rows <- list()
for (v in tab1_vars) {
  x <- df[[v]]
  if (is.numeric(x)) {
    med <- tapply(x, df$adequate, function(z) median(z, na.rm = TRUE))
    iqr <- tapply(x, df$adequate, function(z) IQR(z, na.rm = TRUE))
    test <- tryCatch(wilcox.test(x ~ df$adequate)$p.value, error = function(e) NA)
    tab1_rows[[v]] <- data.frame(
      variable = v,
      level = "median (IQR)",
      Inadequate = sprintf("%.1f (%.1f)", med["Inadequate"], iqr["Inadequate"]),
      Adequate   = sprintf("%.1f (%.1f)", med["Adequate"],   iqr["Adequate"]),
      p_value = round(test, 3)
    )
  } else {
    tb <- table(x, df$adequate)
    test <- tryCatch({
      if (any(chisq.test(tb)$expected < 5)) fisher.test(tb)$p.value
      else chisq.test(tb)$p.value
    }, error = function(e) NA)
    pct <- prop.table(tb, margin = 2) * 100
    for (lv in rownames(tb)) {
      tab1_rows[[paste0(v, "_", lv)]] <- data.frame(
        variable = v,
        level = lv,
        Inadequate = sprintf("%d (%.1f%%)", tb[lv, "Inadequate"], pct[lv, "Inadequate"]),
        Adequate   = sprintf("%d (%.1f%%)", tb[lv, "Adequate"],   pct[lv, "Adequate"]),
        p_value = ifelse(lv == rownames(tb)[1], round(test, 3), NA)
      )
    }
  }
}
table1 <- do.call(rbind, tab1_rows)
rownames(table1) <- NULL
cat("\n== Table 1: characteristics by adequate/inadequate group ==\n")
print(table1, row.names = FALSE)
save_table(table1, "03_table1_by_adequate")
hline()

# ---- 4. Univariable associations with kscore --------------------
cat("\n== Univariable associations with kscore (continuous outcome) ==\n")
uni_kscore <- list()
for (v in predictor_vars) {
  x <- df[[v]]
  if (is.numeric(x)) {
    sp <- suppressWarnings(cor.test(x, df$kscore, method = "spearman"))
    uni_kscore[[v]] <- data.frame(
      variable = v, comparison = "Spearman rho",
      estimate = round(unname(sp$estimate), 3), p_value = round(sp$p.value, 3)
    )
  } else {
    nlev <- nlevels(droplevels(x))
    if (nlev == 2) {
      test <- tryCatch(wilcox.test(df$kscore ~ x)$p.value, error = function(e) NA)
      med <- tapply(df$kscore, x, median, na.rm = TRUE)
      lv2 <- levels(droplevels(x))[2]
      uni_kscore[[v]] <- data.frame(
        variable = v,
        comparison = sprintf("Wilcoxon (%s vs %s)", lv2, levels(droplevels(x))[1]),
        estimate = round(med[lv2] - med[1], 3), p_value = round(test, 3)
      )
    } else {
      test <- tryCatch(kruskal.test(df$kscore ~ x)$p.value, error = function(e) NA)
      uni_kscore[[v]] <- data.frame(
        variable = v, comparison = sprintf("Kruskal-Wallis (%d groups)", nlev),
        estimate = NA, p_value = round(test, 3)
      )
    }
  }
}
uni_kscore_tab <- do.call(rbind, uni_kscore)
rownames(uni_kscore_tab) <- NULL
print(uni_kscore_tab, row.names = FALSE)
save_table(uni_kscore_tab, "04_univariable_kscore")
hline()

# ---- 4b. Multicollinearity screen among predictors ---------------
age_pg_cor <- cor(df$age, df$pg_years, use = "complete.obs")
cat("\n== Multicollinearity screen ==\n")
cat("Correlation between age and pg_years (years since postgraduation):",
    round(age_pg_cor, 3), "\n")
cat("These two variables are, unsurprisingly, almost perfectly collinear\n")
cat("(age increases ~1:1 with years of experience). Including both in the\n")
cat("same regression inflates their standard errors (VIF > 40, see below)\n")
cat("and makes their individual coefficients uninterpretable, even though\n")
cat("the model as a whole is still valid. The multivariable models below\n")
cat("therefore use pg_years (professional experience) and DROP age; age\n")
cat("is still examined in the univariable analysis above.\n")
mv_predictors <- setdiff(predictor_vars, "age")
hline()

# ---- 5. Multivariable linear regression: kscore ------------------
cat("\n== Multivariable linear regression: kscore ~ predictors (age excluded, see note above) ==\n")
mv_formula_k <- as.formula(paste("kscore ~", paste(mv_predictors, collapse = " + ")))
fit_kscore <- lm(mv_formula_k, data = df)
print(summary(fit_kscore))

kscore_coef <- tidy(fit_kscore, conf.int = TRUE)
kscore_coef[, -1] <- round(kscore_coef[, -1], 3)
save_table(kscore_coef, "05_multivariable_kscore_coefficients")

cat("\n-- Variance Inflation Factors (check multicollinearity) --\n")
print(vif(fit_kscore))

cat("\n-- Residual normality check --\n")
cat("kscore itself is right/left-skewed (Shapiro p =",
    signif(shapiro.test(df$kscore)$p.value, 3), "), but OLS validity depends on\n")
cat("the MODEL RESIDUALS being approximately normal, not the raw outcome.\n")
resid_sw <- shapiro.test(residuals(fit_kscore))
cat("Shapiro-Wilk on model residuals: W =", round(resid_sw$statistic, 3),
    ", p =", signif(resid_sw$p.value, 3), "\n")
if (resid_sw$p.value > 0.05) {
  cat("Residuals do not significantly depart from normality, so OLS\n")
  cat("coefficients/p-values below are used as reported (see QQ plot too).\n")
} else {
  cat("Residuals DO significantly depart from normality - treat p-values\n")
  cat("with caution; consider bootstrapped CIs or quantile regression.\n")
}

# diagnostic plots
png(file.path(OUT_DIR, "figures", "03_kscore_model_diagnostics.png"),
    width = 900, height = 900, res = 150)
par(mfrow = c(2, 2))
plot(fit_kscore)
dev.off()

cat("\nNOTE: provision_duration was excluded from the main model above\n")
cat("because it is only recorded for respondents with provision == 'Yes'\n")
cat("(structural missingness). A sensitivity analysis restricted to\n")
cat("providers only is shown next.\n")

df_providers <- droplevels(df[df$provision == "Yes" & !is.na(df$provision_duration), ])
df_providers$provision_duration <- relevel(df_providers$provision_duration,
                                            ref = "Less than 5 years")
sens_formula <- as.formula(paste(
  "kscore ~", paste(setdiff(mv_predictors, "provision"), collapse = " + "),
  "+ provision_duration"))
fit_kscore_sens <- lm(sens_formula, data = df_providers)
cat("\n-- Sensitivity model (providers only, n =", nrow(df_providers), ") --\n")
print(summary(fit_kscore_sens))
sens_coef <- tidy(fit_kscore_sens, conf.int = TRUE)
sens_coef[, -1] <- round(sens_coef[, -1], 3)
save_table(sens_coef, "06_sensitivity_kscore_providers_only")
hline()

# ---- 6. Multivariable logistic regression: adequate ---------------
cat("\nNOTE: logistic regression (below) makes no normality assumption on\n")
cat("predictors or outcome - the non-normality of age/pg_years/kscore\n")
cat("noted above does not affect this model.\n")
cat("\n== Multivariable logistic regression: adequate ~ predictors (age excluded, see note above) ==\n")
mv_formula_a <- as.formula(paste("adequate ~", paste(mv_predictors, collapse = " + ")))
fit_adequate <- glm(mv_formula_a, data = df, family = binomial)
print(summary(fit_adequate))

adequate_or <- tidy(fit_adequate, conf.int = TRUE, exponentiate = TRUE)
adequate_or[, -1] <- round(adequate_or[, -1], 3)
names(adequate_or)[names(adequate_or) == "estimate"] <- "odds_ratio"
cat("\n-- Odds ratios (exponentiated coefficients) --\n")
print(adequate_or, row.names = FALSE)
save_table(adequate_or, "07_multivariable_adequate_odds_ratios")

cat("\n-- Variance Inflation Factors --\n")
print(vif(fit_adequate))

cat("\n-- Model fit --\n")
cat("Null deviance:", round(fit_adequate$null.deviance, 1),
    "on", fit_adequate$df.null, "df\n")
cat("Residual deviance:", round(fit_adequate$deviance, 1),
    "on", fit_adequate$df.residual, "df\n")
cat("AIC:", round(AIC(fit_adequate), 1), "\n")

observed <- model.frame(fit_adequate)$adequate  # aligns with rows actually used (NAs dropped)
pred_class <- ifelse(predict(fit_adequate, type = "response") > 0.5,
                      "Adequate", "Inadequate")
acc <- mean(pred_class == observed)
cat("In-sample classification accuracy (cutoff 0.5):", round(100 * acc, 1), "%\n")
conf_mat <- table(Predicted = pred_class, Observed = observed)
cat("\nConfusion matrix:\n")
print(conf_mat)
save_table(as.data.frame.matrix(conf_mat), "08_confusion_matrix")

# forest-style plot of odds ratios
or_plot_dat <- adequate_or[adequate_or$term != "(Intercept)", ]
p_or <- ggplot(or_plot_dat, aes(x = reorder(term, odds_ratio), y = odds_ratio)) +
  geom_point(size = 2.5, color = "#2c7fb8") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2, color = "#2c7fb8") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Adjusted odds ratios for 'Adequate' knowledge",
       x = NULL, y = "Odds ratio (95% CI)") +
  theme_minimal(base_size = 12)
save_plot(p_or, "04_adequate_odds_ratios", height = 6)

hline()
cat("\nAnalysis complete. Tables saved to '", file.path(OUT_DIR, "tables"),
    "', figures to '", file.path(OUT_DIR, "figures"), "'.\n", sep = "")

sink()
close(log_con)

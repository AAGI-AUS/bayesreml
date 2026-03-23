## ============================================================================
##  source("05_compare_barrero.R")
##
##  Compare variance-component effect sizes between greta and ASReml for
##  the barrero.maize MET model:
##
##    fixed:    yield ~ loc + yearf + loc:yearf
##    random:   ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf
##    residual: ~ dsum(~ units|env)   (heterogeneous error by environment)
##
##  All comparisons are on the standard-deviation (sigma) scale.
##  ASReml reports variances (sigma^2) -> sqrt() applied.
##  greta reports sigma directly (posterior mean +/- posterior SD).
##
##  MCMC tuning rationale
##  ---------------------
##  This is a very large model: ~847 gen + 428 rep:loc:yearf + 3426 gen:loc:yearf
##  random-effect levels. Baseline run (500 samples, 1000 warmup) gave
##  Rhat >> 1.1 for all VC parameters (worst: 56 for sigma_yearf_in_loc_in_rep).
##  Two improvements applied here:
##    1. Longer warmup (WARM = 3000) so HMC step-size adaptation finishes.
##    2. Tighter lognormal prior (PRIOR_VC_SD = 0.5 vs default 1.0) to
##       regularise the large sigma space and reduce prior-pull toward exp(0.5)~1.65.
##       lognormal(0, 0.5): median = 1, mean = 1.13, P5/P95 = [0.44, 2.30].
## ============================================================================

suppressPackageStartupMessages({
  library(agridat)
  library(data.table)
  library(asreml)
  library(greta)
  library(coda)
  library(ggplot2)
  library(patchwork)
})
source("02_translator.R")

## ── 0. Settings ───────────────────────────────────────────────────────────────
SMPL         <- 2000L   # posterior samples per chain
WARM         <- 3000L   # warmup (HMC step-size adaptation)
CH           <- 2L      # chains (run in parallel)
PRIOR_VC_SD  <- 0.5     # lognormal(0, PRIOR_VC_SD) prior on sigma parameters
                        # smaller -> tighter prior -> less pull toward median=1
set.seed(2024)

## ── 1. Data ───────────────────────────────────────────────────────────────────
data(barrero.maize)
dat <- data.table(barrero.maize)
dat[, yearf := factor(year)]
dat[, env   := factor(env)]
dat[, loc   := factor(loc)]
dat[, gen   := factor(gen)]
dat[, rep   := factor(rep)]
dat <- dat[order(env)][!is.na(yield)]

env_levels <- levels(dat$env)
n_env      <- nlevels(dat$env)

cat(sprintf("Data: %d rows | %d gen | %d env | %d loc | %d year\n",
            nrow(dat), nlevels(dat$gen), n_env,
            nlevels(dat$loc), nlevels(dat$yearf)))
cat(sprintf("Random-effect levels: %d rep:loc:yearf | %d gen:loc:yearf\n",
            length(unique(paste(dat$rep, dat$loc, dat$yearf))),
            length(unique(paste(dat$gen, dat$loc, dat$yearf)))))

## ── 2. ASReml fit ─────────────────────────────────────────────────────────────
cat("\n--- Fitting ASReml ---\n")
m_asreml <- asreml(
  yield    ~ loc + yearf + loc:yearf,
  random   = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
  residual = ~ dsum(~ units|env),
  data     = dat,
  workspace = "500mb",
  trace    = FALSE
)
cat("ASReml converged:", m_asreml$converge, "\n")
vc_asreml_raw <- summary(m_asreml)$varcomp

## ── 3. greta fit ──────────────────────────────────────────────────────────────
cat(sprintf(
  "\n--- Fitting greta (prior_vc_sd=%.1f | warmup=%d | samples=%d | chains=%d) ---\n",
  PRIOR_VC_SD, WARM, SMPL, CH))

m_greta <- asreml_to_greta(
  fixed         = yield ~ loc + yearf + loc:yearf,
  random        = ~ gen + rep:loc:yearf + gen:yearf + gen:loc + gen:loc:yearf,
  rcov          = ~ dsum(~ units|env),
  data          = dat,
  n_samples     = SMPL,
  warmup        = WARM,
  chains        = CH,
  prior_vc_sd   = PRIOR_VC_SD,
  mcmc_verbose  = TRUE,
  verbose       = FALSE
)
cat("greta run time:", round(m_greta$run_time, 1), "s\n")

## ── 4. Convergence diagnostics ────────────────────────────────────────────────
vc_params <- c("sigma_gen", "sigma_yearf_in_loc_in_rep",
               "sigma_yearf_in_gen", "sigma_loc_in_gen",
               "sigma_yearf_in_loc_in_gen")

st <- m_greta$summary$statistics

diag_df <- do.call(rbind, lapply(vc_params, function(p) {
  neff <- effectiveSize(m_greta$draws[, p, drop = FALSE])
  rhat <- gelman.diag(m_greta$draws[, p, drop = FALSE],
                      autoburnin = FALSE)$psrf[1, "Point est."]
  data.frame(param = p,
             mean  = st[p, "Mean"],
             sd    = st[p, "SD"],
             n_eff = round(neff),
             Rhat  = round(rhat, 3),
             stringsAsFactors = FALSE)
}))

cat("\n")
cat(strrep("=", 72), "\n")
cat("  MCMC CONVERGENCE DIAGNOSTICS\n")
cat(strrep("=", 72), "\n")
cat(sprintf("  %-34s %8s %8s %8s %8s\n",
            "param", "mean", "sd", "n_eff", "Rhat"))
cat(strrep("-", 72), "\n")
for (i in seq_len(nrow(diag_df))) {
  flag <- if (diag_df$Rhat[i] > 1.1) " <<" else ""
  cat(sprintf("  %-34s %8.4f %8.4f %8.0f %8.3f%s\n",
              diag_df$param[i], diag_df$mean[i], diag_df$sd[i],
              diag_df$n_eff[i], diag_df$Rhat[i], flag))
}
cat(strrep("=", 72), "\n")
conv_ok <- all(diag_df$Rhat < 1.1)
cat(sprintf("  Convergence (Rhat < 1.1 for all VC): %s\n",
            if (conv_ok) "PASS" else "FAIL -- interpret estimates with caution"))

## ── 5. Extract and align variance components ──────────────────────────────────
vc_map <- data.frame(
  asreml_term = c("gen", "rep:loc:yearf", "gen:yearf",
                  "gen:loc", "gen:loc:yearf"),
  greta_sigma = c("sigma_gen", "sigma_yearf_in_loc_in_rep",
                  "sigma_yearf_in_gen", "sigma_loc_in_gen",
                  "sigma_yearf_in_loc_in_gen"),
  label       = c("gen", "rep:loc:yearf", "gen:yearf",
                  "gen:loc", "gen:loc:yearf"),
  stringsAsFactors = FALSE
)

## ASReml: sigma^2 -> sigma via sqrt()
vc_random <- data.frame(
  label        = vc_map$label,
  asreml_sigma = sqrt(vc_asreml_raw[vc_map$asreml_term, "component"]),
  asreml_se    = 0.5 * vc_asreml_raw[vc_map$asreml_term, "std.error"] /
                       sqrt(vc_asreml_raw[vc_map$asreml_term, "component"]),
  greta_mean   = st[vc_map$greta_sigma, "Mean"],
  greta_sd     = st[vc_map$greta_sigma, "SD"],
  Rhat         = diag_df$Rhat,
  stringsAsFactors = FALSE
)
vc_random$ratio <- vc_random$greta_mean / vc_random$asreml_sigma

## Residual VCs by environment
asreml_resid_rows  <- grep("!R$", rownames(vc_asreml_raw), value = TRUE)
asreml_resid_env   <- sub("^env_", "", sub("!R$", "", asreml_resid_rows))
asreml_resid_sigma <- sqrt(vc_asreml_raw[asreml_resid_rows, "component"])
asreml_resid_se    <- 0.5 * vc_asreml_raw[asreml_resid_rows, "std.error"] /
                            asreml_resid_sigma

greta_resid_rows <- grep("^sigma_e_atg\\[", rownames(st), value = TRUE)
k_idx            <- as.integer(sub("^sigma_e_atg\\[(\\d+),1\\]$", "\\1",
                                   greta_resid_rows))
greta_resid_env  <- env_levels[k_idx]
greta_resid_mean <- st[greta_resid_rows, "Mean"]
greta_resid_sd   <- st[greta_resid_rows, "SD"]

resid_df <- merge(
  data.frame(env          = asreml_resid_env,
             asreml_sigma = asreml_resid_sigma,
             asreml_se    = asreml_resid_se,
             stringsAsFactors = FALSE),
  data.frame(env        = greta_resid_env,
             greta_mean = greta_resid_mean,
             greta_sd   = greta_resid_sd,
             stringsAsFactors = FALSE),
  by = "env", all = TRUE
)
resid_df <- resid_df[order(resid_df$asreml_sigma), ]
resid_df$env <- factor(resid_df$env, levels = resid_df$env)

## ── 6. Print comparison table ─────────────────────────────────────────────────
cat("\n")
cat(strrep("=", 78), "\n")
cat("  RANDOM EFFECT VARIANCE COMPONENTS (sigma scale) | Rhat flagged if > 1.1\n")
cat(strrep("=", 78), "\n")
cat(sprintf("  %-22s  %9s %8s  %9s %8s  %6s  %5s\n",
            "Term", "ASReml-sd", "(SE)", "greta-sd", "(SD)", "ratio", "Rhat"))
cat(strrep("-", 78), "\n")
for (i in seq_len(nrow(vc_random))) {
  flag <- if (vc_random$Rhat[i] > 1.1) " *" else ""
  cat(sprintf("  %-22s  %9.4f %8.4f  %9.4f %8.4f  %6.2f  %5.3f%s\n",
              vc_random$label[i],
              vc_random$asreml_sigma[i], vc_random$asreml_se[i],
              vc_random$greta_mean[i],   vc_random$greta_sd[i],
              vc_random$ratio[i],        vc_random$Rhat[i], flag))
}
cat("  * Rhat > 1.1 -> not converged\n")
cat(strrep("=", 78), "\n")

r_resid  <- cor(resid_df$asreml_sigma, resid_df$greta_mean, use = "complete.obs")
mad_resid <- mean(abs(resid_df$asreml_sigma - resid_df$greta_mean), na.rm = TRUE)
cat(sprintf("\nResidual sd across %d environments:\n", nrow(resid_df)))
cat(sprintf("  Pearson r (ASReml vs greta):  %.4f\n", r_resid))
cat(sprintf("  Mean absolute difference:     %.4f\n", mad_resid))

## ── 7. Plots ──────────────────────────────────────────────────────────────────

## 7a. Random effects — dot plot
vc_long <- rbind(
  data.frame(label   = vc_random$label, package = "ASReml",
             sigma   = vc_random$asreml_sigma, se = vc_random$asreml_se,
             stringsAsFactors = FALSE),
  data.frame(label   = vc_random$label, package = "greta",
             sigma   = vc_random$greta_mean,   se = vc_random$greta_sd,
             stringsAsFactors = FALSE)
)
vc_long$label   <- factor(vc_long$label, levels = rev(vc_map$label))
vc_long$package <- factor(vc_long$package, levels = c("ASReml", "greta"))

## Annotate convergence
conv_ann <- vc_random[, c("label", "Rhat")]
conv_ann$conv_label <- sprintf("Rhat=%.2f%s", conv_ann$Rhat,
                                ifelse(conv_ann$Rhat > 1.1, "*", ""))
conv_ann$x_pos <- max(c(vc_long$sigma + vc_long$se), na.rm = TRUE) * 1.05
conv_ann$label <- factor(conv_ann$label, levels = rev(vc_map$label))

p_vc <- ggplot(vc_long, aes(x = sigma, y = label, colour = package,
                             xmin = sigma - se, xmax = sigma + se)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.6) +
  geom_text(data = conv_ann,
            aes(x = x_pos, y = label, label = conv_label),
            inherit.aes = FALSE, size = 2.8, hjust = 0,
            colour = ifelse(conv_ann$Rhat > 1.1, "red", "grey40")) +
  scale_colour_manual(values = c(ASReml = "#1f78b4", greta = "#e31a1c")) +
  labs(title    = sprintf("Random effect sigma -- ASReml vs greta  (prior_vc_sd=%.1f, warmup=%d)",
                          PRIOR_VC_SD, WARM),
       subtitle = "+/- 1 SE (ASReml) or +/- 1 posterior SD (greta). * = Rhat > 1.1 (not converged)",
       x = "Standard deviation (sigma)", y = NULL, colour = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom") +
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5, 80, 5, 5))

## 7b. Residual sigma scatter
p_resid_scatter <- ggplot(resid_df, aes(x = asreml_sigma, y = greta_mean)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(ymin = greta_mean - greta_sd,
                    ymax = greta_mean + greta_sd),
                colour = "#e31a1c", width = 0, alpha = 0.35) +
  geom_linerange(aes(y = greta_mean,
                     xmin = asreml_sigma - asreml_se,
                     xmax = asreml_sigma + asreml_se),
                 colour = "#1f78b4", alpha = 0.35) +
  geom_point(size = 1.8, colour = "#555555") +
  annotate("text", x = -Inf, y = Inf,
           label = sprintf("r = %.3f\nMAD = %.3f", r_resid, mad_resid),
           hjust = -0.1, vjust = 1.2, size = 3.5) +
  labs(title    = "Residual sigma by environment (107 envs)",
       subtitle = "Dashed line = 1:1 agreement",
       x = "ASReml sigma  (sqrt variance component)",
       y = "greta posterior mean sigma") +
  theme_bw(base_size = 11)

## 7c. Forest plot of residual sigma by environment
resid_long <- rbind(
  data.frame(env = resid_df$env, package = "ASReml",
             sigma = resid_df$asreml_sigma, se = resid_df$asreml_se,
             stringsAsFactors = FALSE),
  data.frame(env = resid_df$env, package = "greta",
             sigma = resid_df$greta_mean,   se = resid_df$greta_sd,
             stringsAsFactors = FALSE)
)
resid_long$package <- factor(resid_long$package, levels = c("ASReml", "greta"))

p_resid_forest <- ggplot(resid_long,
                          aes(x = sigma, y = env, colour = package,
                              xmin = sigma - se, xmax = sigma + se)) +
  geom_pointrange(position = position_dodge(width = 0.5),
                  size = 0.2, linewidth = 0.3) +
  scale_colour_manual(values = c(ASReml = "#1f78b4", greta = "#e31a1c")) +
  labs(title = "Residual sigma by environment (ordered by ASReml estimate)",
       x = "sigma", y = "Environment", colour = NULL) +
  theme_bw(base_size = 8) +
  theme(axis.text.y  = element_text(size = 4.5),
        legend.position = "bottom")

## Save
pdf_name <- sprintf("compare_barrero_greta_asreml_pvc%.1f_w%d.pdf",
                    PRIOR_VC_SD, WARM)
pdf(pdf_name, width = 10, height = 14)
print(
  p_vc /
  p_resid_scatter /
  p_resid_forest +
  plot_layout(heights = c(1.3, 1.5, 3))
)
invisible(dev.off())
cat(sprintf("\nPlots saved to %s\n", pdf_name))

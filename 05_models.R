# 05_models.R ----------------------------------------------------------------------
# Nesting occurrence models for the five species.
#   Study area : segments in the species' RMUs (buffered 50 km) whose mean SST lies within the min-max range of
#                the species' nesting segments.
#   Predictors : seaward moonlight, onshore wave exposure (log(x+1)), onshore wind exposure (log(x+1)),
#                proximity to a strong current (-log distance), cross-shore current; all standardised within the
#                species' study area.
#   (a) Binomial GLM with RMU fixed effects; standard errors clustered by 1-degree cell.
#   (b) Conditional logistic regression stratified by 1-degree cell.
#   (c) Leave-one-RMU-out refits of (a)                                                      
#   (d) Alternative constructions of the moonlight score in (a) and (b)       

library(data.table)
library(sandwich)
library(survival)
library(car)
dir.create("results", showWarnings = FALSE)

b <- fread("data/beaches.csv", na.strings = c("", "NA"))
b <- merge(b, fread("data/moonlight_score.csv"), by = "id")
b <- merge(b, fread("data/waves_beach.csv")[, .(id, wave_pow_on)], by = "id")
b <- merge(b, fread("data/wind_beach.csv")[, .(id, wind_on)], by = "id")
b <- merge(b, fread("data/currents_beach.csv")[, .(id, dist_strong05, cur_cross)], by = "id")

z  <- function(x) (x - mean(x)) / sd(x)
sp_names <- c(CC = "C. caretta", CM = "C. mydas", LO = "L. olivacea", EI = "E. imbricata", DC = "D. coriacea")
lab <- c(moonlight = "Seaward moonlight", waves = "Onshore wave exposure", wind = "Onshore wind exposure",
         near_current = "Proximity to strong current", cross_current = "Cross-shore current")
preds <- names(lab)
variants <- c(ml_3h = "0-3 h after sunset", ml_hatch = "Hatching season by hemisphere", ml_lowmoon = "Low moon weighted (cos altitude)")

fit_rmu  <- function(d, terms) { f <- glm(reformulate(c(terms, "rmu"), "y"), binomial, d)
  list(coef = coef(f), vcov = vcovCL(f, cluster = ~cell1, type = "HC0")) }
fit_cell <- function(d, terms) { k <- clogit(as.formula(paste("y ~", paste(terms, collapse = " + "), "+ strata(cell1) + cluster(sub025)")),
                                            d, method = "approximate")
  list(coef = coef(k), vcov = vcov(k)) }
ci <- function(m, term) { e <- m$coef[term]; s <- sqrt(m$vcov[term, term])
  data.table(est = e, lo = e - 1.96 * s, hi = e + 1.96 * s, p = 2 * pnorm(-abs(e / s))) }

T1 <- T2 <- S1 <- S2 <- R <- list()
for (sp in names(sp_names)) {
  d <- b[!is.na(get(paste0("rmu_", sp))) & !is.na(wave_pow_on) & !is.na(cur_cross)]
  n_rmu <- nrow(d)
  r <- range(d[get(paste0("nest_", sp)) == 1, sst_mean], na.rm = TRUE)          # temperature filter
  d <- d[!is.na(sst_mean) & sst_mean >= r[1] & sst_mean <= r[2]]
  d[, `:=`(y = get(paste0("nest_", sp)), rmu = factor(get(paste0("rmu_", sp))),
           cell1 = paste(floor(lat), floor(lon)), sub025 = paste(floor(lat * 4), floor(lon * 4)),
           moonlight = z(moonlight), waves = z(log1p(wave_pow_on)), wind = z(log1p(wind_on)),
           near_current = z(-log(dist_strong05)), cross_current = z(cur_cross))]
  for (v in names(variants)) set(d, j = v, value = z(d[[v]]))
  inf <- d[, .(k = sum(y), n = .N), by = cell1][k > 0 & k < n]$cell1             # informative cells
  dc  <- d[cell1 %in% inf]

  T1[[sp]] <- data.table(species = sp_names[sp], sst_range = sprintf("%.1f-%.1f", r[1], r[2]), segments_in_rmus = n_rmu,
                         segments_after_sst = nrow(d), nesting_segments = sum(d$y), rmus = uniqueN(d[y == 1]$rmu), cells = length(inf))
  # (a, b) main models
  mr <- fit_rmu(d, preds); mc <- fit_cell(dc, preds)
  T2[[sp]] <- rbindlist(lapply(preds, function(t) cbind(data.table(species = sp_names[sp], variable = lab[t]),
                                                         setnames(ci(mr, t), paste0("rmu_", names(ci(mr, t)))),
                                                         setnames(ci(mc, t), paste0("cell_", names(ci(mc, t)))))))
  # predictor correlations and variance inflation factors
  cc <- cor(d[, ..preds]); R[[sp]] <- data.table(species = sp_names[sp], var1 = rep(preds, each = 5), var2 = rep(preds, 5), r = as.vector(cc))
  T2[[sp]][, VIF := vif(lm(reformulate(preds, "y"), d))[preds]]
  # (c) leave one RMU out
  S1[[sp]] <- rbind(cbind(data.table(species = sp_names[sp], rmu_left_out = "None", nests_removed = NA_integer_), ci(mr, "moonlight")),
    rbindlist(lapply(levels(d$rmu), function(R0) { if (sum(d[rmu == R0]$y) == 0) return(NULL)
      dd <- d[rmu != R0]; dd[, rmu := droplevels(rmu)]
      cbind(data.table(species = sp_names[sp], rmu_left_out = R0, nests_removed = sum(d[rmu == R0]$y)), ci(fit_rmu(dd, preds), "moonlight")) })))
  # (d) alternative constructions of the moonlight score
  S2[[sp]] <- rbindlist(lapply(c("moonlight", names(variants)), function(v) {
    tt <- c(v, preds[-1]); a <- ci(fit_rmu(d, tt), v); c2 <- ci(fit_cell(dc, tt), v)
    data.table(species = sp_names[sp], score = if (v == "moonlight") "Main score (0-6 h, all year)" else variants[v],
               rmu_est = a$est, rmu_lo = a$lo, rmu_hi = a$hi, cell_est = c2$est, cell_lo = c2$lo, cell_hi = c2$hi) }))
  message(sp_names[sp], " done")
}
fwrite(rbindlist(T1), "results/table1_sample_sizes.csv")
fwrite(rbindlist(T2), "results/table2.csv")
fwrite(rbindlist(S1), "results/tableS1_leave_one_rmu_out.csv")
fwrite(rbindlist(S2), "results/tableS2_score_variants.csv")
fwrite(rbindlist(R),  "results/predictor_correlations.csv")
print(rbindlist(T2)[variable == "Seaward moonlight", .(species, rmu_est, rmu_p, cell_est, cell_p)], digits = 2)

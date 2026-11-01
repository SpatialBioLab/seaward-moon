# 04_surface_currents.R ---------------------------------------------------------
# Surface-current variables per beach segment, 1994-2023 monthly means.

library(data.table)
glo_file <- "data/currents_glo_clim_monthly.rds"
med_file <- "data/currents_med_clim_monthly.rds"

# ---- download and monthly climatologies ----------------------------------------------------
clim_from_cmems <- function(product, layer, region, depth) {
  library(CopernicusMarine); library(stars); acc <- NULL
  for (y in 1994:2023) {
    s <- cms_download_subset(product = product, layer = layer, variable = c("uo", "vo"), region = region,
                             timerange = c(sprintf("%d-01-01", y), sprintf("%d-12-31", y)), verticalrange = depth, progress = FALSE)
    uo <- drop(s[["uo"]]); vo <- drop(s[["vo"]])                                     # lon x lat x 12
    if (inherits(uo, "units")) { uo <- units::drop_units(uo); vo <- units::drop_units(vo) }
    sp <- sqrt(uo^2 + vo^2); ok <- !is.na(uo); uo[!ok] <- 0; vo[!ok] <- 0; sp[!ok] <- 0
    if (is.null(acc)) { acc <- list(uo = uo, vo = vo, speed = sp); n <- ok * 1
      lon <- st_get_dimension_values(s, 1); lat <- st_get_dimension_values(s, 2)
    } else { acc$uo <- acc$uo + uo; acc$vo <- acc$vo + vo; acc$speed <- acc$speed + sp; n <- n + ok }
    rm(s); gc(); message(product, " ", y)
  }
  n[n == 0] <- NA
  list(lon = lon, lat = lat, uo = acc$uo / n, vo = acc$vo / n, speed = acc$speed / n)
}
# options(CopernicusMarine_uid = "USERNAME", CopernicusMarine_pwd = "PASSWORD")
if (!file.exists(glo_file)) saveRDS(clim_from_cmems("MULTIOBS_GLO_PHY_MYNRT_015_003", "cmems_obs-mob_glo_phy-cur_my_0.25deg_P1M-m",
                                                    c(-180, -45, 180, 55), c(0, -1)), glo_file)
if (!file.exists(med_file)) saveRDS(clim_from_cmems("MEDSEA_MULTIYEAR_PHY_006_004", "med-cmcc-cur-rean-m",
                                                    c(-6, 30, 36.5, 46), c(0, -1.1)), med_file)

# ---- fill the Mediterranean into the global grid ----------------------------------------------------
cu <- readRDS(glo_file); md <- readRDS(med_file)
gx <- round((md$lon - cu$lon[1]) / 0.25) + 1; gy <- round((md$lat - cu$lat[1]) / 0.25) + 1
for (m in 1:12) for (v in c("uo", "vo", "speed")) {
  df <- data.frame(gx = rep(gx, times = length(gy)), gy = rep(gy, each = length(gx)), val = as.vector(md[[v]][, , m]))
  df <- aggregate(val ~ gx + gy, df[!is.na(df$val), ], mean)                       # block average to 0.25 deg
  idx <- cbind(df$gx, df$gy); A <- cu[[v]][, , m]; empty <- is.na(A[idx])
  A[idx[empty, , drop = FALSE]] <- df$val[empty]; cu[[v]][, , m] <- A
}

# ---- annual fields and per-beach extraction -------------------------------------------------
F <- lapply(cu[c("uo", "vo", "speed")], function(a) { x <- apply(a, c(1, 2), mean, na.rm = TRUE); x[!is.finite(x)] <- NA; x })
lon0 <- cu$lon[1]; lat0 <- cu$lat[1]; dg <- 0.25; nx <- length(cu$lon); ny <- length(cu$lat)
look <- function(fld, lon, lat) {
  ix <- (round(((lon - lon0) %% 360) / dg) %% nx) + 1; iy <- round((lat - lat0) / dg) + 1
  out <- rep(NA_real_, length(lon)); ok <- iy >= 1 & iy <= ny; out[ok] <- fld[cbind(ix[ok], iy[ok])]; out
}
dest <- function(lon, lat, r, g) list(lon = lon + r / (111.32 * cos(lat * pi / 180)) * sin(g * pi / 180), lat = lat + r / 111.32 * cos(g * pi / 180))

b <- fread("data/beaches.csv", select = c("id", "lon", "lat", "angle")); n <- nrow(b); a <- b$angle * pi / 180
r_off <- seq(10, 100, 10); r_pol <- seq(12.5, 300, 12.5); brg <- seq(0, 345, 15)
U <- V <- matrix(NA_real_, n, length(r_off))
for (j in seq_along(r_off)) { p <- dest(b$lon, b$lat, r_off[j], b$angle); U[, j] <- look(F$uo, p$lon, p$lat); V[, j] <- look(F$vo, p$lon, p$lat) }
cur_cross <- rowMeans(U * sin(a) + V * cos(a), na.rm = TRUE)                     # projection on the seaward normal
miss <- which(!is.finite(cur_cross)); fb <- rep(FALSE, n)
for (r in r_pol) { if (!length(miss)) break                                         # fallback: nearest ring, any direction
  uu <- sapply(brg, function(g) { p <- dest(b$lon[miss], b$lat[miss], r, g); look(F$uo, p$lon, p$lat) })
  vv <- sapply(brg, function(g) { p <- dest(b$lon[miss], b$lat[miss], r, g); look(F$vo, p$lon, p$lat) })
  if (length(miss) == 1) { uu <- matrix(uu, 1); vv <- matrix(vv, 1) }
  hit <- rowSums(!is.na(uu)) > 0
  if (any(hit)) { i <- miss[hit]
    cur_cross[i] <- rowMeans(uu[hit, , drop = FALSE], na.rm = TRUE) * sin(a[i]) + rowMeans(vv[hit, , drop = FALSE], na.rm = TRUE) * cos(a[i])
    fb[i] <- TRUE; miss <- miss[!hit] } }
d05 <- rep(300, n)                                                                  # distance to strong current
for (r in r_pol) {
  mx <- suppressWarnings(apply(sapply(brg, function(g) { p <- dest(b$lon, b$lat, r, g); look(F$speed, p$lon, p$lat) }), 1, max, na.rm = TRUE))
  d05 <- ifelse(d05 == 300 & is.finite(mx) & mx >= 0.5, r, d05)
}
cur_cross[!is.finite(cur_cross)] <- NA
fwrite(data.table(id = b$id, dist_strong05 = d05, cur_cross = cur_cross, cur_fallback = fb), "data/currents_beach.csv")
cat("beaches:", n, "| fallback used:", sum(fb), "| no strong current within 300 km:", sum(d05 == 300), "\n")

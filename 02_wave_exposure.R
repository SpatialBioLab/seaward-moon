# 02_wave_exposure.R ------------------------------------------------------------
# Onshore wave exposure per beach segment from ERA5 monthly means, 1994-2023 (0.5 degree).

library(data.table)
library(ncdf4)
clim_file <- "data/waves_clim_monthly.rds"

# ---- download and monthly climatology -------------------------------------------------------
if (!file.exists(clim_file)) {
  library(ecmwfr)
  dl <- "data/raw_era5"; dir.create(dl, showWarnings = FALSE, recursive = TRUE); acc <- NULL
  for (y in 1994:2023) {
    f <- file.path(dl, sprintf("era5_waves_%d.nc", y))
    if (!file.exists(f)) wf_request(request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means", product_type = "monthly_averaged_reanalysis",
      variable = c("mean_wave_direction", "significant_height_of_combined_wind_waves_and_swell", "mean_wave_period"),
      year = as.character(y), month = sprintf("%02d", 1:12), time = "00:00",
      area = c(55, -180, -45, 180), grid = c(0.5, 0.5), data_format = "netcdf", download_format = "unarchived",
      target = basename(f)), path = dl)
    nc <- nc_open(f); lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude")
    mwd <- ncvar_get(nc, "mwd"); swh <- ncvar_get(nc, "swh"); mwp <- ncvar_get(nc, "mwp"); nc_close(nc)   # lon x lat x 12
    w <- swh^2; th <- mwd * pi / 180
    part <- lapply(list(ws = w * sin(th), wc = w * cos(th), w = w, swh = swh, mwp = mwp), function(a) { a[is.na(a)] <- 0; a })
    ok <- !is.na(swh)
    if (is.null(acc)) { acc <- part; n <- ok * 1 } else { for (k in names(acc)) acc[[k]] <- acc[[k]] + part[[k]]; n <- n + ok }
    file.remove(f); message("ERA5 ", y)
  }
  n[n == 0] <- NA
  clim <- list(lon = lon, lat = lat, mwd = (atan2(acc$ws, acc$wc) * 180 / pi) %% 360,
               cons = sqrt(acc$ws^2 + acc$wc^2) / acc$w, swh = acc$swh / n, mwp = acc$mwp / n)
  clim$mwd[is.na(clim$swh)] <- NA; clim$cons[is.na(clim$swh)] <- NA
  saveRDS(clim, clim_file)
}

# ---- annual fields and per-beach extraction -------------------------------------------------
C <- readRDS(clim_file)
w <- C$swh^2; th <- C$mwd * pi / 180                                  # energy-weighted mean over the 12 months
ws <- apply(w * C$cons * sin(th), 1:2, mean, na.rm = TRUE); wc <- apply(w * C$cons * cos(th), 1:2, mean, na.rm = TRUE)
F <- list(dir = (atan2(ws, wc) * 180 / pi) %% 360, swh = apply(C$swh, 1:2, mean, na.rm = TRUE), mwp = apply(C$mwp, 1:2, mean, na.rm = TRUE))
for (k in names(F)) F[[k]][!is.finite(F[[k]])] <- NA
dg <- 0.5; lon0 <- C$lon[1]; nx <- length(C$lon); lat_desc <- C$lat[1] > C$lat[2]
look <- function(fld, lo, la) {                                       # nearest grid cell; longitude wraps
  ix <- (round(((lo - lon0) %% 360) / dg) %% nx) + 1
  iy <- if (lat_desc) round((C$lat[1] - la) / dg) + 1 else round((la - C$lat[1]) / dg) + 1
  out <- rep(NA_real_, length(lo)); ok <- iy >= 1 & iy <= length(C$lat); out[ok] <- fld[cbind(ix[ok], iy[ok])]; out
}
dest <- function(lo, la, r, g) list(lon = lo + r / (111.32 * cos(la * pi / 180)) * sin(g * pi / 180), lat = la + r / 111.32 * cos(g * pi / 180))

b <- fread("data/beaches.csv", select = c("id", "lon", "lat", "angle")); n <- nrow(b)
acc <- list(sn = numeric(n), cs = numeric(n), swh = numeric(n), mwp = numeric(n), k = numeric(n))
add <- function(i, p) {
  d <- look(F$dir, p$lon, p$lat); ok <- !is.na(d); j <- i[ok]
  acc$sn[j]  <<- acc$sn[j]  + sin(d[ok] * pi / 180); acc$cs[j] <<- acc$cs[j] + cos(d[ok] * pi / 180)
  acc$swh[j] <<- acc$swh[j] + look(F$swh, p$lon, p$lat)[ok]; acc$mwp[j] <<- acc$mwp[j] + look(F$mwp, p$lon, p$lat)[ok]
  acc$k[j]   <<- acc$k[j] + 1
}
for (r in seq(10, 100, 10)) add(seq_len(n), dest(b$lon, b$lat, r, b$angle))     # along the seaward normal
fb <- rep(FALSE, n)
for (r in seq(25, 200, 25)) {                                                    # fallback: nearest ring, any direction
  miss <- which(acc$k == 0); if (!length(miss)) break
  for (g in seq(0, 345, 15)) add(miss, dest(b$lon[miss], b$lat[miss], r, g)); fb[miss[acc$k[miss] > 0]] <- TRUE
}
k <- ifelse(acc$k > 0, acc$k, NA)
dir_m <- (atan2(acc$sn, acc$cs) * 180 / pi) %% 360; dir_m[is.na(k)] <- NA
inc <- cos((dir_m - b$angle) * pi / 180)
out <- data.table(id = b$id, wave_inc = inc, wave_pow_on = (acc$swh / k)^2 * (acc$mwp / k) * pmax(0, inc), wave_fallback = fb)
fwrite(out, "data/waves_beach.csv")
cat("beaches:", n, "| fallback used:", sum(fb), "| missing:", sum(is.na(out$wave_pow_on)), "\n")

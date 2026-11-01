# 03_wind_exposure.R -------------------------------------------------------------
# Onshore wind exposure per beach segment from ERA5 monthly means (0.25 degree).

library(data.table)
library(ncdf4)
clim_file <- "data/wind_clim.rds"

# ---- download and long-term means ----------------------------------------------------------
if (!file.exists(clim_file)) {
  library(ecmwfr)
  dl <- "data/raw_era5"; dir.create(dl, showWarnings = FALSE, recursive = TRUE); acc <- NULL; n <- 0
  for (y in 1994:2023) {
    f <- file.path(dl, sprintf("era5_wind_%d.nc", y))
    if (!file.exists(f)) wf_request(request = list(
      dataset_short_name = "reanalysis-era5-single-levels-monthly-means", product_type = "monthly_averaged_reanalysis",
      variable = c("10m_u_component_of_wind", "10m_v_component_of_wind", "10m_wind_speed"),
      year = as.character(y), month = sprintf("%02d", 1:12), time = "00:00",
      area = c(48, -180, -39, 180), data_format = "netcdf", download_format = "unarchived",
      target = basename(f)), path = dl)
    nc <- nc_open(f); lon <- ncvar_get(nc, "longitude"); lat <- ncvar_get(nc, "latitude")
    part <- list(u = ncvar_get(nc, "u10"), v = ncvar_get(nc, "v10"), s = ncvar_get(nc, "si10")); nc_close(nc)   # lon x lat x 12
    part <- lapply(part, function(a) apply(a, 1:2, sum))
    if (is.null(acc)) acc <- part else for (k in names(acc)) acc[[k]] <- acc[[k]] + part[[k]]
    n <- n + 12; file.remove(f); message("ERA5 wind ", y)
  }
  saveRDS(list(lon = lon, lat = lat, u = acc$u / n, v = acc$v / n, s = acc$s / n), clim_file)
}

# ---- per-beach extraction -------------------------------------------------------------------
C <- readRDS(clim_file)
dg <- 0.25; lon0 <- C$lon[1]; nx <- length(C$lon); lat_desc <- C$lat[1] > C$lat[2]
look <- function(fld, lo, la) {                                       # nearest grid cell; longitude wraps
  ix <- (round(((lo - lon0) %% 360) / dg) %% nx) + 1
  iy <- if (lat_desc) round((C$lat[1] - la) / dg) + 1 else round((la - C$lat[1]) / dg) + 1
  out <- rep(NA_real_, length(lo)); ok <- iy >= 1 & iy <= length(C$lat); out[ok] <- fld[cbind(ix[ok], iy[ok])]; out
}
dest <- function(lo, la, r, g) list(lon = lo + r / (111.32 * cos(la * pi / 180)) * sin(g * pi / 180), lat = la + r / 111.32 * cos(g * pi / 180))

b <- fread("data/beaches.csv", select = c("id", "lon", "lat", "angle")); n <- nrow(b)
acc <- list(u = numeric(n), v = numeric(n), s = numeric(n), k = numeric(n))
add <- function(i, p) {
  u <- look(C$u, p$lon, p$lat); ok <- !is.na(u); j <- i[ok]
  acc$u[j] <<- acc$u[j] + u[ok]; acc$v[j] <<- acc$v[j] + look(C$v, p$lon, p$lat)[ok]
  acc$s[j] <<- acc$s[j] + look(C$s, p$lon, p$lat)[ok]; acc$k[j] <<- acc$k[j] + 1
}
for (r in seq(10, 100, 10)) add(seq_len(n), dest(b$lon, b$lat, r, b$angle))     # along the seaward normal
for (r in seq(25, 200, 25)) {                                                    # fallback: nearest ring, any direction
  miss <- which(acc$k == 0); if (!length(miss)) break
  for (g in seq(0, 345, 15)) add(miss, dest(b$lon[miss], b$lat[miss], r, g))
}
k <- ifelse(acc$k > 0, acc$k, NA)
u <- acc$u / k; v <- acc$v / k; U <- acc$s / k
dir_from <- (atan2(-u, -v) * 180 / pi) %% 360
inc <- cos((dir_from - b$angle) * pi / 180)
out <- data.table(id = b$id, wind_inc = inc, wind_on = U^2 * pmax(0, inc))
fwrite(out, "data/wind_beach.csv")
cat("beaches:", n, "| missing:", sum(is.na(out$wind_on)), "\n")

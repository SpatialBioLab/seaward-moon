# 01_seaward_moonlight.R -------------------------------------------------------
# Long-term seaward moonlight score for every beach segment.
#
# For every night 1900-2024 and every 2-degree latitude band (longitude 0), the moon's azimuth, altitude and
# illuminated fraction are computed every 30 min from sunset to 6 h after sunset (or to sunrise, if earlier).

library(data.table)
library(suncalc)
library(parallel)

bands <- seq(-40, 44, 2)                                      # band centres; beaches outside go to the nearest band
days  <- seq(as.Date("1900-01-01"), as.Date("2024-12-31"), by = "day")
mon   <- as.integer(format(days, "%m"))
azdeg <- function(a) (a * 180 / pi + 180) %% 360              # suncalc azimuth (from south, radians) -> degrees from north
h360  <- function(az, w) { h <- numeric(360); if (length(az)) { t <- tapply(w, floor(az) %% 360, sum); h[as.integer(names(t)) + 1] <- t }; h }

# Monthly azimuth histograms per band: 360 one-degree bins x 12 months, for three weightings and two windows
one_band <- function(lat) {
  st <- getSunlightTimes(date = days,     lat = lat, lon = 0, keep = "sunset",  tz = "UTC")$sunset
  sr <- getSunlightTimes(date = days + 1, lat = lat, lon = 0, keep = "sunrise", tz = "UTC")$sunrise
  night <- as.numeric(difftime(sr, st, units = "hours"))
  H <- list(w6 = matrix(0, 12, 360), w3 = matrix(0, 12, 360), low = matrix(0, 12, 360)); N6 <- N3 <- numeric(12)
  for (s in seq(0.25, 5.75, 0.5)) {                           # sample at the middle of each 30-min step
    ok <- s <= night; t <- st[ok] + s * 3600; mm <- mon[ok]
    N6 <- N6 + tabulate(mm, 12); if (s < 3) N3 <- N3 + tabulate(mm, 12)
    p  <- getMoonPosition(t, lat, 0); up <- p$altitude > 0
    fr <- getMoonIllumination(t[up], keep = "fraction")$fraction
    az <- azdeg(p$azimuth[up]); alt <- p$altitude[up]; mu <- mm[up]
    for (m in 1:12) { i <- mu == m
      H$w6[m, ]  <- H$w6[m, ]  + h360(az[i], fr[i])
      H$low[m, ] <- H$low[m, ] + h360(az[i], fr[i] * cos(alt[i]))
      if (s < 3) H$w3[m, ] <- H$w3[m, ] + h360(az[i], fr[i]) }
  }
  list(H = H, N6 = N6, N3 = N3)
}
res <- mclapply(bands, one_band, mc.cores = max(1, detectCores() - 1)); names(res) <- bands

# Score per beach: histogram weighted by max(0, cos(normal - azimuth)), divided by the number of samples
b <- fread("data/beaches.csv", select = c("id", "lat", "angle"))
b[, band := as.character(pmin(max(bands), pmax(min(bands), 2 * round(lat / 2))))]
cosm <- function(a) { m <- cos(outer(a, 0:359 + 0.5, "-") * pi / 180); m[m < 0] <- 0; m }
hatch_months <- function(lat) if (lat >= 10) 7:10 else if (lat <= -10) 1:4 else 1:12
score <- function(hist, nkey, by_season = FALSE) {
  v <- numeric(nrow(b))
  for (B in unique(b$band)) {
    i <- which(b$band == B); x <- res[[B]]
    if (!by_season) { v[i] <- as.vector(cosm(b$angle[i]) %*% colSums(x$H[[hist]])) / sum(x[[nkey]]); next }
    grp <- split(i, sapply(b$lat[i], function(l) paste(hatch_months(l), collapse = ",")))
    for (g in grp) { m <- hatch_months(b$lat[g[1]])
      v[g] <- as.vector(cosm(b$angle[g]) %*% colSums(x$H[[hist]][m, , drop = FALSE])) / sum(x[[nkey]][m]) }
  }
  v
}
b[, `:=`(moonlight = score("w6", "N6"), ml_3h = score("w3", "N3"),
         ml_hatch = score("w6", "N6", by_season = TRUE), ml_lowmoon = score("low", "N6"))]
cat("Correlation of variants with the main score:\n"); print(round(cor(b[, .(moonlight, ml_3h, ml_hatch, ml_lowmoon)])[1, ], 3))
fwrite(b[, .(id, moonlight, ml_3h, ml_hatch, ml_lowmoon)], "data/moonlight_score.csv")

## Generate the forecastbaselines hex logo.
##
## Run from the package root:
##   Rscript data-raw/logo.R
##   inkscape man/figures/logo.svg -o man/figures/logo.png -w 240
## Concept: the three baseline families the package provides, all leaving the
## same last observation -- persistence holds the level, drift carries the
## trend, seasonal repeats last season's shape. The three colours are Julia's,
## since this is the R interface to ForecastBaselines.jl.

W <- 1732L; H <- 2000L                   # pointy-top hex canvas

## --- palette ---------------------------------------------------------------
bg       <- "#FFFFFF"
bg_foot  <- "#EEF3F8"
border   <- "#3D4E63"
ink      <- "#16212A"
axis_col <- "#CBD6E2"

## Julia's three dots, one per baseline
persist_col <- "#CB3C33"
drift_col   <- "#389826"
season_col  <- "#9558B2"

## --- geometry (data coords, y up; svg y flipped about `ybase`) --------------
x0 <- 215; xf <- 880; x1 <- 1470         # start, forecast date, horizon
ybase <- 1190
flip <- function(y) ybase - y
pt <- function(x, y) sprintf("%.1f,%.1f", x, flip(y))

## the series: a rising level with a clear season riding on it
level <- function(u) 250 + 235 * u
season <- function(u) 140 * sin(2 * pi * 2.25 * u)

uh <- seq(0, 1, length.out = 150)
hx <- x0 + uh * (xf - x0)
hy <- level(uh) + season(uh)

last_u <- 1
last_y <- level(last_u) + season(last_u)
slope <- 235 / (xf - x0)                 # the trend, per x unit

uf <- seq(0, 1, length.out = 90)
fx <- xf + uf * (x1 - xf)

## persistence: hold the last value
persistence <- rep(last_y, length(uf))
## drift: carry the trend on, without the season
drift <- last_y + slope * (fx - xf)
## seasonal: repeat last season's shape from the same level
season_fc <- last_y + (season(last_u + uf * 0.62) - season(last_u))

## every baseline is probabilistic: the interval opens with the horizon
cone <- function(scale) 18 + scale * sqrt(uf)

## --- svg helpers ------------------------------------------------------------
line <- function(x, y, col, w, dash = NULL) sprintf(
  '<polyline points="%s" fill="none" stroke="%s" stroke-width="%.0f"%s stroke-linecap="round" stroke-linejoin="round"/>',
  paste(pt(x, y), collapse = " "), col, w,
  if (is.null(dash)) "" else sprintf(' stroke-dasharray="%s"', dash))

band <- function(x, mid, half, col, op) sprintf(
  '<path d="M %s L %s Z" fill="%s" fill-opacity="%.2f"/>',
  paste(pt(x, mid - half), collapse = " L "),
  paste(rev(pt(x, mid + half)), collapse = " L "), col, op)

hex <- sprintf("%d,0 %d,500 %d,1500 %d,%d 0,1500 0,500", W %/% 2L, W, W, W %/% 2L, H)

svg <- c(
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">', W, H, W, H),
  '<defs><linearGradient id="g" x1="0" y1="0" x2="0" y2="1">',
  sprintf('<stop offset="0" stop-color="%s"/><stop offset="1" stop-color="%s"/>', bg, bg_foot),
  '</linearGradient>',
  sprintf('<clipPath id="clip"><polygon points="%s"/></clipPath>', hex),
  '</defs>',
  sprintf('<polygon points="%s" fill="url(#g)"/>', hex),
  '<g clip-path="url(#clip)">',
  sprintf('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="%s" stroke-width="8" stroke-linecap="round"/>',
          x0 - 16, flip(40), x1 + 16, flip(40), axis_col),
  sprintf('<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="%s" stroke-width="8" stroke-dasharray="24 26"/>',
          xf, flip(870), xf, flip(60), axis_col),
  ## the intervals, widest first so the narrower ones stay readable
  band(fx, season_fc, cone(132), season_col, 0.18),
  band(fx, drift, cone(150), drift_col, 0.18),
  band(fx, persistence, cone(118), persist_col, 0.18),
  line(fx, season_fc, season_col, 26),
  line(fx, drift, drift_col, 26),
  line(fx, persistence, persist_col, 26),
  line(hx, hy, ink, 30),
  sprintf('<circle cx="%.0f" cy="%.0f" r="30" fill="%s"/>', xf, flip(last_y), ink),
  sprintf('<text x="%d" y="1600" text-anchor="middle" font-family="Fira Sans, DejaVu Sans, Helvetica, Arial, sans-serif" font-size="158" font-weight="600" letter-spacing="1" fill="%s">forecast<tspan fill="%s">baselines</tspan></text>',
          W %/% 2L, ink, season_col),
  '</g>',
  sprintf('<polygon points="%s" fill="none" stroke="%s" stroke-width="44" stroke-linejoin="round"/>', hex, border),
  '</svg>')

writeLines(svg, "man/figures/logo.svg")

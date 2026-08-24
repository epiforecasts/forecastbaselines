## Generate the forecastbaselines hex logo.
##
## Run from the package root:
##   Rscript data-raw/logo.R
##   inkscape man/figures/logo.svg -o man/figures/logo.png -w 240
##
## The graphic is the three baseline families, all leaving the same last
## observation: persistence holds the level, drift carries the trend, and
## seasonal repeats last season's shape. Each carries an interval that opens
## with the horizon. The three colours are Julia's, since this is the R
## interface to ForecastBaselines.jl.

## --- canvas ----------------------------------------------------------------

width <- 1732L
height <- 2000L
hex <- sprintf(
  "%d,0 %d,500 %d,1500 %d,%d 0,1500 0,500",
  width %/% 2L, width, width, width %/% 2L, height
)

## --- palette ---------------------------------------------------------------

bg <- "#FFFFFF"
bg_foot <- "#EEF3F8"
border <- "#3D4E63"
ink <- "#16212A"
axis_col <- "#CBD6E2"

## Julia's three dots, one per baseline
persist_col <- "#CB3C33"
drift_col <- "#389826"
season_col <- "#9558B2"

## --- geometry --------------------------------------------------------------

x_start <- 215
x_split <- 880
x_end <- 1470
ybase <- 1190

flip <- function(y) {
  ybase - y
}

pt <- function(x, y) {
  sprintf("%.1f,%.1f", x, flip(y))
}

## the series: a rising level with a clear season riding on it
level <- function(u) {
  250 + 235 * u
}

season <- function(u) {
  140 * sin(2 * pi * 2.25 * u)
}

u_hist <- seq(0, 1, length.out = 150)
hist_x <- x_start + u_hist * (x_split - x_start)
hist_y <- level(u_hist) + season(u_hist)

last_y <- level(1) + season(1)
slope <- 235 / (x_split - x_start)

u_fc <- seq(0, 1, length.out = 90)
fc_x <- x_split + u_fc * (x_end - x_split)

## persistence holds the last value; drift carries the trend on without the
## season; seasonal repeats last season's shape from the same level
persistence <- rep(last_y, length(u_fc))
drift <- last_y + slope * (fc_x - x_split)
seasonal <- last_y + (season(1 + u_fc * 0.62) - season(1))

## every baseline is probabilistic: the interval opens with the horizon
cone <- function(scale) {
  18 + scale * sqrt(u_fc)
}

## --- svg parts -------------------------------------------------------------

line <- function(x, y, col, w) {
  sprintf(
    paste0(
      '<polyline points="%s" fill="none" stroke="%s" ',
      'stroke-width="%.0f" stroke-linecap="round" stroke-linejoin="round"/>'
    ),
    paste(pt(x, y), collapse = " "), col, w
  )
}

band <- function(x, mid, half, col) {
  sprintf(
    '<path d="M %s L %s Z" fill="%s" fill-opacity="0.18"/>',
    paste(pt(x, mid - half), collapse = " L "),
    paste(rev(pt(x, mid + half)), collapse = " L "),
    col
  )
}

rule <- function(x1, y1, x2, y2, dash) {
  sprintf(
    paste0(
      '<line x1="%.0f" y1="%.0f" x2="%.0f" y2="%.0f" stroke="%s" ',
      'stroke-width="8" stroke-dasharray="%s" stroke-linecap="round"/>'
    ),
    x1, flip(y1), x2, flip(y2), axis_col, dash
  )
}

## --- assemble --------------------------------------------------------------

svg <- c(
  sprintf(
    paste0(
      '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" ',
      'viewBox="0 0 %d %d">'
    ),
    width, height, width, height
  ),
  "<defs>",
  '<linearGradient id="g" x1="0" y1="0" x2="0" y2="1">',
  sprintf('<stop offset="0" stop-color="%s"/>', bg),
  sprintf('<stop offset="1" stop-color="%s"/>', bg_foot),
  "</linearGradient>",
  sprintf('<clipPath id="clip"><polygon points="%s"/></clipPath>', hex),
  "</defs>",
  sprintf('<polygon points="%s" fill="url(#g)"/>', hex),
  ## clipped so an interval can never bleed past the border
  '<g clip-path="url(#clip)">',
  rule(x_start - 16, 40, x_end + 16, 40, "none"),
  rule(x_split, 870, x_split, 60, "24 26"),
  ## the intervals, widest first so the narrower ones stay readable
  band(fc_x, seasonal, cone(132), season_col),
  band(fc_x, drift, cone(150), drift_col),
  band(fc_x, persistence, cone(118), persist_col),
  line(fc_x, seasonal, season_col, 26),
  line(fc_x, drift, drift_col, 26),
  line(fc_x, persistence, persist_col, 26),
  ## observed history, on top
  line(hist_x, hist_y, ink, 30),
  ## where all three leave from
  sprintf(
    '<circle cx="%.0f" cy="%.0f" r="30" fill="%s"/>',
    x_split, flip(last_y), ink
  ),
  "</g>",
  sprintf(
    paste0(
      '<text x="%d" y="1600" text-anchor="middle" ',
      'font-family="Fira Sans, DejaVu Sans, Helvetica, Arial, sans-serif" ',
      'font-size="158" font-weight="600" letter-spacing="1" fill="%s">',
      'forecast<tspan fill="%s">baselines</tspan></text>'
    ),
    width %/% 2L, ink, season_col
  ),
  sprintf(
    paste0(
      '<polygon points="%s" fill="none" stroke="%s" ',
      'stroke-width="44" stroke-linejoin="round"/>'
    ),
    hex, border
  ),
  "</svg>"
)

writeLines(svg, file.path("man", "figures", "logo.svg"))

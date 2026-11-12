# =====================================================================
#  CROWDING SENSITIVITY (Table 4.9)
#  The crowding factor f in Stages 2-3 is assumed, so the what-if scenarios
#  are repeated with f = 1.5, 2 and 3. Uses the live reading logged at 20:00
#  on 25 September 2026 (the reading shown in Figures 4.11 and 4.12).
#
# AI USE STATEMENT: Generative AI (Claude) was used to help debug, refine and
# structure the R code in this file. Generative AI was not used to conceive the
# project or to choose its aims. All modelling decisions, analyses and
# interpretations are the author's own work, and the author takes full
# responsibility for this code.
# =====================================================================

source("app.R")      # settings and the simulate() engine; the dashboard itself is not started

# The 25 September 2026 reading and that day's hourly average, from the twin's own logs
live_log  <- read.csv("data/live_log.csv")
busy_log  <- read.csv("data/busyness_log.csv")
live <- list(live     = live_log[live_log$scrape_time == "2026-09-25 20:00:14", ][1, ],
             busyness = busy_log[busy_log$date == "2026-09-25", ])
f <- arrival_rates(live, as.POSIXct("2026-09-25 21:00:00", tz = "Australia/Brisbane"))
cat(sprintf("Live ratio applied: %.2f\n", f$ratio))

scenarios <- list("Device use 0%"       = list(device_mult = 0),
                  "Device use 50%"      = list(device_mult = 0.5),
                  "Device use 150%"     = list(device_mult = 1.5),
                  "ED arrivals 50%"     = list(arrivals_mult = 0.5),
                  "ED arrivals 150%"    = list(arrivals_mult = 1.5),
                  "ICU transfers 0%"    = list(icu_mult = 0),
                  "ICU transfers 200%"  = list(icu_mult = 2))

# Expected infections over 72 h (mean over 1,000 runs), and the % change of each scenario vs baseline
expected <- function(...) mean(rowSums(simulate(f$rate, f$times, ...)$exp_inf))
change_table <- function(stage, crowd_factor) {
  base <- expected(stage = stage, crowd_factor = crowd_factor)
  changes <- sapply(scenarios, function(s)
    (do.call(expected, c(list(stage = stage, crowd_factor = crowd_factor), s)) / base - 1) * 100)
  c(Baseline = base, changes)
}

runs <- list("Stage 1"          = c(1, 2),
             "Stage 2 (f = 1.5)" = c(2, 1.5), "Stage 2 (f = 2)" = c(2, 2), "Stage 2 (f = 3)" = c(2, 3),
             "Stage 3 (f = 1.5)" = c(3, 1.5), "Stage 3 (f = 2)" = c(3, 2), "Stage 3 (f = 3)" = c(3, 3))
result <- sapply(runs, function(r) change_table(r[1], r[2]))
print(round(result, 2))
write.csv(round(result, 2), "data/crowding_sensitivity.csv")

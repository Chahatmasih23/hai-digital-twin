# =====================================================================
# LIVE DATA SCRAPER - Redcliffe Hospital ED (Queensland Open Hospitals)
# =====================================================================
# What it collects:
#   - median waiting time (minutes)
#   - number of patients currently waiting
#   - treatment spaces
#   - "busyness" (the purple chart) = 24 hourly averages of patients waiting, for the same
#     weekday and month, over the past 2 years (official definition:
#     https://www.health.qld.gov.au/open-hospitals-about-the-data)
#
# Used in two ways:
#   1. app.R calls scrape_redcliffe() to show the live numbers.
#   2. Run this file on a schedule to build up a history of readings:
#        Rscript scrape_live.R
#
# AI USE STATEMENT: Generative AI (Claude) was used to help debug, refine and structure the R code in this
# file. Generative AI was not used to conceive the project or to choose its aims. All modelling decisions,
# analyses and interpretations are the author's own work, and the author takes full responsibility for this code.
# =====================================================================

URL <- "https://openhospitals.health.qld.gov.au/facility/redcliffe-hospital"

# Folder where the logs are saved (same folder as this script, /data)
DATA_DIR <- "data"


# ---- Small helper: find a value in the page text ----
# The website hides its data inside the page as text like  "patients_waiting":11
grab <- function(page, field) {
  pattern <- paste0('"', field, '":"?([^",}]*)')
  hit <- regmatches(page, regexec(pattern, page))[[1]]
  if (length(hit) < 2) return(NA)
  hit[2]
}


# ---- Main function: read the page once and return one row of data ----
scrape_redcliffe <- function() {
  page <- paste(readLines(URL, warn = FALSE, encoding = "UTF-8"), collapse = "")
  page <- gsub("\\\\", "", page)          # remove escape backslashes  \"  ->  "

  # Median wait time comes as text, e.g. "8 minutes" -> keep the number only
  wait_text <- grab(page, "wait_time_all_patients")
  wait_min  <- as.numeric(gsub("[^0-9]", "", wait_text))

  # The purple chart: 24 numbers inside  "busyness":[3.1,2.1, ... ]
  busy_text <- regmatches(page, regexec('"busyness":\\[([^]]*)\\]', page))[[1]][2]
  busyness  <- as.numeric(strsplit(busy_text, ",")[[1]])

  # Always use Brisbane time (same clock as the website), whatever the laptop's timezone
  now <- as.POSIXlt(Sys.time(), tz = "Australia/Brisbane")
  list(
    live = data.frame(
      scrape_time       = format(now, "%Y-%m-%d %H:%M:%S"),
      website_updated   = sub("Last Updated: ", "", grab(page, "data_update_time")),
      weekday           = format(now, "%A"),
      hour              = as.integer(format(now, "%H")),
      median_wait_min   = wait_min,
      patients_waiting  = as.numeric(grab(page, "patients_waiting")),
      treatment_spaces  = as.numeric(grab(page, "treatment_spaces"))
    ),
    busyness = data.frame(
      date     = format(now, "%Y-%m-%d"),
      weekday  = format(now, "%A"),
      month    = format(now, "%B"),
      hour     = 0:23,
      avg_waiting = busyness
    )
  )
}


# ---- Save one scrape to the log files (adds new rows, never overwrites) ----
save_scrape <- function(result, data_dir = DATA_DIR) {
  dir.create(data_dir, showWarnings = FALSE)

  live_file <- file.path(data_dir, "live_log.csv")

  # skip if the website hasn't updated since the last saved row (avoids duplicates)
  if (file.exists(live_file)) {
    last <- tail(read.csv(live_file), 1)
    if (last$website_updated == result$live$website_updated) return(invisible(NULL))
  }
  write.table(result$live, live_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(live_file), append = file.exists(live_file))

  # The purple chart only changes once a day, so save it once per date
  busy_file <- file.path(data_dir, "busyness_log.csv")
  old <- if (file.exists(busy_file)) read.csv(busy_file) else NULL
  if (is.null(old) || !(result$busyness$date[1] %in% old$date)) {
    write.table(result$busyness, busy_file, sep = ",", row.names = FALSE,
                col.names = is.null(old), append = !is.null(old))
  }
}


# ---- When run as a script (Rscript scrape_live.R): scrape once and save ----
if (sys.nframe() == 0) {
  # make paths relative to this script's folder
  args <- commandArgs(trailingOnly = FALSE)
  script <- sub("--file=", "", args[grep("--file=", args)])
  if (length(script) == 1) setwd(dirname(normalizePath(script)))

  result <- tryCatch(scrape_redcliffe(), error = function(e) NULL)
  if (is.null(result)) {
    cat("Scrape failed (website down or no internet)\n")
  } else {
    save_scrape(result)
    print(result$live)
  }
}

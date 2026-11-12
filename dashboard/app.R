# =====================================================================
#  HAI DIGITAL TWIN  -  E. coli treated with ceftriaxone  (R Shiny dashboard)
#  Reads the files saved by ABM_DT.ipynb (data/) and the live Redcliffe ED feed (scrape_live.R),
#  and forecasts ward occupancy and infections for the next 72 hours.
# =====================================================================
#
# AI USE STATEMENT: Generative AI (Claude) was used to help debug, refine and structure the R code in this
# file. Generative AI was not used to conceive the project or to choose its aims. All modelling decisions,
# analyses and interpretations are the author's own work, and the author takes full responsibility for this code.
# =====================================================================

library(shiny)
library(plotly)
source("scrape_live.R")        # gives us scrape_redcliffe() and save_scrape()


# ---------------------------------------------------------------------
# 1. SETTINGS
# ---------------------------------------------------------------------
HOURS        <- 72     # forecast length: 3 days, in hourly steps
N_RUNS       <- 1000    # how many times each forecast is simulated (more = smoother, but slower)
LIVE_HOURS   <- 6      # how many hours ahead the live "busier / quieter than usual" signal is used
N_TRACK      <- 20     # how many runs are recorded patient-by-patient for the hospital animation

# ---------------------------------------------------------------------
# 2. LOAD WHAT THE NOTEBOOK SAVED
# ---------------------------------------------------------------------
# Stage 1: Poisson Model 2 (ward device utilisation)
coef_table  <- read.csv("data/coefficients.csv")
coef        <- setNames(coef_table$coef, coef_table$term)          # Intercept, device_utilisation

# Stage 2: crowding factor and the scaling k that keeps the overall rate unchanged
crowd_table <- read.csv("data/stage2_crowding.csv", row.names = 1)
crowd       <- setNames(crowd_table$value, rownames(crowd_table))  # crowd_factor, k, data_irr, crowded_share

# k for any crowding factor f: k = 1 / (1 + (f - 1) * share of crowded patient-hours in training)
crowd_k <- function(f) 1 / (1 + (f - 1) * crowd["crowded_share"])

# Stage 3: patient-level model (own device, age 65+, male)
s3_table    <- read.csv("data/stage3_coefficients.csv", row.names = 1)
coef3       <- setNames(s3_table$coef, rownames(s3_table))         # Intercept, own_device, age65, male

# Share of arriving patients who are 65+ and male (to give new agents a profile)
profile_tbl <- read.csv("data/arrival_profile.csv", row.names = 1)
arrival_profile <- setNames(profile_tbl$share, rownames(profile_tbl))   # age65, male

ward_params <- read.csv("data/ward_params.csv", row.names = 1)     # device_prob, crowd_threshold per ward
WARDS       <- rownames(ward_params)                               # "ED" "ICU" "Medical" "Surgical"

transitions <- as.matrix(read.csv("data/transitions.csv", row.names = 1, check.names = FALSE))

stays       <- read.csv("data/stay_hours.csv")
stay_hours  <- split(stays$hours, stays$ward)                      # a list: real stay lengths for each ward

arrivals    <- read.csv("data/arrivals_by_hour.csv")               # average arrivals per weekday and hour
history     <- read.csv("data/daily_infections_history.csv")       # E. coli infections per day, past 2 years


# ---------------------------------------------------------------------
# 3. LIVE DATA HOOKS
# ---------------------------------------------------------------------

# Who is in the hospital right now: their ward and hours left in that ward.
# The live feed has no ward-level data, so this uses a snapshot of the synthetic hospital.
get_current_state <- function() {
  read.csv("data/start_state.csv")
}

# Current time in Brisbane (same clock as the website)
brisbane_now <- function() {
  as.POSIXct(format(Sys.time(), tz = "Australia/Brisbane"), tz = "Australia/Brisbane")
}

# Expected ED arrivals for each of the next 72 hours
arrival_rates <- function(live, start_time) {
  times   <- start_time + 3600 * (0:(HOURS - 1))       # the 72 forecast hours
  hour    <- as.integer(format(times, "%H"))
  weekday <- weekdays(times)

  # (a) daily level: the synthetic hospital's average arrivals for that weekday
  per_day <- tapply(arrivals$mean_arrivals, arrivals$weekday, sum)[weekday]

  # (b) shape across the day: taken from the website's hourly average (the purple chart),
  #     because the synthetic arrivals are flat across the day
  shape <- rep(1, 24)
  ratio <- 1
  if (!is.null(live)) {
    avg_by_hour <- live$busyness$avg_waiting
    if (mean(avg_by_hour) > 0) shape <- avg_by_hour / mean(avg_by_hour)

    # (c) live signal: how busy the ED is right now compared with usual for this hour
    usual_now <- avg_by_hour[live$live$hour + 1]
    if (usual_now > 0) ratio <- live$live$patients_waiting / usual_now
    ratio <- min(max(ratio, 0.5), 2)                    # keep it between half and double
  }

  rate <- per_day / 24 * shape[hour + 1]
  rate[1:LIVE_HOURS] <- rate[1:LIVE_HOURS] * ratio      # live signal only for the next few hours
  list(rate = as.numeric(rate), ratio = ratio, times = times)
}


# ---------------------------------------------------------------------
# 4. THE DIGITAL TWIN ENGINE  (hour-by-hour agent-based simulation)
# ---------------------------------------------------------------------
# Same infection rule as the notebook:
#   hourly chance   p = 1 - exp(-lambda * crowd_multiplier / hours_today)
#   where hours_today = hours the patient spends in that ward on that day, and
#   Stage 1: lambda = exp(Intercept + b * device utilisation of the ward), no crowding
#   Stage 2: Stage 1 + crowding (chance x crowd_factor in crowded hours, scaled by k)
#   Stage 3: lambda from each patient's own device, age 65+ and sex, + crowding
# crowd_factor is only changed by crowding_sensitivity.R; the dashboard always uses the default (2)
simulate <- function(rates, times, stage = 1, arrivals_mult = 1, icu_mult = 1, device_mult = 1,
                     crowd_factor = crowd["crowd_factor"], track = FALSE) {
  # track = TRUE also records every patient, every hour, in the first N_TRACK runs (for the animation)
  frames <- list()                                          # one table per run and hour: who is where
  gone   <- list()                                          # discharged patients
  set.seed(1)    # same random numbers every time -> baseline and scenario are compared fairly

  # What-if changes
  trans <- transitions
  trans[, "ICU"] <- trans[, "ICU"] * icu_mult                    # more / fewer ICU transfers
  trans <- trans / rowSums(trans)                                # rows add up to 1 again
  device_prob <- pmin(ward_params$device_prob * device_mult, 1)  # more / less device use
  names(device_prob) <- WARDS

  hour_of_day <- as.integer(format(times, "%H"))
  start <- get_current_state()

  # Results: one row per run, one column per hour
  new_inf <- matrix(0, N_RUNS, HOURS)     # simulated new infections
  exp_inf <- matrix(0, N_RUNS, HOURS)     # expected new infections (sum of everyone's chance)
  occupancy <- array(0, c(N_RUNS, HOURS, length(WARDS)), dimnames = list(NULL, NULL, WARDS))
  occ_male  <- occupancy                  # male patients per ward
  occ_age65 <- occupancy                  # patients aged 65+ per ward
  exp_male  <- matrix(0, N_RUNS, HOURS)   # expected new infections in male patients
  exp_age65 <- matrix(0, N_RUNS, HOURS)   # expected new infections in patients aged 65+

  for (run in 1:N_RUNS) {

    # ---- the patients in hospital at the start ----
    ward     <- start$ward                                   # which ward each patient is in
    left     <- start$hours_left                             # hours left in that ward
    device   <- rbinom(length(ward), 1, device_prob[ward])   # 1 = has a device
    age65    <- start$age65                                  # 1 = aged 65 or over
    male     <- start$male                                   # 1 = male
    infected <- rep(FALSE, length(ward))
    id       <- seq_along(ward)                              # a number for each patient (for the animation)
    next_id  <- length(ward)

    for (h in 1:HOURS) {
      hours_left_in_day <- 24 - hour_of_day[h]

      # hours each patient is in their ward today
      if (h == 1) {
        # patients already in hospital: we assume they have been in this ward since midnight
        today <- pmin(pmax(ceiling(left), 1), hours_left_in_day) + hour_of_day[1]
      }
      if (h > 1 && hour_of_day[h] == 0) {
        # a new day starts at midnight
        today <- pmin(pmax(ceiling(left), 1), 24)
      }

      # ---- (a) new patients arrive in ED ----
      n_new <- rpois(1, rates[h] * arrivals_mult)
      if (n_new > 0) {
        new_left <- sample(stay_hours$ED, n_new, replace = TRUE)
        ward     <- c(ward, rep("ED", n_new))
        left     <- c(left, new_left)
        device   <- c(device, rbinom(n_new, 1, device_prob["ED"]))
        age65    <- c(age65, rbinom(n_new, 1, arrival_profile["age65"]))
        male     <- c(male,  rbinom(n_new, 1, arrival_profile["male"]))
        infected <- c(infected, rep(FALSE, n_new))
        today    <- c(today, pmin(pmax(ceiling(new_left), 1), hours_left_in_day))
        id       <- c(id, next_id + 1:n_new)
        next_id  <- next_id + n_new
      }

      # ---- (b) infection: each uninfected patient may get infected this hour ----
      patients_in_ward <- table(factor(ward, levels = WARDS))

      if (stage < 3) {
        # Stages 1-2: everyone in a ward shares the ward's rate (from its device utilisation)
        device_util <- tapply(device, factor(ward, levels = WARDS), mean)    # share with a device, per ward
        device_util[is.na(device_util)] <- 0
        ward_lambda <- exp(coef["Intercept"] + coef["device_utilisation"] * device_util)
        names(ward_lambda) <- WARDS
        lambda <- ward_lambda[ward]
      } else {
        # Stage 3: each patient's own rate (own device, age 65+, sex)
        lambda <- exp(coef3["Intercept"] + coef3["own_device"] * device +
                      coef3["age65"] * age65 + coef3["male"] * male)
      }

      # Stages 2-3: higher chance in crowded hours, scaled by k so the overall rate stays the same
      crowd_mult <- 1
      if (stage >= 2) {
        crowded <- patients_in_ward[ward] > ward_params[ward, "crowd_threshold"]
        crowd_mult <- crowd_k(crowd_factor) * ifelse(crowded, crowd_factor, 1)
      }

      p <- 1 - exp(-lambda * crowd_mult / today)                             # each patient's hourly chance

      got_infected <- !infected & (runif(length(ward)) < p)
      exp_inf[run, h] <- sum(p[!infected])
      exp_male[run, h]  <- sum(p[!infected & male == 1])
      exp_age65[run, h] <- sum(p[!infected & age65 == 1])
      infected <- infected | got_infected
      new_inf[run, h] <- sum(got_infected)
      occupancy[run, h, ] <- patients_in_ward
      occ_male[run, h, ]  <- tapply(male,  factor(ward, levels = WARDS), sum, default = 0)
      occ_age65[run, h, ] <- tapply(age65, factor(ward, levels = WARDS), sum, default = 0)

      # animation: remember where everyone is this hour
      if (track && run <= N_TRACK) {
        frames[[length(frames) + 1]] <- data.frame(run = run, hour = h, id = id, ward = ward,
                                                   infected = infected, device = device, new = got_infected,
                                                   age65 = age65, male = male)
      }

      # ---- (c) one hour passes: patients whose ward stay is over move on, or go home ----
      left <- left - 1
      for (i in which(left <= 0)) {
        next_ward <- sample(colnames(trans), 1, prob = trans[ward[i], ])
        ward[i] <- next_ward
        if (next_ward != "DISCHARGE") {
          left[i]   <- sample(stay_hours[[next_ward]], 1)
          device[i] <- rbinom(1, 1, device_prob[next_ward])
          today[i]  <- min(max(ceiling(left[i]), 1), max(hours_left_in_day - 1, 1))
        }
      }

      # remove discharged patients
      stay <- ward != "DISCHARGE"
      if (track && run <= N_TRACK && any(!stay)) {
        gone[[length(gone) + 1]] <- data.frame(run = run, id = id[!stay], infected = infected[!stay], hour = h,
                                               age65 = age65[!stay], male = male[!stay])
      }
      ward <- ward[stay]; left <- left[stay]; device <- device[stay]
      age65 <- age65[stay]; male <- male[stay]; id <- id[stay]
      infected <- infected[stay]; today <- today[stay]
    }
  }
  list(new_inf = new_inf, exp_inf = exp_inf, occupancy = occupancy,
       occ_male = occ_male, occ_age65 = occ_age65, exp_male = exp_male, exp_age65 = exp_age65,
       track = list(frames = do.call(rbind, frames), gone = do.call(rbind, gone)))
}

# Summary for plots. Infections are rare, so:
#   line = expected number (mean of everyone's risk), band = 2.5% and 97.5% quantiles of simulated counts
summarise_cumulative <- function(sim) {
  cum_sim <- t(apply(sim$new_inf, 1, cumsum))      # running total, for each run
  cum_exp <- t(apply(sim$exp_inf, 1, cumsum))
  data.frame(mean = colMeans(cum_exp),
             low  = apply(cum_sim, 2, quantile, 0.025),
             high = apply(cum_sim, 2, quantile, 0.975))
}



# ---------------------------------------------------------------------
# HOSPITAL ANIMATION  (NetLogo-style view of one simulated run)
# ---------------------------------------------------------------------
# Where each ward is drawn on the screen (x0-x1 = left-right, y0-y1 = bottom-top)
BOXES <- data.frame(
  ward = c("ED", "ICU", "Medical", "Surgical", "Discharged"),
  x0   = c(0,    5,     5,         5,          9.5),
  x1   = c(4,    8,     8,         8,          13),
  y0   = c(0,    3.4,   1.7,       0,          0),
  y1   = c(4,    4.6,   2.9,       1.2,        4)
)

make_animation <- function(track, times, run_no, shape_by = "Device") {
  hosp <- track$frames[track$frames$run == run_no, ]      # the chosen example run
  gone <- track$gone[track$gone$run == run_no, ]

  # every hour: patients in hospital + everyone discharged so far
  all_frames <- lapply(1:HOURS, function(h) {
    in_hosp <- hosp[hosp$hour == h, ]
    out <- gone[gone$hour < h, ]
    rbind(
      data.frame(hour = h, id = in_hosp$id, ward = in_hosp$ward,
                 infected = in_hosp$infected | in_hosp$new, device = in_hosp$device,
                 age65 = in_hosp$age65, male = in_hosp$male),
      data.frame(hour = rep(h, nrow(out)), id = out$id, ward = rep("Discharged", nrow(out)),
                 infected = out$infected, device = rep(0, nrow(out)),
                 age65 = out$age65, male = out$male)
    )
  })
  df <- do.call(rbind, all_frames)

  # position inside the ward box: a fixed spot for each patient number, so dots only move when the ward changes
  box <- BOXES[match(df$ward, BOXES$ward), ]
  df$x <- box$x0 + 0.2 + ((df$id * 0.618034) %% 1) * (box$x1 - box$x0 - 0.4)
  df$y <- box$y0 + 0.2 + ((df$id * 0.754878) %% 1) * (box$y1 - box$y0 - 0.4)
  df$status <- ifelse(df$infected, "Infected", "Not infected")
  df$device <- ifelse(df$device == 1, "Device", "No device")
  # dot shape: device, sex or age group (chosen above the animation)
  if (shape_by == "Sex") {
    df$shape <- ifelse(df$male == 1, "Male", "Female")
    symbols  <- c("Male" = "square", "Female" = "circle")
  } else if (shape_by == "Age group") {
    df$shape <- ifelse(df$age65 == 1, "Aged 65+", "Under 65")
    symbols  <- c("Aged 65+" = "triangle-up", "Under 65" = "circle")
  } else {
    df$shape <- df$device
    symbols  <- c("Device" = "diamond", "No device" = "circle")
  }
  df$hover  <- paste0("Patient ", df$id, "<br>", df$ward, "<br>",
                      ifelse(df$male == 1, "Male", "Female"), ", ", ifelse(df$age65 == 1, "aged 65+", "under 65"), "<br>",
                      df$status, ", ", df$device)

  # plotly needs every colour/shape group in every frame, so add 4 invisible points per hour
  dummy <- expand.grid(hour = 1:HOURS, status = c("Infected", "Not infected"), shape = names(symbols),
                       stringsAsFactors = FALSE)
  dummy <- cbind(dummy, id = -(1:nrow(dummy)), ward = "", infected = NA, device = "", age65 = NA, male = NA,
                 x = -5, y = -5, hover = "")
  df <- rbind(df, dummy[, names(df)])

  # ward titles with patient counts, and a heading with the time and infections so far
  labels <- do.call(rbind, lapply(1:HOURS, function(h) {
    d <- df[df$hour == h & df$id > 0, ]
    n <- table(factor(d$ward, levels = BOXES$ward))
    rbind(
      data.frame(hour = h, x = BOXES$x0, y = BOXES$y1 + 0.18,
                 text = paste0("<b>", BOXES$ward, "</b> (", as.numeric(n), ")")),
      data.frame(hour = h, x = 4.5, y = 5.2,
                 text = paste0("<b>", format(times[h], "%a %d %b, %H:00"), "</b>   |   infected so far: ",
                               sum(d$infected)))
    )
  }))

  # grey rectangles for the wards
  shapes <- lapply(1:nrow(BOXES), function(i) list(
    type = "rect", x0 = BOXES$x0[i], x1 = BOXES$x1[i], y0 = BOXES$y0[i], y1 = BOXES$y1[i],
    fillcolor = "rgba(200,200,200,0.15)", line = list(color = "grey")))

  plot_ly(df, x = ~x, y = ~y, frame = ~hour, ids = ~id, type = "scatter", mode = "markers",
          color = ~status, colors = c("Infected" = "#d62728", "Not infected" = "#2ca02c"),
          symbol = ~shape, symbols = symbols,
          marker = list(size = 9, line = list(width = 0.5, color = "white")),
          text = ~hover, hoverinfo = "text") |>
    add_text(data = labels, x = ~x, y = ~y, text = ~text, frame = ~hour, textposition = "middle right",
             inherit = FALSE, showlegend = FALSE, hoverinfo = "skip") |>
    layout(shapes = shapes,
           xaxis = list(range = c(-0.3, 13.3), visible = FALSE),
           yaxis = list(range = c(-0.3, 5.5), visible = FALSE),
           legend = list(orientation = "h", y = -0.05)) |>
    animation_opts(frame = 400, transition = 250, redraw = FALSE) |>
    animation_slider(currentvalue = list(prefix = "Hour ")) |>
    animation_button(label = "Play")
}

# ---------------------------------------------------------------------
# 5. USER INTERFACE  (layout based on the REDCAST dashboard)
# ---------------------------------------------------------------------
about_page <- fluidRow(column(8, offset = 2,
  h2("Project overview"),
  p("This dashboard is a digital twin of a hospital. It forecasts hospital-acquired ", em("E. coli"),
    " infections treated with Ceftriaxone for the next 3 days, hour by hour and ward by ward, and lets the user test what-if scenarios (changes to ED arrivals, ICU transfers and device use). Other ",
    em("E. coli"), " infections (treated with other antibiotics) are not included."),
  p("This web-based app was built as part of an honours thesis in Biomedical Engineering (The University of Sydney), completed in collaboration with CSIRO: ",
    em("Agent-based digital twins for antimicrobial resistance: a case study of hospital-acquired ", em("E. coli"), " infections.")),

  h3("Why ", em("E. coli"), " treated with Ceftriaxone?"),
  p("The clinical importance of this project reflects ", em("E. coli"), " resistant to third-generation cephalosporins (such as Ceftriaxone) being one of the two antimicrobial resistance indicators the World Health Organization tracks under Sustainable Development Goal 3.d.2. The model therefore forecasts this one drug-bug combination rather than all infections."),

  h3("What it forecasts (and what it does not)"),
  p("The data contain no susceptibility (lab) results, so resistance cannot be confirmed. The twin forecasts how many infections caused by this antimicrobial-resistance priority organism and drug combination are expected, in which ward, over the next 72 hours. It does not yet forecast resistant infections specifically (see Future work)."),

  h3("The data"),
  tags$ul(
    tags$li(strong("Synthetic hospital: "), "six synthetic data files simulated from real hospital data (admissions, patients, ward movements, infections, medications and devices), covering January 2023 to December 2024. No real patient data are used."),
    tags$li(strong("Size: "), "20,000 admissions from 12,000 patients, 33,110 ward movements across 4 wards (ED, ICU, Medical, Surgical), 2,000 infections and 10,032 device episodes (urinary catheters and central lines)."),
    tags$li(strong("Live feed: "), "Redcliffe Hospital Emergency Department, from ",
            tags$a(href = "https://openhospitals.health.qld.gov.au/facility/redcliffe-hospital", target = "_blank", "Queensland Health Open Hospitals"),
            ". The website updates every 15-30 minutes.")
  ),

  h3("How the model was built"),
  tags$ol(
    tags$li(strong("Poisson regression: "), "infections per ward per day were modelled with Poisson regression, using the number of admissions present in the ward that day (n_at_risk) as the exposure and the ward's device utilisation as the predictor. Data were split by calendar date (70% train, 15% validation, 15% test) and the model was chosen on the validation set. For ", em("E. coli"), " treated with Ceftriaxone, Model 2 (device utilisation only) was kept: each one-unit rise in device utilisation multiplies the infection rate by about 6 (IRR 5.97, p = 0.052)."),
    tags$li(strong("Agent-based model: "), "every admission is an agent that moves through the wards hour by hour. Each agent's daily infection rate from the Poisson model is spread evenly over the hours it spends in the ward that day."),
    tags$li(strong("Validation: "), "the agent-based model replayed the 112 test days with the real patients and their real ward movements. In one validation run of 1,000 simulations, the median was 23 infections (95% range 14 to 32) against 27 actual infections, for all three stages. Repeating the run gives slightly different values.")
  ),

  h3("Model stages (choose one in the sidebar)"),
  tags$ul(
    tags$li(strong("Stage 1 - Baseline: "), "a patient's chance depends on their ward's device utilisation (Poisson Model 2). Everyone in the same ward on the same day has the same chance."),
    tags$li(strong("Stage 2 - Crowding: "), "Stage 1, plus a higher chance in hours when a ward has more patients than usual (above its 90th percentile), as a proxy for transmission. The factor (2) is assumed, because the synthetic data shows no crowding effect; the overall rate is kept the same."),
    tags$li(strong("Stage 3 - Agent profiles: "), "each patient's own device, age (65+) and sex set their chance (patient-level Poisson model), plus Stage 2 crowding.")
  ),

  h3("How the digital twin runs"),
  tags$ol(
    tags$li("Loads the synthetic hospital's settings: infection rates, crowding thresholds, ward transitions, stay lengths and arrival patterns."),
    tags$li("Places 57 admissions from a snapshot of the synthetic hospital (noon, 16 December 2024) in their wards, with their recorded age, sex and hours left in that ward. Each gets a random device draw and starts uninfected."),
    tags$li("Reads the live Redcliffe ED data: the number of patients waiting now and the usual number waiting for each hour of today."),
    tags$li("Patients waiting now divided by the usual number for this hour (limited to 0.5-2) tells the twin whether the ED is busier or quieter than normal. This scales new arrivals for the next 6 hours only."),
    tags$li("Starting at the next full hour, simulates 72 hours in one-hour steps: new patients arrive in ED, patients move between wards when their stay ends, and infections are simulated."),
    tags$li("Simulates 1,000 possible 72-hour scenarios, each with different random outcomes."),
    tags$li("Every 15 minutes (or when 'Refresh live data' is clicked) the forecast is recalculated from the same starting snapshot (always the same 57 patients), using the current time as the start and the newest live reading. The twin does not advance in time."),
    tags$li("If the website cannot be read, the app falls back to normal arrival patterns and shows a message that live data is unavailable.")
  ),

  h3("Reading the outputs"),
  tags$ul(
    tags$li(strong("Infection forecast: "), "the line is the expected number of infections (the average over 1,000 runs of the sum of every agent's hourly infection probability). The average is used because infections are rare and the median simulated count is usually zero. Bands show the 2.5% to 97.5% range of the simulated counts."),
    tags$li(strong("Ward occupancy: "), "median number of patients per ward, with 95% uncertainty bands."),
    tags$li(strong("Patient profile: "), "the share of patients who are male and aged 65+ in each ward, and the expected infections split by sex and age group. Sex and age only change a patient's chance in Stage 3."),
    tags$li(strong("What-if: "), "the scenario is rerun with the same random numbers as the baseline, so the difference is due to the change, not chance. Results vary slightly at each refresh."),
    tags$li(strong("Animation: "), "the NetLogo-style view shows one example run only, not the forecast. Hover over a patient to see their sex and age group; the dot shape can show device use, sex or age group.")
  ),

  h3("Colours"),
  tags$ul(
    tags$li(span(style = "color:grey", "● Grey"), " - historical records"),
    tags$li(span(style = "color:red", "● Red"), " - historical average"),
    tags$li(span(style = "color:green", "● Green"), " - current (live) value"),
    tags$li(span(style = "color:blue", "● Blue"), " - forecast (with 95% range)")
  ),

  h3("Key terms"),
  tags$ul(
    tags$li(strong("n_at_risk: "), "the number of admissions present in a ward at any time on a calendar day."),
    tags$li(strong("Device utilisation: "), "device days divided by n_at_risk for a ward on a day; every agent in that ward on that day has the same value."),
    tags$li(strong("Daily rate (λ): "), "the expected number of infections per admission per day from the Poisson model."),
    tags$li(strong("Hourly chance (p): "), "p = 1 − exp(−λ × crowd multiplier ÷ hours in the ward that day)."),
    tags$li(strong("Crowded hour: "), "an hour when a ward has more patients than its 90th percentile over two years (ED 42, ICU 6, Medical 14, Surgical 13). It reflects occupancy, not ED capacity or arrivals."),
    tags$li(strong("Scaling factor (k = 0.897): "), "keeps the hospital's total infection risk unchanged after crowding is added, so crowded hours get 1.79× and normal hours 0.90× the Stage 1 chance.")
  ),

  h3("Metric definitions (Queensland Health)"),
  tags$ul(
    tags$li(strong("Patients currently waiting: "), "triaged and waiting for treatment to commence."),
    tags$li(strong("Median waiting time: "), "time from triage to treatment within which 50% of patients started treatment, over the last 2 hours."),
    tags$li(strong("Treatment spaces: "), "physical capacity of the treatment area, not real-time availability."),
    tags$li(strong("Hourly average (purple chart on the website): "), "patients waiting in the same hour, same weekday and same month, averaged over the past 2 years.")
  ),

  h3("Assumptions"),
  tags$ul(
    tags$li("The ", em("E. coli"), " infections treated with Ceftriaxone were considered hospital-acquired and assigned to the ward the patient was in when the infection was detected. The usual surveillance rule is onset 48 hours or more after admission, but the data do not record onset time or whether the infection was present on admission."),
    tags$li(HTML("Specific to one drug-bug combination (<em>E. coli</em> treated with Ceftriaxone); other infections and medications are treated the same as in non-infected patients.")),
    tags$li("Each admission is an independent agent; readmissions are not linked."),
    tags$li("The chance of infection is constant within a day. Stage 1 has no person-to-person transmission; Stage 2 uses crowding as a proxy to simulate possible patient-to-patient transmission."),
    tags$li("The crowding effect (factor 2) is assumed, not estimated from the data, because the synthetic data showed no crowding effect."),
    tags$li("For clinicians, the ED is not technically a ward. In this prototype (built without clinician input) it is treated as one; future work will incorporate clinician input and terminology."),
    tags$li("A patient can be infected at most once per admission; infection duration is not modelled."),
    tags$li("Stages 1-2: all patients have identical susceptibility. Stage 3: susceptibility depends only on own device, age 65+ and sex."),
    tags$li("Length of stay is unrelated to patient profile and infection. Principal diagnoses (ICD-10) are not linked to infections and are not used."),
    tags$li("Patients already in a ward at the start are assumed to have been there since midnight, which gives them a realistic share of the day's risk."),
    tags$li("Device use for each agent is drawn at random from historical shares, because a real live feed does not report device use."),
    tags$li("No snapshot patient is infected at the start (6 of the 57 had an infection before the snapshot, none with ", em("E. coli"), " treated with Ceftriaxone)."),
    tags$li("The live feed only scales ED arrivals; ward moves, stay lengths and device use come from the synthetic data."),
    tags$li("The data are synthetic, so weak associations with ward and device use are expected.")
  ),

  h3("Limitations"),
  tags$ul(
    tags$li("Redcliffe Hospital is larger and has more than 4 wards, so live data is used only as a relative signal (busier or quieter than usual), not as raw counts."),
    tags$li("In about a third of admissions the last ward record ends before discharge. 312 infections detected in this gap are assigned to the patient's last ward, which raises the Medical and Surgical rates; a sensitivity analysis was carried out."),
    tags$li("Spreading each daily rate over the hours present causes slightly more simulated infections around midnight; this is an artefact, not a real time-of-day effect."),
    tags$li("A real crowding / transmission effect needs real data with contact or occupancy information; for this reason, a data-based crowding effect was not estimated in this work (Stage 2 uses an assumed factor).")
  ),

  h3("Future work"),
  tags$ul(
    tags$li("Add antibiotic susceptibility lab results (susceptible / intermediate / resistant) from hospital pathology so the twin forecasts resistant infections specifically."),
    tags$li("Add real prescribing data, so prior antibiotic use (the main real-world driver of resistance) raises an agent's chance of acquiring a resistant organism."),
    tags$li("Future work will incorporate clinical input to design antibiotic stewardship scenarios, for example: what if first-line Ceftriaxone use dropped by 20%? How many resistant infections and isolation rooms would that save?")
  ),

  h3("Contact and feedback"),
  p("Built by Chahat Masih (BME Honours, The University of Sydney; in collaboration with CSIRO Biosecurity), supervised by Dr Sandhya Clement and Dr Marcela Cespedes."),
  p("Questions, ideas or feedback on this dashboard are very welcome: ",
    tags$a(href = "mailto:cmas0826@uni.sydney.edu.au?subject=HAI%20Digital%20Twin%20feedback", "cmas0826@uni.sydney.edu.au")),

  h3("AI use statement"),
  p("AI assisted with debugging, refining and structuring the code. AI was not used to conceive the project, design the study, choose the methods or interpret the results; these, and all modelling decisions and validation, are the work of Chahat Masih (The University of Sydney).")
))

ui <- navbarPage("HAI Digital Twin", selected = "Digital Twin",
  tabPanel("About", about_page),
  tabPanel("Digital Twin", sidebarLayout(

    sidebarPanel(width = 3,
      h4("Live: Redcliffe ED"),
      uiOutput("live_box"),
      actionButton("refresh", "Refresh live data", icon = icon("rotate")),
      hr(),
      h4("Model stage"),
      radioButtons("stage", NULL, selected = 1, choices = c(
        "1 - Baseline (ward device use)" = 1,
        "2 - + Crowding" = 2,
        "3 - + Agent profiles" = 3)),
      hr(),
      h4("What-if scenario"),
      sliderInput("arr", "ED arrivals (% of normal)",   50, 150, 100, step = 10, post = "%"),
      sliderInput("icu", "ICU transfers (% of normal)",  0, 200, 100, step = 10, post = "%"),
      sliderInput("dev", "Device use (% of normal)",     0, 150, 100, step = 10, post = "%"),
      actionButton("run", "Run what-if", class = "btn-primary", icon = icon("play"))
    ),

    mainPanel(width = 9, tabsetPanel(
      tabPanel("Live ED data",
        plotlyOutput("live_plot", height = "420px"),
        helpText("Grey: every saved scrape for this weekday. Red: the website's 2-year hourly average. ",
                 "Orange: our own hourly average from the saved scrapes (should match red once enough data builds up). Green: now.")),
      tabPanel("Infection forecast",
        plotlyOutput("day_plot", height = "380px"),
        plotlyOutput("hour_plot", height = "320px")),
      tabPanel("Ward occupancy",
        plotlyOutput("occupancy_plot", height = "480px"),
        helpText("Median patients per ward (95% band). Dotted lines: each ward's usual busy level (90th percentile), used by Stage 2.")),
      tabPanel("Patient profile",
        h4("Who is in the hospital - next 72 hours (baseline)"),
        tableOutput("profile_table"),
        helpText("Average number of patients in each ward over the next 72 hours (all 1,000 runs), and the share who are male and aged 65+. ",
                 "Patients already in hospital keep their recorded sex and age; new arrivals are drawn from the synthetic hospital's history (49% male, 36% aged 65+)."),
        plotlyOutput("profile_inf_plot", height = "380px"),
        uiOutput("profile_note")),
      tabPanel("What-if simulation",
        uiOutput("whatif_text"),
        plotlyOutput("whatif_plot", height = "420px"),
        hr(),
        h4("Watch the hospital (one example run)"),
        fluidRow(
          column(3, radioButtons("anim_which", NULL, choices = c("Baseline", "Scenario"), selected = "Baseline", inline = TRUE)),
          column(5, radioButtons("anim_shape", "Dot shape shows:", choices = c("Device", "Sex", "Age group"),
                                 selected = "Device", inline = TRUE)),
          column(4, uiOutput("anim_run_choice"))
        ),
        plotlyOutput("anim_plot", height = "560px"),
        helpText(strong("Hover over a patient to see their sex and age group (under 65 / 65+), ward, infection status and device. "),
                 "Press Play. Each dot is a patient: green = not infected, red = infected. ",
                 "The dot shape shows device use, sex or age group, as chosen above (see the legend). ",
                 "Patients arrive in ED, move between wards and leave to Discharged. ",
                 "Infections are rare, so pick a run from the list to see different possible outcomes; ",
                 "the numbers above use all runs."))
    ))
  ))
)


# ---------------------------------------------------------------------
# 6. SERVER  (what happens when the app runs)
# ---------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- Live data: read on start, when the button is pressed, and every 15 minutes ----
  live <- reactiveVal(NULL)
  get_live <- function() {
    result <- tryCatch(scrape_redcliffe(), error = function(e) NULL)
    if (!is.null(result)) {
      save_scrape(result)      # adds to data/live_log.csv, so a history builds up
      live(result)
    }
  }
  observe({ invalidateLater(15 * 60 * 1000); get_live() })
  observeEvent(input$refresh, get_live())

  output$live_box <- renderUI({
    l <- live()
    if (is.null(l)) return(p(style = "color:#b00", "Live data unavailable - using normal arrival patterns."))
    tagList(
      p(strong(l$live$patients_waiting), " patients waiting"),
      p(strong(l$live$median_wait_min), " min median wait"),
      p(strong(l$live$treatment_spaces), " treatment spaces"),
      p(style = "font-size:85%; color:grey", "Website updated: ", l$live$website_updated)
    )
  })

  # ---- Baseline forecast (re-runs when the live data changes) ----
  forecast_input <- reactive({
    live()
    start_time <- trunc(brisbane_now(), "hours") + 3600          # next full hour
    arrival_rates(live(), start_time)
  })
  baseline <- reactive({
    f <- forecast_input()
    withProgress(message = "Simulating the next 72 hours...",
                 simulate(f$rate, f$times, stage = as.numeric(input$stage), track = TRUE))
  })

  # ---- What-if scenario (runs when the button is pressed) ----
  scenario <- eventReactive(input$run, {
    f <- forecast_input()
    result <- withProgress(message = "Simulating what-if...",
      simulate(f$rate, f$times, stage = as.numeric(input$stage),
               arrivals_mult = input$arr / 100, icu_mult = input$icu / 100, device_mult = input$dev / 100,
               track = TRUE))
    result$settings <- c(input$arr, input$icu, input$dev, as.numeric(input$stage))   # remember what was used
    result
  })

  # ---- Tab 1: live ED data (REDCAST style) ----
  output$live_plot <- renderPlotly({
    l <- live(); req(l)
    today <- l$live$weekday
    log <- if (file.exists("data/live_log.csv")) read.csv("data/live_log.csv") else l$live
    log <- log[log$weekday == today, ]
    our_average <- aggregate(patients_waiting ~ hour, data = log, FUN = mean)

    plot_ly() |>
      add_markers(data = log, x = ~hour, y = ~patients_waiting, name = "Historical (scraped)",
                  marker = list(color = "lightgrey", size = 7)) |>
      add_markers(x = 0:23, y = l$busyness$avg_waiting, name = "Average (website, 2 yrs)",
                  marker = list(color = "red", size = 11)) |>
      add_trace(data = our_average, x = ~hour, y = ~patients_waiting, name = "Our average (scraped)",
                type = "scatter", mode = "lines+markers",
                line = list(color = "orange", dash = "dot"), marker = list(color = "orange", size = 8)) |>      
      add_markers(x = l$live$hour, y = l$live$patients_waiting, name = "Current",
                  marker = list(color = "green", size = 14)) |>
      layout(title = paste("Redcliffe ED - patients waiting by hour,", today),
             xaxis = list(title = "Hour of day", dtick = 1),
             yaxis = list(title = "No. of waiting patients"),
             legend = list(orientation = "h", y = -0.2))
  })

  # ---- Tab 2a: infections per day, next 3 days vs history (REDCAST style) ----
  output$day_plot <- renderPlotly({
    sim   <- baseline()
    day   <- as.Date(forecast_input()$times, tz = "Australia/Brisbane")
    days  <- unique(day)
    label <- format(days, "%a %d %b")

    # forecast per day: expected number and 95% range of simulated counts
    forecast <- data.frame(x = label, mean = NA, low = NA, high = NA)
    for (i in seq_along(days)) {
      counts   <- rowSums(sim$new_inf[, day == days[i], drop = FALSE])
      expected <- rowSums(sim$exp_inf[, day == days[i], drop = FALSE])
      forecast[i, 2:4] <- c(mean(expected), quantile(counts, 0.025), quantile(counts, 0.975))
      if (sum(day == days[i]) < 24) forecast$x[i] <- paste(label[i], "(part day)")
    }

    # history: all past days with the same weekday
    past <- do.call(rbind, lapply(seq_along(days), function(i)
      data.frame(x = forecast$x[i], y = history$infections[history$weekday == weekdays(days[i])])))
    past_average <- aggregate(y ~ x, past, mean)

    plot_ly() |>
      add_markers(data = past, x = ~x, y = ~jitter(y, 0.6), name = "Historical (same weekday)",
                  marker = list(color = "lightgrey", size = 6), hoverinfo = "skip") |>
      add_markers(data = past_average, x = ~x, y = ~y, name = "Average",
                  marker = list(color = "red", size = 13)) |>
      add_markers(data = forecast, x = ~x, y = ~mean, name = "Forecast",
                  marker = list(color = "blue", size = 13),
                  error_y = list(type = "data", symmetric = FALSE, color = "blue",
                                 array = ~high - mean, arrayminus = ~mean - low),
                  text = ~sprintf("Expected: %.2f<br>95%% range: %.0f - %.0f", mean, low, high),
                  hoverinfo = "text") |>
      layout(title = "<i>E. coli</i> treated with Ceftriaxone: infections per day - next 3 days",
             xaxis = list(title = "", categoryorder = "array", categoryarray = forecast$x),
             yaxis = list(title = "Infections per day", rangemode = "tozero"),
             legend = list(orientation = "h", y = -0.15))
  })

  # ---- Tab 2b: running total of infections, hour by hour ----
  output$hour_plot <- renderPlotly({
    s <- summarise_cumulative(baseline())
    s$time <- forecast_input()$times
    plot_ly(s, x = ~time) |>
      add_ribbons(ymin = ~low, ymax = ~high, name = "95% range",
                  fillcolor = "rgba(0,0,255,0.15)", line = list(width = 0)) |>
      add_lines(y = ~mean, name = "Forecast (expected)", line = list(color = "blue")) |>
      layout(title = "<i>E. coli</i> treated with Ceftriaxone: cumulative infections, hour by hour", xaxis = list(title = ""),
             yaxis = list(title = "Cumulative infections"), legend = list(orientation = "h", y = -0.2))
  })

  # ---- Tab 3: patients per ward ----
  output$occupancy_plot <- renderPlotly({
    sim   <- baseline()
    times <- forecast_input()$times
    colours <- c(ED = "#1f77b4", ICU = "#d62728", Medical = "#2ca02c", Surgical = "#9467bd")
    fig <- plot_ly()
    for (w in WARDS) {
      patients <- sim$occupancy[, , w]                    # runs x hours
      fig <- fig |>
        add_ribbons(x = times, ymin = apply(patients, 2, quantile, 0.025), ymax = apply(patients, 2, quantile, 0.975),
                    fillcolor = adjustcolor(colours[w], 0.15), line = list(width = 0),
                    showlegend = FALSE, hoverinfo = "skip") |>
        add_lines(x = times, y = apply(patients, 2, median), name = w, line = list(color = colours[w])) |>
        add_lines(x = range(times), y = rep(ward_params[w, "crowd_threshold"], 2),
                  line = list(color = colours[w], dash = "dot"), showlegend = FALSE, hoverinfo = "skip")
    }
    fig |> layout(title = "Patients per ward - next 72 hours", yaxis = list(title = "Patients"),
                  legend = list(orientation = "h", y = -0.15))
  })

  # ---- Tab 3b: patient profile (sex and age group) ----
  output$profile_table <- renderTable({
    sim <- baseline()
    rows <- lapply(c(WARDS, "All wards"), function(w) {
      pick <- if (w == "All wards") WARDS else w
      total <- sum(sim$occupancy[, , pick]); n <- sum(sim$occupancy[, , pick]) / (N_RUNS * HOURS)
      data.frame(Ward = w,
                 `Average patients` = sprintf("%.1f", n),
                 `Male` = sprintf("%.0f%%", 100 * sum(sim$occ_male[, , pick]) / total),
                 `Female` = sprintf("%.0f%%", 100 * (1 - sum(sim$occ_male[, , pick]) / total)),
                 `Aged 65+` = sprintf("%.0f%%", 100 * sum(sim$occ_age65[, , pick]) / total),
                 `Under 65` = sprintf("%.0f%%", 100 * (1 - sum(sim$occ_age65[, , pick]) / total)),
                 check.names = FALSE)
    })
    do.call(rbind, rows)
  }, striped = TRUE, align = "lccccc")

  # expected infections in the next 72 h, split by sex and by age group
  output$profile_inf_plot <- renderPlotly({
    sim   <- baseline()
    total <- mean(rowSums(sim$exp_inf))
    male  <- mean(rowSums(sim$exp_male))
    old   <- mean(rowSums(sim$exp_age65))
    d <- data.frame(split = c("Sex", "Sex", "Age group", "Age group"),
                    group = factor(c("Male", "Female", "Aged 65+", "Under 65"),
                                   levels = c("Male", "Female", "Aged 65+", "Under 65")),
                    value = c(male, total - male, old, total - old))
    plot_ly(d, x = ~split, y = ~value, color = ~group, type = "bar",
            colors = c("Male" = "#1f77b4", "Female" = "#ff7f0e", "Aged 65+" = "#7f7f7f", "Under 65" = "#bcbd22"),
            text = ~sprintf("%s: %.2f", group, value), textposition = "inside", hoverinfo = "text") |>
      layout(barmode = "stack",
             title = "<i>E. coli</i> treated with Ceftriaxone: expected infections in the next 72 h, by sex and age group",
             xaxis = list(title = "", categoryorder = "array", categoryarray = c("Sex", "Age group")),
             yaxis = list(title = "Expected infections"), legend = list(orientation = "h", y = -0.15))
  })

  output$profile_note <- renderUI({
    if (as.numeric(input$stage) < 3) {
      helpText("Stages 1-2: sex and age do not change a patient's chance of infection, so the split above only reflects who is in hospital. ",
               "Choose Stage 3 in the sidebar to let sex and age affect each patient's chance.")
    } else {
      helpText("Stage 3: each patient's chance depends on their own device, sex and age. In the synthetic data, males had a 1.48 times higher ",
               "infection rate (95% CI 1.02-2.14, p = 0.037); being aged 65+ showed no clear effect (0.92 times, p = 0.66).")
    }
  })

  # ---- Tab 4: what-if results ----
  output$whatif_text <- renderUI({
    s <- scenario()
    b <- baseline()
    expected_base     <- mean(rowSums(b$exp_inf))
    expected_scenario <- mean(rowSums(s$exp_inf))
    range_text <- function(sim) {
      totals <- rowSums(sim$new_inf)
      sprintf("95%% range %.0f to %.0f", quantile(totals, 0.025), quantile(totals, 0.975))
    }
    wellPanel(
      h4(sprintf("Stage %d | scenario: arrivals %d%%, ICU transfers %d%%, device use %d%%",
                 s$settings[4], s$settings[1], s$settings[2], s$settings[3])),
      p("Expected infections in the next 72 h - baseline: ", strong(sprintf("%.2f", expected_base)),
        " (", range_text(b), ")  |  scenario: ", strong(sprintf("%.2f", expected_scenario)),
        " (", range_text(s), ")"),
      # positive = more infections than baseline, negative = infections avoided
      p(style = "font-size:120%", "Change in expected infections vs baseline: ",
        strong(sprintf("%+.2f (%+.0f%%)", expected_scenario - expected_base,
                       (expected_scenario / expected_base - 1) * 100))),
      p(style = "color:grey", em("E. coli"), " infections treated with Ceftriaxone are rare (about 1 every 4-5 days), so over 72 hours the expected numbers are small.")
    )
  })

  # ---- Hospital animation: scenario if it has been run, otherwise baseline ----
  anim_sim <- reactive({
    if (input$anim_which == "Scenario" && input$run > 0) scenario() else baseline()
  })

  # dropdown of example runs, showing how many infections each one had
  output$anim_run_choice <- renderUI({
    n_inf  <- rowSums(anim_sim()$new_inf)[1:N_TRACK]
    runs   <- setNames(1:N_TRACK, sprintf("Run %d (%d infection%s)", 1:N_TRACK, n_inf, ifelse(n_inf == 1, "", "s")))
    first  <- which(n_inf > 0)[1]                                # start on a run that has an infection
    selectInput("anim_run", NULL, choices = runs, selected = ifelse(is.na(first), 1, first))
  })

  output$anim_plot <- renderPlotly({
    req(input$anim_run)
    make_animation(anim_sim()$track, forecast_input()$times, as.numeric(input$anim_run), input$anim_shape)
  })

  output$whatif_plot <- renderPlotly({
    base  <- summarise_cumulative(baseline())
    scen  <- summarise_cumulative(scenario())
    times <- forecast_input()$times
    plot_ly(x = times) |>
      add_ribbons(ymin = base$low, ymax = base$high, name = "Baseline 95%",
                  fillcolor = "rgba(128,128,128,0.2)", line = list(width = 0)) |>
      add_lines(y = base$mean, name = "Baseline", line = list(color = "grey")) |>
      add_ribbons(ymin = scen$low, ymax = scen$high, name = "Scenario 95%",
                  fillcolor = "rgba(0,0,255,0.15)", line = list(width = 0)) |>
      add_lines(y = scen$mean, name = "Scenario", line = list(color = "blue")) |>
      layout(title = "<i>E. coli</i> treated with Ceftriaxone: cumulative infections, baseline vs what-if", xaxis = list(title = ""),
             yaxis = list(title = "Cumulative infections"), legend = list(orientation = "h", y = -0.2))
  })
}

shinyApp(ui, server)

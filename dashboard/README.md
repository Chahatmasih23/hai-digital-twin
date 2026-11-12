# Dashboard (R Shiny digital twin)

| File | Purpose |
|---|---|
| `app.R` | The dashboard. Tabs: About, Live ED, Infection forecast, Ward occupancy, Patient profile, What-if |
| `scrape_live.R` | Reads the Redcliffe Hospital ED page every 15 minutes (`scrape_redcliffe()`) and logs readings (`Rscript scrape_live.R`) |
| `crowding_sensitivity.R` | Repeats the what-if scenarios with crowding factor f = 1.5, 2 and 3 |
| `data/` | Parameter files exported by `notebooks/03_ABM_digital_twin.ipynb`, plus the live logs |

Run from this folder:

```r
install.packages(c("shiny", "plotly"))
shiny::runApp()
```

### Files in `data/`

| File | Contents |
|---|---|
| `coefficients.csv` | Stage 1 Poisson coefficients (*E. coli*/ceftriaxone, Model 2) |
| `stage2_crowding.csv` | Crowding factor f, rescaling constant k, data-estimated crowding IRR and crowded share |
| `stage3_coefficients.csv` | Stage 3 patient-level coefficients (own device, age ≥ 65, male) |
| `ward_params.csv` | Device probability and crowding threshold (90th percentile occupancy) per ward |
| `transitions.csv` | Ward-to-ward transfer probabilities |
| `stay_hours.csv` | Observed ward stay lengths (hours), resampled for new stays |
| `arrivals_by_hour.csv` | Mean ED arrivals by weekday and hour |
| `arrival_profile.csv` | Share of arrivals aged 65+ and male |
| `start_state.csv` | Starting snapshot (57 admissions at noon, 16 December 2024) |
| `daily_infections_history.csv` | Daily *E. coli*/ceftriaxone infections over the two years (for comparison plots) |
| `live_log.csv`, `busyness_log.csv` | Logged Redcliffe ED readings |
| `crowding_sensitivity.csv` | Output of `crowding_sensitivity.R` |

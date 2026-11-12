# HAI Digital Twin

An interpretable agent-based **digital twin** for forecasting hospital-acquired *Escherichia coli* infections in patients treated with ceftriaxone. The twin is validated on held-out data and driven by a live emergency department (ED) feed.

**Live dashboard:** https://chahat23masih.shinyapps.io/hai-digital-twin/

This repository holds the code behind the University of Sydney honours thesis *"Digital twin agent-based modelling for simulating and predicting antimicrobial resistance in hospital environments"* (Chahat Masih, 2026) and the related manuscript.

---

## How it works

The pipeline has four layers. Each layer is checked before the next one is built.

```
 Synthetic hospital records ──► Statistical layer ──► Agent-based layer ──► Digital twin ──► Users
 (2 years, 4 wards,             (Poisson regression     (hourly ABM, 3 stages,  (R Shiny dashboard,   (forecasts and
  20,000 admissions)             of ward-day infection   verified and validated  live Redcliffe ED     what-if scenarios)
                                 counts, temporal split) on 112 held-out days)    feed every 15 min)
```

| Step | What it does | Where |
|---|---|---|
| 1. Exploratory analysis | Data quality checks, cohort description, infection and device patterns | `notebooks/01_EDA.ipynb` |
| 2. Statistical layer | Builds the ward-day panel, fits three Poisson models (baseline, + device utilisation, + ward), selects on validation error, evaluates once on the test period, checks calibration and overdispersion, and repeats the analysis for *E. coli*/ceftriaxone | `notebooks/02_Poisson_regression.ipynb` |
| 3. Agent-based layer | Turns the Poisson rates into hourly infection risks for each admission. Three stages: baseline, crowding, agent profiles. Verifies against the regression, validates by 1,000 simulations of the test period, and exports the parameter files the dashboard uses | `notebooks/03_ABM_digital_twin.ipynb` |
| 4. Digital twin | Reads the live Redcliffe Hospital ED page, runs 1,000 simulations of the next 72 hours, and shows forecasts and what-if scenarios | `dashboard/` |

## Repository structure

```
hai-digital-twin/
├── notebooks/
│   ├── 01_EDA.ipynb                    exploratory data analysis
│   ├── 02_Poisson_regression.ipynb     statistical layer (Poisson regression)
│   └── 03_ABM_digital_twin.ipynb       agent-based model, validation, dashboard inputs
├── dashboard/
│   ├── app.R                           R Shiny digital twin dashboard
│   ├── scrape_live.R                   live scraper for the Redcliffe Hospital ED page
│   ├── crowding_sensitivity.R          what-if scenarios with crowding factor f = 1.5, 2, 3
│   └── data/                           model parameters exported by notebook 03, plus live logs
├── results/
│   ├── poisson/                        model comparison, coefficients, scores, pipeline flowchart
│   └── abm/                            ABM validation summary and figures
├── data/
│   ├── synthetic_hospital/             the six synthetic hospital CSV files
│   └── README.md                       data description
├── requirements.txt                    Python packages
├── CITATION.cff
└── LICENSE
```

## Getting started

### 1. Data
The synthetic hospital dataset is included in `data/synthetic_hospital/` (see [`data/README.md`](data/README.md)).

### 2. Run the notebooks (Python)
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
cd notebooks
jupyter lab
```
Run the notebooks in order (01 → 02 → 03), starting each one from inside `notebooks/`. Notebook 02 saves its tables and figures in `results/poisson/`. Notebook 03 reads those outputs and saves the dashboard's input files in `dashboard/data/`. The saved notebooks already contain their outputs, so you can read the results without running anything.

### 3. Run the dashboard (R)
```r
install.packages(c("shiny", "plotly"))
setwd("dashboard")
shiny::runApp()
```
The dashboard runs from the parameter files already in `dashboard/data/`, so it does not need the raw dataset. To record live readings on a schedule, run `Rscript dashboard/scrape_live.R`. To reproduce the crowding sensitivity table, run `Rscript crowding_sensitivity.R` from inside `dashboard/`.

## Key results (synthetic data)

| Layer | Result |
|---|---|
| Poisson (all organisms) | Device utilisation IRR 2.26 (95% CI 1.30–3.95); test error within 2% of the chance-only floor; no overdispersion |
| ABM validation (*E. coli*/ceftriaxone, 112 test days) | Median 23 infections (95% interval 14–32) vs 27 observed; all four wards inside their intervals |
| Live signal | One-hour-ahead error in ED demand reduced by 37% compared with the hourly average |
| What-if | Halving device use: −22% (Stages 1–2), −4% (Stage 3). 50% ED surge: +28% to +80%, depending on the assumed crowding effect |

**Important:** The data are fully synthetic and were not generated with intentional clinical associations. These results show that the pipeline works. They do not describe real hospitals. The model forecasts infections caused by an AMR priority organism, not laboratory-confirmed resistance.

## Data sources
- **Synthetic hospital records:** provided through a collaboration with CSIRO (see `data/README.md`).
- **Live ED data:** [Queensland Health Open Hospitals – Redcliffe Hospital](https://openhospitals.health.qld.gov.au/facility/redcliffe-hospital), public aggregate statistics ([about the data](https://www.health.qld.gov.au/open-hospitals-about-the-data)).

## Generative AI use
Generative AI (Claude, Anthropic) was used to help debug, refine and structure the code, and to organise this repository. It was not used to conceive the project, choose its aims, derive the methodology or implement the analyses. All modelling decisions, analyses and interpretations are the author's own.

## Citation
If you use this code, please cite it using the metadata in [`CITATION.cff`](CITATION.cff) (on GitHub, click **"Cite this repository"**).

## Authors and contact
**Chahat Masih**, School of Biomedical Engineering, The University of Sydney: cmas0826@uni.sydney.edu.au
Supervisors: Dr Sandhya Clement (The University of Sydney) and Dr Marcela Cespedes (CSIRO).

## License
Code is released under the MIT License (see [`LICENSE`](LICENSE)).

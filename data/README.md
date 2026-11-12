# Data

## Synthetic hospital dataset

The analyses use a fully synthetic hospital dataset provided through a collaboration with the Commonwealth Scientific and Industrial Research Organisation (CSIRO). No record corresponds to a real person. The six CSV files are included in `data/synthetic_hospital/`, which is where the notebooks expect them:

```
data/synthetic_hospital/
├── patients.csv           12,000 synthetic patients (sex, date of birth, age, postcode)
├── admissions.csv         20,000 admissions (admit and discharge times)
├── ward_movements.csv     33,110 ward stays (ED, ICU, Medical, Surgical)
├── device_exposure.csv    10,032 device episodes (urinary catheters, central lines)
├── infection_events.csv   2,000 infections (organism, time)
└── medications.csv        2,000 antibiotic courses
```

Coverage: 1 January 2023 to 30 December 2024; 9,746 of the patients have at least one admission. The data were not generated with intentional clinical associations between variables, so results describe the synthetic hospital, not real hospitals.

## Live emergency department data

`dashboard/scrape_live.R` reads the public Redcliffe Hospital ED page on Queensland Health's [Open Hospitals](https://openhospitals.health.qld.gov.au/facility/redcliffe-hospital) website. It saves readings to `dashboard/data/live_log.csv` (patients waiting, median wait, treatment spaces) and `dashboard/data/busyness_log.csv` (the two-year hourly average for the same weekday and month). These are aggregate public statistics with no patient-level information.

# SBE 19plus V2 — Universal CTD Processing & Caribbean Water Mass Analysis

**Research-oriented Python workflow** for Sea-Bird SBE 19plus / 19plus V2 CTD profile processing, quality control, and subsequent hydrographic water mass analysis.

**Author / maintainer:** Lorenzo Portillo Cogollo  
**Research interests:** physical oceanography, scientific data processing, reproducible marine research.  
**Language:** Python · Interface/documentation: Spanish/English

> **Development status:** research software. The two recovered scripts passed a Python syntax check. An end-to-end run using independently shareable real HEX/XMLCON datasets has **not** been performed for this public release. Review QC decisions and calibration records before scientific use.

## Scientific workflow

```text
SBE 19plus V2 .hex + matching .xmlcon + optional station log (.csv/.xlsx)
                       |
             [1] CTD processor v2.6
                 - sensor decoding / instrument coefficients
                 - cast detection, soak / surface equilibration
                 - adaptive loop editing and pressure-binned downcast
                 - physical variables and TEOS-10
                 - QC/QARTOD-style flags, exceptions & reporting
                 - static records handled separately
                       |
         CTD_BASE_FINAL.xlsx + QC reports + 300-dpi figures
                       |
            [2] Water mass analysis v3.4
                 - MLD / ILD / barrier layer thickness (BLT)
                 - hydrographic core detection (SUW, AAIW)
                 - theta-SA diagrams, sections, profiles
                 - local optimum multiparameter analysis (OMP)
                   with mSUW / mAAIW / mNADW reference types
                       |
        Analysis tables + figures + methods / diagnostic metadata
```

## Scope and limitations

- The raw-file decoder targets **SBE 19plus / 19plus V2, profiling mode, strain-gauge pressure, OutputFormat=0 HEX**. “Universal” means reusable across campaigns meeting these assumptions; it does **not** mean any CTD brand/model or HEX encoding.
- Requires the correct, instrument-specific **XMLCON calibration/configuration**. Do not reuse another instrument's coefficients.
- Additional channels, including oxygen, depend on actual sensor configuration and available calibration information. Automated oxygen timing adjustments **do not constitute oxygen calibration**.
- Horizontal comparisons depend on verified station time, latitude/longitude and data quality.
- OMP fractions are interpreted **locally at observed stations and depths**, not as percentages of the entire Caribbean Sea. SUW/AAIW core identification and OMP need profiles with sufficient depth and suitable properties; shallow casts cannot support deep-water assignments.
- Quality flags are aids to review, **not a certificate** that the measurements meet a formal oceanographic quality standard.
- Original HEX/XMLCON and institutional datasets are intentionally excluded from this repository.

## Installation

Python 3.10+ recommended. In the same environment used by Spyder or the terminal:

```bash
python -m pip install -r requirements.txt
```

Dependencies include GSW (TEOS-10) and the Sea-Bird Scientific Python processing toolkit. For desktop file selectors run in an interactive environment.

## Step 1 — Process raw CTD data

Open `src/ctd_processor_v2_6.py` in **Spyder** and run the entire file (F5), or in a terminal:

```bash
python src/ctd_processor_v2_6.py
```

The program asks for or obtains paths to the HEX folder, matching `.xmlcon`, station table when relevant, and a results folder. It offers **three input modes**:

1. Complete station log, including `archivo_hex` associations.
2. Station log with timestamps but missing HEX associations, using reviewable reconciliation.
3. HEX + XMLCON discovery, generating an inventory/template for a later confirmed processing pass.

See `examples/stations_template.csv` for the expected six key columns. The template filenames and positions are fictional, not usable input files.

By default **each batch run creates its own timestamped folder** inside `CTD_CORRIDAS_V26` to avoid mixing results across runs/campaigns. Profiles are distinguished from static series and air tests; the program preserves audit files and marks uncertain cases for manual review. The scientific output includes `CTD_BASE_FINAL.xlsx` and QC tables.

**Before each run:** validate the station log, timezone/clock offsets in `AUTO_CFG`, instrument configuration, and every flagged `REVISAR` cast. For time series classified as static, consult the separate static outputs rather than treating them as downcasts.

## Step 2 — Analyse water masses

Use **the resulting** `CTD_BASE_FINAL.xlsx` from Step 1:

```bash
python scripts/run_water_masses.py --input "path/to/CTD_BASE_FINAL.xlsx" --output "outputs/water_masses"
```

The analyzer is `src/water_masses_v3_4.py`; the small wrapper above supplies portable input/output paths. It writes OMP diagnostics, water-mass classification tables and publication-oriented figures. The method applies selected source-water types and thresholds based on **Torres et al. (2023)**; appropriate interpretation requires independent evaluation for each region/cruise.

This is a **two-stage workflow**. CTD decoding does not itself identify source-water fractions, and the OMP stage does not directly decode HEX files.

## Repository layout

```text
src/ctd_processor_v2_6.py     # Raw SBE CTD processing (v2.6-DISCOVERY-VF-ESTATICOS)
src/water_masses_v3_4.py      # Hydrography, MLD/ILD/BLT, water masses and OMP
scripts/run_water_masses.py   # Portable CLI for the second stage
examples/stations_template.csv
docs/METODOLOGIA.md
docs/PUBLICATION_CHECKLIST.md
tests/test_repository_structure.py
requirements.txt
.gitignore
```

## Reproducibility & privacy

- Keep archived original raw files **outside the public repository**; do not publish deployment coordinates or proprietary calibration/metadata without permission.
- The station template contains **synthetic metadata only**.
- Verify the contents of `git status` before every push; do not use `git add --force` for ignored data.
- Real-campaign verification, independent measurements, sensor serials, and calibration-date review are essential before release of scientific results.

## Reference

Torres, R. R., Latandret, S., Salon, J., & Dagua, C. (2023). Water masses in the Caribbean Sea and sub-annual variability in the Guajira upwelling region. *Ocean Dynamics, 73*, 39–57. https://doi.org/10.1007/s10236-022-01529-5

## License and reuse

License **not yet assigned**. Public visibility alone does not grant an open-source license. Confirm ownership/authorization and choose a license before inviting third-party reuse, especially for code developed in connection with institutional work.

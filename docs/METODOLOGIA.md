# Methodology and quality-control overview / Metodología

## Scope of the implementation

This repository preserves two recovered research scripts:
**CTD processor v2.6-DISCOVERY-VF-ESTATICOS** and **water-mass analysis v3.4**.

### 1. From HEX to oceanographic variables

- Read and inspect proprietary SBE-19plus/V2 profiling HEX data using the *corresponding* XMLCON sensor coefficients and acquisition settings (specified decoder format restrictions apply).
- Diagnose cast shape, initial immersion/equilibration, pressure loops and static/air-test cases.
- Calculate physically meaningful variables; use TEOS-10 for Absolute Salinity, Conservative Temperature, potential-density anomaly, depth, and related quantities where coordinates are available.
- Retain original recordings. Report applied transformations, decisions, reasons for review, derived QC flags and figures.
- Pressure-bin selected downcasts at the configured 1-dbar default. Flag rather than silently discard ambiguous events.

### 2. Hydrographic interpretation

- Identify MLD with sigma0 increment of 0.03 kg/m3 relative to 10 m; ILD with temperature threshold 0.2 °C relative to 10 m; compute BLT according to the implementation.
- Identify SUW subsurface salinity maxima (below the configurable threshold) and AAIW salinity minima where the vertical sampling supports them.
- Diagnose water-mass interpretation with theta–SA and oxygen information when trustworthy.
- Perform local three-reference-type OMP (mSUW, mAAIW, mNADW), inspecting mass residuals and the sensitivity to omission of oxygen.

These regional end-member choices **must not** be interpreted as a globally universal taxonomy. They are oceanographic assumptions specific to the analysis.

### 3. Validation status

Python source was syntax-checked when the repository was assembled. **No fully documented external example pair (HEX and XMLCON) was available to perform an end-to-end integration test in this release.** Users must test each campaign and validate instrument/sensor-dependent processing before publication.

### 4. References

- Torres, R. R., Latandret, S., Salon, J., & Dagua, C. (2023). *Water masses in the Caribbean Sea and sub-annual variability in the Guajira upwelling region*. Ocean Dynamics, 73, 39–57. https://doi.org/10.1007/s10236-022-01529-5
- TEOS-10 / Gibbs SeaWater: https://www.teos-10.org/
- Sea-Bird Scientific toolkit: https://github.com/Sea-BirdScientific/seabirdscientific
- IOOS QARTOD manuals: https://ioos.noaa.gov/project/qartod/

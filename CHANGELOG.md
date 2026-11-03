# Changelog

## Repository packaging — 2026-11
- Recovered original CTD v2.6 (`ctd_crucero_auto_v2_6_VF_directorios_estaticos.py`) with per-run folders and separated static data processing.
- Recovered hydrographic analysis v3.4 (`analisis_masas_agua_caribe_v3_4_publicacion.py`).
- Made the water-mass stage's default input path portable and added a CLI wrapper; scientific algorithms were not rewritten.
- Added a documented synthetic station-log template, dependency manifest, gitignore and preliminary syntax tests.
- No end-to-end physical/instrument validation has been performed for this packaging step.

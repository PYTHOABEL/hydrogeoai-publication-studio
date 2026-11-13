# HydroGeoAI Publication Studio v6

Scientific plotting **local prototype**. Gibbs removed from the interface and standalone ternary plots replaced with Durov-style composition fields. Retains Piper, Schoeller, PCA, HCA, quality control and correlation.

## Windows setup (Python 3.11 or 3.12)

Extract the ZIP. Open the extracted **hydrogeoai** folder, type `cmd` in File Explorer's address bar, press Enter. Then:

```bat
python -m pip install --default-timeout=1000 --retries=10 -r requirements.txt
python -m streamlit run app.py
```

Open http://localhost:8501. If another Streamlit server occupies 8501, stop it with Ctrl+C in its terminal.

## Features

- **Piper**: 20% interior grid, stronger facies region fills and boundary lines, labels on/off, legend position options.
- **Durov**: top cation triangle, left anion triangle, central nine colored compositional cells. Composition axes use normalized ion charge percentages, not raw mg/L.
- **Chadha**: quadrants based on group differences of meq%.
- **Stiff**: single-sample and comparative overlays, cations left and anions right, meq/L.
- **Schoeller**, **PCA** (scores/loadings/biplot), **HCA** Ward, **Pearson/Spearman correlation**, **facies/hardness**, **CBE**.
- **Theme presets**: Journal, Thesis, Conference, Consultancy, Grayscale. User-selectable fonts, bold, line/axis weights, markers, colors, legend inside/outside, SVG/PDF/PNG exports.
- **PCA labels**: per-sample and per-variable visibility, X/Y offsets and rotation in an editable table. This is **not direct dragging**.
- **Plot preview**: mouse zoom, pan, reset in a Plotly-rendered *raster image*. Exact print-quality figures should be exported to SVG/PDF; not pixel screenshots.
- Download complete outputs and figure archive.

## Mandatory scientific caveats

1. **Durov zones:** Colored nine cells indicate a descriptive 3 × 3 grouping of normalized major-ion compositions. Do **not** interpret them as validated hydrogeochemical processes, fixed reaction pathways, or canonical published Durov classification without specifying and verifying the convention.
2. **Piper fields:** Colors are compositional screening fields, not a replacement for quantified charge-balance or a mineral saturation model. Separating facies groups is convention-dependent.
3. **TDS:** Gibbs has been removed; the `TDS_est_mgL` calculated in chemistry.py is simply the sum of included ions, not measured TDS.
4. **PCA:** Standardized concentration variables; arrows are rescaled for visualization. Verify loadings, choice of variables, sample size and compositional bias before drawing geochemical conclusions.
5. **Stiff:** The same scale is applied to all samples in a dataset for meaningful comparison.
6. **Upload:** Only unqualified numeric ion concentrations, in mg/L. Zero is assigned to absent optional CO3/NO3 *columns*. Detection-limited values like `<0.01` are rejected. **Do not deploy publicly with sensitive customer data without security, privacy controls and independent validation.**
7. Fonts not installed on the Windows computer can be substituted by Matplotlib. Use a local installation of the selected font; no fonts are distributed.

## Status

This version has no accounts, subscription billing, project database, confidential uploads handling, or regulatory validation. It is a functional local research prototype for hands-on evaluation, not yet a ready commercial SaaS.

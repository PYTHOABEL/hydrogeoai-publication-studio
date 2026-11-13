# v6 release notes — 13 November 2026

- Removed Gibbs UI and separate ternary tab.
- Added Durov-style ternary-square visualization, Chadha quadrant diagram and Stiff polygon workflows.
- Enhanced Piper compositional coloring with explicit field border lines.
- Added publication-oriented theme presets, bold fonts, adjustable legends and vector exports.
- Added a raster preview with browser mouse zoom/pan/reset backed by Plotly. The vector export remains Matplotlib-based.
- Added row-by-row PCA label positions, visibility and rotation with `st.data_editor`.
- Preserved v5.1 chemistry, CBE, PCA, HCA, correlation, facies/hardness, and Schoeller.
- Corrected warning language: facies polygons do not establish mineral reactions or water evolution mechanisms.

**Testing:** automated numerical and figure-export tests are included. Full Windows browser smoke test and independent validation against established geochemical software are still needed.

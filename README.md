# Israel Terrain Lab

Static web app (GitHub Pages ready): 3D map of Israel with real elevation, plus terrain-analysis tools that run in the browser.

## Data
Elevation = **Copernicus GLO-30** DSM (ESA, TanDEM-X 2011–2015, 30 m, release 2022_1). It is the most accurate open, near-global elevation model (absolute vertical accuracy about ±2–4 m); it is a surface model, so tree tops and buildings are included. Israeli national 5 m DTMs are not openly downloadable. The public AWS bucket has no browser (CORS) access, so `tools/build_tiles.py` cuts it once into Web-Mercator Terrarium PNG tiles (`tiles/7..12`, about 110 MB, z12 ≈ 32 m/pixel) which are served from the same site.

Rebuild: `python tools/build_tiles.py` (needs numpy, scipy, pillow, requests, tifffile, imagecodecs).

### DTM (bare earth), default
`python tools/build_tiles.py fabdem` builds `tiles_dtm/` from **FABDEM V1-2** (Hawker et al. 2022, University of Bristol / Fathom): Copernicus GLO-30 with trees and buildings removed by machine learning, 30 m. Licence **CC BY-NC-SA 4.0**: non-commercial use only, attribution required, derivatives must carry the same licence. Commercial use needs a licence from the FABDEM authors. The Layers tab switches between DTM (default) and DSM; if `tiles_dtm/` is absent the app falls back to DSM.

### Higher resolution?
No open 10 m model exists for Israel: Copernicus EEA-10 is Europe only; Survey of Israel's national DTM (25 m) is restricted; Mapterhorn's "Israel 10 m" layer is derived from ASTER GDEM (30 m native) and is not more accurate.

## Tools
- **Pick point** – WGS84 + Israeli Transverse Mercator (ITM) coordinates and elevation; copy button; draggable.
- **Measure** – air distance, distance along the terrain, climb and descent.
- **Hill line** – two points: air distance, bearing, elevation profile, the highest ground between them, its height above the straight A–B line (vertical and perpendicular), and whether the line of sight is clear (with earth curvature and refraction k = 0.13, observer/target heights).
- **Line of sight** – viewshed from a point (radius up to 30 km): visible area green, hidden red, % and km².
- **Focus area** – circle (coordinate + radius) or polygon; analyses run only inside it.
- **Find terrain features** – filter by elevation, slope and facing direction (highlighted on the map), contour lines, peaks list, area statistics.
- **Layers** – satellite / topo / streets / dark, 3D terrain with exaggeration and pitch, hillshade, per-layer toggles.
- **Projects** – new/open/delete, autosave in the browser, export/import JSON file.

## Run locally
`python -m http.server 8765` in this folder, open http://localhost:8765 (the page needs a web server, not `file://`).

# Israel Terrain Lab

Static web app (GitHub Pages ready): 3D map of Israel with real elevation, plus terrain-analysis tools that run in the browser.

## Data
Elevation = **Copernicus GLO-30** DSM (ESA, TanDEM-X 2011–2015, 30 m, release 2022_1). It is the most accurate open, near-global elevation model (absolute vertical accuracy about ±2–4 m); it is a surface model, so tree tops and buildings are included. Israeli national 5 m DTMs are not openly downloadable. The public AWS bucket has no browser (CORS) access, so `tools/build_tiles.py` cuts it once into Web-Mercator Terrarium PNG tiles (`tiles/7..12`, about 110 MB, z12 ≈ 32 m/pixel) which are served from the same site.

Rebuild: `python tools/build_tiles.py` (needs numpy, scipy, pillow, requests, tifffile, imagecodecs).

### Elevation models (Layers tab)
- **DTM – FABDEM V1-2** (default, bare earth, 30 m): `python tools/build_tiles.py fabdem` builds `tiles_dtm/` (Hawker et al. 2022, Univ. Bristol / Fathom; Copernicus GLO-30 with trees and buildings removed by machine learning). Licence **CC BY-NC-SA 4.0**: non-commercial use only, attribution required, share-alike (see `tiles_dtm/LICENSE.txt`). Commercial use needs a licence from the FABDEM authors.
- **DSM – Copernicus GLO-30** (surface, 30 m): `tiles/`.
- **Mapterhorn** and **AWS Terrain Tiles (SRTM/ASTER mix)**: streamed online, no files needed; handy as independent cross-checks. Models whose tiles are not installed are greyed out.

### Higher resolution?
No open 10 m model exists for Israel: Copernicus EEA-10 is Europe only; Survey of Israel's national DTM (25 m) is restricted; Mapterhorn's "Israel 10 m" layer is derived from ASTER GDEM (30 m native) and is not more accurate.

## Tools
- **Pick point** – WGS84 + Israeli Transverse Mercator (ITM) coordinates and elevation; copy button; draggable.
- **Measure** – air distance, distance along the terrain, climb and descent.
- **Hill line** – two points: air distance, bearing, elevation profile, the highest ground between them, its height above the straight A–B line (vertical and perpendicular), and whether the line of sight is clear (with earth curvature and refraction k = 0.13, observer/target heights).
- **Line of sight** – viewshed from a point (radius up to 30 km): visible area green, hidden red, % and km².
- **Focus area** – circle (coordinate + radius) or polygon; analyses run only inside it.
- **Find terrain features** – filter by elevation, slope and facing direction (highlighted on the map), contour lines, peaks list, area statistics.
- **Tool mode** – a tool opens a small bar on the map: ＋ Add creates, ✎ Move drags the tool's own points, ✕ finishes. With no tool active nothing on the map reacts, so panning never moves points. Undo/redo (Ctrl+Z) covers elements and the focus area; ✎ in the Results list renames points, lines, areas.
- **PDF export** – top-view map picture (whole Israel or the focus area) plus tables of points, lines, measurements and areas.
- **Layers** – satellite / topo / streets / dark, 3D terrain with exaggeration and pitch, hillshade, per-layer toggles.
- **Projects** – new/open/delete, autosave in the browser, export/import JSON file.

## Run locally
`python -m http.server 8765` in this folder, open http://localhost:8765 (the page needs a web server, not `file://`).

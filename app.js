/* Israel Terrain Lab – 3D terrain + analysis on Copernicus GLO-30 tiles (tiles/{z}/{x}/{y}.png, Terrarium encoding) */
'use strict';
const $ = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const RAD = Math.PI / 180, R_EARTH = 6371008.8, K_REF = 0.13;
const ZD = 12, N12 = 2 ** ZD * 256;                       // DEM analysis zoom (pixel ~ 32 m at 32 N)
const BOUNDS = [34.15, 29.35, 36.0, 33.5];
const fmt = (v, d = 0) => Number.isFinite(v) ? v.toLocaleString('en-US', { maximumFractionDigits: d, minimumFractionDigits: d }) : '–';
const km = m => m >= 10000 ? fmt(m / 1000, 1) + ' km' : m >= 1000 ? fmt(m / 1000, 2) + ' km' : fmt(m, 0) + ' m';

proj4.defs('EPSG:2039', '+proj=tmerc +lat_0=31.7343936111111 +lon_0=35.2045169444444 +k=1.0000067 +x_0=219529.584 +y_0=626907.39 +ellps=GRS80 +towgs84=-48,55,52,0,0,0,0 +units=m +no_defs');
const toITM = (lon, lat) => proj4('WGS84', 'EPSG:2039', [lon, lat]);
const fromITM = (e, n) => proj4('EPSG:2039', 'WGS84', [e, n]);

/* ---------- geodesy ---------- */
function hav(a, b) {
  const p1 = a[1] * RAD, p2 = b[1] * RAD, dp = p2 - p1, dl = (b[0] - a[0]) * RAD;
  const h = Math.sin(dp / 2) ** 2 + Math.cos(p1) * Math.cos(p2) * Math.sin(dl / 2) ** 2;
  return 2 * R_EARTH * Math.asin(Math.min(1, Math.sqrt(h)));
}
function bearing(a, b) {
  const p1 = a[1] * RAD, p2 = b[1] * RAD, dl = (b[0] - a[0]) * RAD;
  return (Math.atan2(Math.sin(dl) * Math.cos(p2), Math.cos(p1) * Math.sin(p2) - Math.sin(p1) * Math.cos(p2) * Math.cos(dl)) / RAD + 360) % 360;
}
function dest(p, brg, d) {
  const dr = d / R_EARTH, b = brg * RAD, p1 = p[1] * RAD, l1 = p[0] * RAD;
  const p2 = Math.asin(Math.sin(p1) * Math.cos(dr) + Math.cos(p1) * Math.sin(dr) * Math.cos(b));
  const l2 = l1 + Math.atan2(Math.sin(b) * Math.sin(dr) * Math.cos(p1), Math.cos(dr) - Math.sin(p1) * Math.sin(p2));
  return [((l2 / RAD + 540) % 360) - 180, p2 / RAD];
}
function slerp(a, b, f) {
  const p1 = a[1] * RAD, l1 = a[0] * RAD, p2 = b[1] * RAD, l2 = b[0] * RAD;
  const A = [Math.cos(p1) * Math.cos(l1), Math.cos(p1) * Math.sin(l1), Math.sin(p1)];
  const B = [Math.cos(p2) * Math.cos(l2), Math.cos(p2) * Math.sin(l2), Math.sin(p2)];
  const w = Math.acos(Math.min(1, A[0] * B[0] + A[1] * B[1] + A[2] * B[2]));
  if (w < 1e-9) return [a[0], a[1]];
  const s1 = Math.sin((1 - f) * w) / Math.sin(w), s2 = Math.sin(f * w) / Math.sin(w);
  const x = s1 * A[0] + s2 * B[0], y = s1 * A[1] + s2 * B[1], z = s1 * A[2] + s2 * B[2];
  return [Math.atan2(y, x) / RAD, Math.atan2(z, Math.hypot(x, y)) / RAD];
}
const lonlatToPx = (lon, lat) => { const s = Math.sin(lat * RAD); return [(lon + 180) / 360 * N12, (0.5 - Math.log((1 + s) / (1 - s)) / (4 * Math.PI)) * N12]; };
const circlePoly = (c, r, n = 96) => { const a = []; for (let i = 0; i <= n; i++) a.push(dest(c, i * 360 / n, r)); return a; };
function pointInPoly(x, y, ring) {
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const xi = ring[i][0], yi = ring[i][1], xj = ring[j][0], yj = ring[j][1];
    if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}
function parseCoord(str) {
  const n = str.replace(/[,;]/g, ' ').trim().split(/\s+/).map(Number);
  if (n.length < 2 || n.some(isNaN)) return null;
  let [a, b] = n;
  if (Math.abs(a) <= 180 && Math.abs(b) <= 180) {                 // WGS84, lat first (swap if clearly lon first)
    if (a > 34 && a < 36.5 && b > 29 && b < 34) [a, b] = [b, a];
    return { lat: a, lon: b };
  }
  if (b >= 1e6) b -= 1e6;                                          // old ITM northing offset
  const ll = fromITM(a, b); return { lat: ll[1], lon: ll[0] };
}

/* ---------- DEM (z12 Terrarium tiles decoded to float arrays) ---------- */
const DEM = {
  tiles: new Map(), loading: new Map(),
  key: (x, y) => x * 100000 + y,
  load(tx, ty) {
    const k = this.key(tx, ty);
    if (this.tiles.has(k)) return Promise.resolve();
    if (this.loading.has(k)) return this.loading.get(k);
    const p = fetch(`tiles/${ZD}/${tx}/${ty}.png`).then(r => { if (!r.ok) throw 0; return r.blob(); })
      .then(b => createImageBitmap(b, { premultiplyAlpha: 'none', colorSpaceConversion: 'none' }))
      .then(bm => {
        const c = document.createElement('canvas'); c.width = c.height = 256;
        const x = c.getContext('2d', { willReadFrequently: true }); x.drawImage(bm, 0, 0);
        const d = x.getImageData(0, 0, 256, 256).data, a = new Float32Array(65536);
        for (let i = 0; i < 65536; i++) a[i] = d[i * 4] * 256 + d[i * 4 + 1] + d[i * 4 + 2] / 256 - 32768;
        this.tiles.set(k, a);
      }).catch(() => { this.tiles.set(k, null); }).finally(() => this.loading.delete(k));
    this.loading.set(k, p); return p;
  },
  async ensure(keys) {
    const need = [...keys].filter(k => !this.tiles.has(k.x * 100000 + k.y));
    for (let i = 0; i < need.length; i += 24) await Promise.all(need.slice(i, i + 24).map(k => this.load(k.x, k.y)));
  },
  async ensureLL(pts) {                                            // list of [lon,lat]
    const m = new Map();
    for (const p of pts) { const [X, Y] = lonlatToPx(p[0], p[1]); const x = Math.floor(X / 256), y = Math.floor(Y / 256); m.set(this.key(x, y), { x, y }); }
    await this.ensure(m.values());
  },
  async ensureBBox(w, s, e, n) {
    const [x0, y1] = lonlatToPx(w, s), [x1, y0] = lonlatToPx(e, n), keys = [];
    for (let x = Math.floor(x0 / 256); x <= Math.floor(x1 / 256); x++) for (let y = Math.floor(y0 / 256); y <= Math.floor(y1 / 256); y++) keys.push({ x, y });
    await this.ensure(keys);
  },
  px(ix, iy) { const t = this.tiles.get(this.key(ix >> 8, iy >> 8)); return t ? t[(iy & 255) * 256 + (ix & 255)] : NaN; },
  at(lon, lat) {                                                    // bilinear, NaN if tile missing
    const [X, Y] = lonlatToPx(lon, lat), fx = X - 0.5, fy = Y - 0.5, ix = Math.floor(fx), iy = Math.floor(fy), tx = fx - ix, ty = fy - iy;
    const a = this.px(ix, iy), b = this.px(ix + 1, iy), c = this.px(ix, iy + 1), d = this.px(ix + 1, iy + 1);
    return (a * (1 - tx) + b * tx) * (1 - ty) + (c * (1 - tx) + d * tx) * ty;
  },
  clear() { this.tiles.clear(); }
};

/* ---------- state ---------- */
const DEF = () => ({
  v: 1, name: 'New project', view: { center: [35.0, 31.45], zoom: 7.6, pitch: 60, bearing: 0 },
  roi: null, items: [], viewshed: null,
  set: { base: 'sat', l3d: true, ex: 1.5, pitch: 60, lHs: true, hs: 0.5, lItems: true, lRoi: true, lView: true, lCont: true, lMask: true, lPeaks: true,
    hObs: 2, hTgt: 0, vsR: 8, curv: true, fEmin: '', fEmax: '', fSmin: '', fSmax: '', fAsp: '', fAspTol: 45, cInt: 20, pSep: 600 }
});
let S = DEF(), PID = null;
let map, tool = null, pending = [];
const markers = new Map();                 // item id -> [maplibre markers]
let analysis = null;                       // transient: {grid..., peaks}
let activeLine = null;
const uid = () => Math.random().toString(36).slice(2, 9);

/* ---------- map ---------- */
const hsPaint = k => ({ 'hillshade-exaggeration': k, 'hillshade-shadow-color': '#05080b', 'hillshade-highlight-color': '#ffffff', 'hillshade-accent-color': '#2a3a44', 'hillshade-illumination-direction': 315 });
function style() {
  const raster = (id, tiles, att, max = 19) => ({ [id]: { type: 'raster', tiles, tileSize: 256, maxzoom: max, attribution: att } });
  return {
    version: 8,
    sources: {
      ...raster('sat', ['https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}'], 'Imagery © Esri, Maxar, Earthstar Geographics', 18),
      ...raster('osm', ['https://tile.openstreetmap.org/{z}/{x}/{y}.png'], '© OpenStreetMap contributors'),
      ...raster('topo', ['https://tile.opentopomap.org/{z}/{x}/{y}.png'], '© OpenTopoMap (CC-BY-SA), © OpenStreetMap', 17),
      dem: { type: 'raster-dem', tiles: [new URL('tiles/{z}/{x}/{y}.png', location.href).href.replace(/%7B/g, '{').replace(/%7D/g, '}')], tileSize: 256, minzoom: 7, maxzoom: 12, encoding: 'terrarium', bounds: BOUNDS, attribution: 'Elevation: Copernicus GLO-30 © ESA' },
      hs: { type: 'raster-dem', tiles: [new URL('tiles/{z}/{x}/{y}.png', location.href).href.replace(/%7B/g, '{').replace(/%7D/g, '}')], tileSize: 256, minzoom: 7, maxzoom: 12, encoding: 'terrarium', bounds: BOUNDS },
    },
    layers: [
      { id: 'bg', type: 'background', paint: { 'background-color': '#0b1115' } },
      ...['sat', 'topo', 'osm'].map(id => ({ id: 'b-' + id, type: 'raster', source: id, layout: { visibility: 'none' } })),
      { id: 'hillshade', type: 'hillshade', source: 'hs', paint: hsPaint(0.5) },
    ]
  };
}
const E = () => ({ type: 'FeatureCollection', features: [] });
function initMap() {
  map = new maplibregl.Map({ container: 'map', style: style(), center: S.view.center, zoom: S.view.zoom, pitch: S.view.pitch, bearing: S.view.bearing, maxPitch: 85, doubleClickZoom: false, canvasContextAttributes: { antialias: true }, maxBounds: [[BOUNDS[0] - 0.6, BOUNDS[1] - 0.6], [BOUNDS[2] + 0.6, BOUNDS[3] + 0.6]] });
  map.addControl(new maplibregl.NavigationControl({ visualizePitch: true }), 'top-right');
  map.addControl(new maplibregl.ScaleControl({ unit: 'metric' }), 'bottom-right');
  map.addControl(new maplibregl.GeolocateControl({ trackUserLocation: false }), 'top-right');
  window.__map = map;
  map.on('load', () => {
    for (const id of ['roi', 'lines', 'pending', 'contours', 'peaks']) map.addSource(id, { type: 'geojson', data: E() });
    map.addSource('viewshed', { type: 'image', url: BLANK, coordinates: [[0, 1], [1, 1], [1, 0], [0, 0]] });
    map.addSource('mask', { type: 'image', url: BLANK, coordinates: [[0, 1], [1, 1], [1, 0], [0, 0]] });
    map.addLayer({ id: 'viewshed', type: 'raster', source: 'viewshed', paint: { 'raster-opacity': 0.75, 'raster-fade-duration': 0 } });
    map.addLayer({ id: 'mask', type: 'raster', source: 'mask', paint: { 'raster-opacity': 0.8, 'raster-fade-duration': 0 } });
    map.addLayer({ id: 'roi-fill', type: 'fill', source: 'roi', paint: { 'fill-color': '#40d6c2', 'fill-opacity': 0.07 } });
    map.addLayer({ id: 'roi-line', type: 'line', source: 'roi', paint: { 'line-color': '#40d6c2', 'line-width': 2, 'line-dasharray': [3, 2] } });
    map.addLayer({ id: 'contours', type: 'line', source: 'contours', paint: { 'line-color': ['case', ['==', ['get', 'idx'], 1], '#f3d27a', '#e8c46a'], 'line-width': ['case', ['==', ['get', 'idx'], 1], 1.5, 0.6], 'line-opacity': 0.8 } });
    map.addLayer({ id: 'lines-casing', type: 'line', source: 'lines', paint: { 'line-color': '#000', 'line-width': 5, 'line-opacity': 0.5 }, layout: { 'line-cap': 'round', 'line-join': 'round' } });
    map.addLayer({ id: 'lines', type: 'line', source: 'lines', paint: { 'line-color': ['get', 'color'], 'line-width': 2.5 }, layout: { 'line-cap': 'round', 'line-join': 'round' } });
    map.addLayer({ id: 'pending', type: 'line', source: 'pending', paint: { 'line-color': '#ffffff', 'line-width': 2, 'line-dasharray': [2, 2] } });
    map.addLayer({ id: 'peaks', type: 'circle', source: 'peaks', paint: { 'circle-radius': 5, 'circle-color': '#ef6a5b', 'circle-stroke-color': '#fff', 'circle-stroke-width': 1.5 } });
    applySettings(); renderAll();
    try { map.setSky({ 'sky-color': '#79a8d6', 'horizon-color': '#d7e4ee', 'fog-color': '#d7e4ee', 'sky-horizon-blend': 0.5, 'horizon-fog-blend': 0.6, 'fog-ground-blend': 0.2 }); } catch (e) { }
    window.__ready = true;
  });
  map.on('click', onMapClick);
  map.on('dblclick', e => { if (tool === 'measure' || tool === 'roi-poly') { e.preventDefault(); finish(); } });
  map.on('mousemove', onMove);
  map.on('moveend', () => { S.view = { center: map.getCenter().toArray(), zoom: map.getZoom(), pitch: map.getPitch(), bearing: map.getBearing() }; autosave(); });
  map.on('error', e => { if (e && e.error && !/tiles\//.test(String(e.error.message || e.error.url || ''))) console.warn(e.error); });
}
const BLANK = 'data:image/gif;base64,R0lGODlhAQABAIAAAAAAAP///yH5BAEAAAAALAAAAAABAAEAAAIBRAA7';

function applySettings() {
  const s = S.set; if (!map.getLayer('hillshade')) return;
  for (const id of ['sat', 'topo', 'osm']) map.setLayoutProperty('b-' + id, 'visibility', s.base === id ? 'visible' : 'none');
  map.setLayoutProperty('hillshade', 'visibility', s.lHs ? 'visible' : 'none');
  map.setPaintProperty('hillshade', 'hillshade-exaggeration', +s.hs);
  if (s.l3d) map.setTerrain({ source: 'dem', exaggeration: +s.ex }); else map.setTerrain(null);
  const vis = (id, on) => map.getLayer(id) && map.setLayoutProperty(id, 'visibility', on ? 'visible' : 'none');
  vis('lines', s.lItems); vis('lines-casing', s.lItems); vis('roi-fill', s.lRoi); vis('roi-line', s.lRoi); vis('viewshed', s.lView); vis('contours', s.lCont); vis('mask', s.lMask); vis('peaks', s.lPeaks);
  for (const ms of markers.values()) for (const m of ms) m.getElement().style.display = s.lItems ? '' : 'none';
}

/* ---------- ui helpers ---------- */
function toast(t, ms = 2800) { const d = document.createElement('div'); d.textContent = t; $('#toast').append(d); setTimeout(() => d.remove(), ms); }
const busy = on => $('#busy').hidden = !on;
const tick = () => new Promise(r => setTimeout(r, 15));
function setTool(t) {
  tool = tool === t ? null : t; pending = []; renderPending();
  $$('[data-tool]').forEach(b => b.classList.toggle('on', b.dataset.tool === tool));
  const H = { pick: 'Click the map to drop a point and read its coordinates and elevation.', measure: 'Click points along the path. Double-click or Enter to finish.', line: 'Click point A, then point B. The profile shows the highest ground between them.', view: 'Click the observer position. Green = visible from there.', 'roi-circle': 'Click the centre of the focus circle (set the radius in the Area tab).', 'roi-poly': 'Click the polygon corners. Double-click or Enter to close it.' };
  $('#toolHint').textContent = tool ? H[tool] : 'Choose a tool, then click the map. Esc cancels. Double-click or Enter finishes a measurement.';
  map.getCanvas().style.cursor = tool ? 'crosshair' : '';
}
function renderPending() {
  const multi = tool === 'measure' || tool === 'roi-poly';
  $('#fab').hidden = !(multi && pending.length >= 1);
  $('#fFinish').disabled = pending.length < (tool === 'roi-poly' ? 3 : 2);
  if (!map || !map.getSource('pending')) return;
  const f = []; if (pending.length > 1) f.push({ type: 'Feature', properties: {}, geometry: { type: 'LineString', coordinates: pending } });
  map.getSource('pending').setData({ type: 'FeatureCollection', features: f });
}

/* ---------- items ---------- */
function mkMarker(ll, cls, label, drag, ondrag) {
  const el = document.createElement('div'); el.className = 'pt ' + cls;
  if (label) { const l = document.createElement('div'); l.className = 'lb'; l.textContent = label; el.append(l); }
  const m = new maplibregl.Marker({ element: el, draggable: !!drag }).setLngLat(ll).addTo(map);
  if (drag) m.on('dragend', () => ondrag(m.getLngLat().toArray()));
  return m;
}
function clearMarkers(id) { (markers.get(id) || []).forEach(m => m.remove()); markers.delete(id); }
function itemMarkers(it) {
  clearMarkers(it.id); const ms = [];
  if (it.type === 'point') ms.push(mkMarker(it.ll, '', it.label, true, ll => { it.ll = ll; update(it); }));
  if (it.type === 'measure') it.pts.forEach((p, i) => ms.push(mkMarker(p, 'm' + (i === it.pts.length - 1 ? '' : ' s'), i === it.pts.length - 1 ? it.label : '', true, ll => { it.pts[i] = ll; update(it); })));
  if (it.type === 'line') { ms.push(mkMarker(it.a, '', 'A', true, ll => { it.a = ll; update(it); })); ms.push(mkMarker(it.b, '', 'B', true, ll => { it.b = ll; update(it); })); }
  markers.set(it.id, ms);
  ms.forEach(m => m.getElement().style.display = S.set.lItems ? '' : 'none');
}
async function update(it) { await compute(it); itemMarkers(it); renderAll(false); if (it.type === 'line' && activeLine === it.id) showProfile(it); autosave(); }
async function compute(it) {
  if (it.type === 'point') { await DEM.ensureLL([it.ll]); it.elev = DEM.at(...it.ll); const itm = toITM(...it.ll); it.itm = itm; }
  if (it.type === 'measure') {
    await DEM.ensureLL(it.pts); let flat = 0, surf = 0, up = 0, down = 0;
    for (let i = 1; i < it.pts.length; i++) {
      const a = it.pts[i - 1], b = it.pts[i], d = hav(a, b); flat += d;
      const n = Math.max(1, Math.ceil(d / 30)); let prev = DEM.at(...a);
      for (let k = 1; k <= n; k++) { const q = slerp(a, b, k / n), z = DEM.at(...q), dh = z - prev, dd = d / n; surf += Math.hypot(dd, Number.isFinite(dh) ? dh : 0); if (dh > 0) up += dh; else if (dh < 0) down -= dh; prev = z; }
    }
    it.flat = flat; it.surf = surf; it.up = up; it.down = down; it.label = km(flat);
  }
  if (it.type === 'line') it.res = await lineProfile(it);
}
function lineProfile(it) {
  const { a, b } = it, D = hav(a, b), n = Math.min(4000, Math.max(8, Math.ceil(D / 20)));
  const pts = []; for (let i = 0; i <= n; i++) pts.push(slerp(a, b, i / n));
  return DEM.ensureLL(pts).then(() => {
    const hA = +S.set.hObs, hB = +S.set.hTgt, curv = S.set.curv, Reff = R_EARTH / (1 - K_REF);
    const d = new Float32Array(n + 1), z = new Float32Array(n + 1), eff = new Float32Array(n + 1);
    for (let i = 0; i <= n; i++) { d[i] = D * i / n; z[i] = DEM.at(pts[i][0], pts[i][1]); eff[i] = z[i] + (curv ? d[i] * (D - d[i]) / (2 * Reff) : 0); }
    const zA = z[0], zB = z[n], sA = zA + hA, sB = zB + hB;
    let maxZ = -1e9, iMax = 0, bestEx = -1e9, iEx = 0, minClr = 1e9, iClr = 0;
    for (let i = 0; i <= n; i++) {
      const f = i / n, chord = zA + (zB - zA) * f, sight = sA + (sB - sA) * f;
      if (z[i] > maxZ) { maxZ = z[i]; iMax = i; }
      if (i > 0 && i < n) { const ex = z[i] - chord; if (ex > bestEx) { bestEx = ex; iEx = i; } const c = sight - eff[i]; if (c < minClr) { minClr = c; iClr = i; } }
    }
    const slope = Math.hypot(D, zB - zA), perp = bestEx * D / slope;
    return { D, n, d, z, eff, zA, zB, sA, sB, maxZ, dMax: d[iMax], iMax, bestEx, iEx, dEx: d[iEx], zEx: z[iEx], perp, minClr, iClr, visible: minClr >= 0, pts, hA, hB, curv };
  });
}

/* ---------- map clicks ---------- */
async function onMapClick(e) {
  const ll = [e.lngLat.lng, e.lngLat.lat];
  document.body.classList.remove('open');
  if (!tool) return;
  if (tool === 'pick') { const it = { id: uid(), type: 'point', ll, label: 'P' + (S.items.filter(i => i.type === 'point').length + 1) }; S.items.push(it); await update(it); renderItems(); }
  else if (tool === 'measure') { pending.push(ll); renderPending(); }
  else if (tool === 'roi-poly') { pending.push(ll); renderPending(); }
  else if (tool === 'line') {
    pending.push(ll);
    if (pending.length === 2) { const it = { id: uid(), type: 'line', a: pending[0], b: pending[1] }; pending = []; renderPending(); S.items.push(it); activeLine = it.id; busy(true); await tick(); await update(it); busy(false); renderItems(); }
    else renderPending();
  }
  else if (tool === 'view') { S.viewshed = { ll }; await runViewshed(); }
  else if (tool === 'roi-circle') { S.roi = { type: 'circle', c: ll, r: +$('#roiR').value * 1000 }; await roiChanged(); setTool(null); }
  autosave();
}
async function finish() {
  while (pending.length > 1 && hav(pending[pending.length - 1], pending[pending.length - 2]) < 1) pending.pop();   // a double-click adds the last point twice
  if (tool === 'measure' && pending.length > 1) { const it = { id: uid(), type: 'measure', pts: pending.slice() }; pending = []; renderPending(); S.items.push(it); await update(it); renderItems(); setTool(null); }
  else if (tool === 'roi-poly' && pending.length > 2) { S.roi = { type: 'poly', pts: pending.slice() }; pending = []; renderPending(); await roiChanged(); setTool(null); }
  autosave();
}
let hoverBusy = false;
function onMove(e) {
  const ll = [e.lngLat.lng, e.lngLat.lat]; $('#roLL').textContent = `${ll[1].toFixed(5)}, ${ll[0].toFixed(5)}`;
  try { const t = toITM(...ll); $('#roITM').textContent = `ITM ${fmt(t[0])} / ${fmt(t[1])}`; } catch (_) { }
  if (hoverBusy) return; hoverBusy = true;
  DEM.ensureLL([ll]).then(() => { const z = DEM.at(...ll); $('#roEl').textContent = Number.isFinite(z) ? `${fmt(z, 1)} m` : '–'; }).finally(() => hoverBusy = false);
}

/* ---------- render lists / layers ---------- */
function renderAll(markersToo = true) {
  if (!map || !map.getSource('lines')) return;
  const f = [];
  for (const it of S.items) {
    if (it.type === 'measure') f.push({ type: 'Feature', properties: { color: '#f3b43f' }, geometry: { type: 'LineString', coordinates: it.pts.length > 1 ? densify(it.pts) : [] } });
    if (it.type === 'line') f.push({ type: 'Feature', properties: { color: it.res ? (it.res.visible ? '#5fd38a' : '#ef6a5b') : '#40d6c2' }, geometry: { type: 'LineString', coordinates: it.res ? it.res.pts.filter((_, i) => i % Math.ceil(it.res.n / 200) === 0 || i === it.res.n) : [it.a, it.b] } });
  }
  map.getSource('lines').setData({ type: 'FeatureCollection', features: f });
  const rf = []; if (S.roi) rf.push({ type: 'Feature', properties: {}, geometry: { type: 'Polygon', coordinates: [roiRing()] } });
  map.getSource('roi').setData({ type: 'FeatureCollection', features: rf });
  if (markersToo) { for (const id of [...markers.keys()]) clearMarkers(id); S.items.forEach(itemMarkers); }
  renderItems();
}
const densify = pts => { const o = []; for (let i = 1; i < pts.length; i++) { const n = Math.max(1, Math.ceil(hav(pts[i - 1], pts[i]) / 500)); for (let k = i === 1 ? 0 : 1; k <= n; k++) o.push(slerp(pts[i - 1], pts[i], k / n)); } return o; };
const roiRing = () => !S.roi ? [] : S.roi.type === 'circle' ? circlePoly(S.roi.c, S.roi.r) : [...S.roi.pts, S.roi.pts[0]];
function roiBBox() { const r = roiRing(); let w = 181, s = 91, e = -181, n = -91; for (const p of r) { w = Math.min(w, p[0]); e = Math.max(e, p[0]); s = Math.min(s, p[1]); n = Math.max(n, p[1]); } return [w, s, e, n]; }

function renderItems() {
  const box = $('#items'); box.innerHTML = '';
  const li = [...S.items].reverse();
  if (S.viewshed && S.viewshed.res) li.unshift({ id: 'vs', type: 'view' });
  for (const it of li) {
    const el = document.createElement('div'); el.className = 'item ' + (it.type === 'measure' ? 'm' : it.type === 'view' ? 'v' : '');
    let title = '', body = '';
    if (it.type === 'point') { title = it.label; body = `${it.ll[1].toFixed(5)}, ${it.ll[0].toFixed(5)}\nITM ${fmt(it.itm?.[0])} / ${fmt(it.itm?.[1])}\nElevation ${fmt(it.elev, 1)} m`; }
    if (it.type === 'measure') { title = 'Distance ' + km(it.flat || 0); body = `Air distance  ${fmt(it.flat, 1)} m\nAlong terrain ${fmt(it.surf, 1)} m\nClimb ${fmt(it.up, 0)} m · descent ${fmt(it.down, 0)} m\nBearing first leg ${fmt(bearing(it.pts[0], it.pts[1]), 1)}°`; }
    if (it.type === 'line' && it.res) {
      const r = it.res; title = `Hill line ${km(r.D)}`;
      body = `Air distance ${fmt(r.D, 1)} m · bearing ${fmt(bearing(it.a, it.b), 1)}°\nA ${fmt(r.zA, 1)} m  →  B ${fmt(r.zB, 1)} m\nHighest ground on line ${fmt(r.maxZ, 1)} m at ${km(r.dMax)}\nHeight above the A–B line ${fmt(r.bestEx, 1)} m (perpendicular ${fmt(r.perp, 1)} m)\nat ${km(r.dEx)}, ground ${fmt(r.zEx, 1)} m\nLine of sight ${r.visible ? '<span class="ok">CLEAR</span> (margin ' + fmt(r.minClr, 1) + ' m)' : '<span class="no">BLOCKED</span> (by ' + fmt(-r.minClr, 1) + ' m at ' + km(r.d[r.iClr]) + ')'}`;
    }
    if (it.type === 'view') { const v = S.viewshed.res; title = 'Line of sight area'; body = `Observer ${S.viewshed.ll[1].toFixed(5)}, ${S.viewshed.ll[0].toFixed(5)}\nGround ${fmt(v.z0, 1)} m + ${v.hObs} m\nRadius ${km(v.R)} · visible ${fmt(v.pct, 1)}% (${fmt(v.areaKm2, 2)} km²)`; }
    el.innerHTML = `<div class="t"><span>${title}</span><span></span></div><div class="d">${body}</div>`;
    const bt = el.querySelector('.t span:last-child');
    const add = (txt, fn, tip) => { const b = document.createElement('button'); b.className = 'mini'; b.textContent = txt; b.title = tip; b.onclick = fn; bt.append(b); };
    if (it.type === 'point') add('⧉', () => { navigator.clipboard?.writeText(`${it.ll[1].toFixed(6)}, ${it.ll[0].toFixed(6)}`); toast('Coordinates copied'); }, 'Copy lat, lon');
    if (it.type === 'line') add('▤', () => { activeLine = it.id; showProfile(it); }, 'Show profile');
    add('⌖', () => fly(it), 'Zoom to');
    add('✕', () => { if (it.type === 'view') { S.viewshed = null; setImage('viewshed', null); } else { S.items = S.items.filter(x => x !== it); clearMarkers(it.id); if (activeLine === it.id) $('#profile').hidden = true; } renderAll(false); autosave(); }, 'Delete');
    box.append(el);
  }
}
function fly(it) {
  let pts = it.type === 'point' ? [it.ll] : it.type === 'measure' ? it.pts : it.type === 'line' ? [it.a, it.b] : [S.viewshed.ll];
  if (pts.length === 1) return map.easeTo({ center: pts[0], zoom: Math.max(map.getZoom(), 12) });
  const b = pts.reduce((a, p) => [Math.min(a[0], p[0]), Math.min(a[1], p[1]), Math.max(a[2], p[0]), Math.max(a[3], p[1])], [181, 91, -181, -91]);
  map.fitBounds([[b[0], b[1]], [b[2], b[3]]], { padding: 80, pitch: map.getPitch() });
}

/* ---------- profile chart ---------- */
function showProfile(it) {
  const r = it.res; if (!r) return;
  $('#profile').hidden = false; $('#pTitle').textContent = 'Hill line profile';
  $('#pStats').textContent = `${fmt(r.D, 0)} m · top of hill ${fmt(r.bestEx, 1)} m above line · LOS ${r.visible ? 'clear' : 'blocked'}`;
  requestAnimationFrame(() => drawProfile(r));
}
let lastRes = null;
function drawProfile(r, hoverI) {
  lastRes = r;
  const cv = $('#pCanvas'), W = cv.clientWidth, H = cv.clientHeight, dpr = devicePixelRatio || 1;
  cv.width = W * dpr; cv.height = H * dpr; const g = cv.getContext('2d'); g.scale(dpr, dpr); g.clearRect(0, 0, W, H);
  const L = 52, Rm = 14, T = 14, B = 26, w = W - L - Rm, h = H - T - B;
  let lo = Infinity, hi = -Infinity; for (let i = 0; i <= r.n; i++) { lo = Math.min(lo, r.eff[i]); hi = Math.max(hi, r.eff[i]); }
  lo = Math.min(lo, r.sA, r.sB); hi = Math.max(hi, r.sA, r.sB); const pad = (hi - lo) * 0.12 + 5; lo -= pad; hi += pad;
  const X = d => L + d / r.D * w, Y = z => T + (1 - (z - lo) / (hi - lo)) * h;
  const css = getComputedStyle(document.documentElement);
  g.font = '11px JetBrains Mono, monospace'; g.fillStyle = '#8ba0ab'; g.strokeStyle = '#24343e'; g.lineWidth = 1;
  const step = niceStep((hi - lo) / 4);
  for (let v = Math.ceil(lo / step) * step; v < hi; v += step) { g.beginPath(); g.moveTo(L, Y(v)); g.lineTo(L + w, Y(v)); g.stroke(); g.textAlign = 'right'; g.fillText(fmt(v), L - 6, Y(v) + 4); }
  const dstep = niceStep(r.D / (W < 500 ? 3 : 6)); g.textAlign = 'center';
  for (let v = 0; v <= r.D; v += dstep) g.fillText(v >= 1000 ? fmt(v / 1000, v % 1000 ? 1 : 0) + ' km' : fmt(v) + ' m', X(v), H - 8);
  const gr = g.createLinearGradient(0, T, 0, T + h); gr.addColorStop(0, '#7a6a3a'); gr.addColorStop(1, '#2c3a30');
  g.beginPath(); g.moveTo(X(0), Y(lo)); for (let i = 0; i <= r.n; i += Math.max(1, Math.floor(r.n / (w * 1.5)))) g.lineTo(X(r.d[i]), Y(r.eff[i])); g.lineTo(X(r.D), Y(r.eff[r.n])); g.lineTo(X(r.D), Y(lo)); g.closePath(); g.fillStyle = gr; g.fill();
  g.strokeStyle = '#e8c46a'; g.lineWidth = 1.5; g.beginPath(); for (let i = 0; i <= r.n; i += Math.max(1, Math.floor(r.n / (w * 1.5)))) { const x = X(r.d[i]), y = Y(r.eff[i]); i ? g.lineTo(x, y) : g.moveTo(x, y); } g.stroke();
  g.setLineDash([5, 4]); g.strokeStyle = '#8ba0ab'; g.lineWidth = 1; g.beginPath(); g.moveTo(X(0), Y(r.zA)); g.lineTo(X(r.D), Y(r.zB)); g.stroke(); g.setLineDash([]);
  g.strokeStyle = r.visible ? '#5fd38a' : '#ef6a5b'; g.lineWidth = 2; g.beginPath(); g.moveTo(X(0), Y(r.sA)); g.lineTo(X(r.D), Y(r.sB)); g.stroke();
  const f = r.dEx / r.D, cz = r.zA + (r.zB - r.zA) * f;                // hill height marker (vertical above chord)
  g.strokeStyle = '#ef6a5b'; g.lineWidth = 2; g.beginPath(); g.moveTo(X(r.dEx), Y(cz)); g.lineTo(X(r.dEx), Y(r.eff[r.iEx])); g.stroke();
  g.fillStyle = '#ef6a5b'; g.beginPath(); g.arc(X(r.dEx), Y(r.eff[r.iEx]), 4, 0, 7); g.fill();
  const lab = `${fmt(r.bestEx, 1)} m above line · ${fmt(r.zEx, 0)} m asl`, tw = g.measureText(lab).width;
  g.fillStyle = '#dbe7ec'; g.textAlign = 'left'; g.fillText(lab, Math.max(4, Math.min(W - 4 - tw, X(r.dEx) - tw / 2)), Math.max(T + 10, Y(r.eff[r.iEx]) - 9));
  g.fillStyle = '#40d6c2'; g.beginPath(); g.arc(X(0), Y(r.sA), 4, 0, 7); g.arc(X(r.D), Y(r.sB), 4, 0, 7); g.fill();
  if (hoverI != null) { const x = X(r.d[hoverI]); g.strokeStyle = '#ffffff88'; g.lineWidth = 1; g.beginPath(); g.moveTo(x, T); g.lineTo(x, T + h); g.stroke(); }
  cv._geo = { L, w, r };
}
const niceStep = s => { const p = 10 ** Math.floor(Math.log10(s)), m = s / p; return (m < 1.5 ? 1 : m < 3.5 ? 2 : m < 7.5 ? 5 : 10) * p; };
let hoverMk = null;
$('#pCanvas').addEventListener('mousemove', ev => {
  const cv = ev.target, G = cv._geo; if (!G) return;
  const x = ev.offsetX, f = Math.min(1, Math.max(0, (x - G.L) / G.w)), i = Math.round(f * G.r.n), r = G.r;
  drawProfile(r, i); const tip = $('#pTip'); tip.hidden = false; tip.style.left = Math.min(x + 10, cv.clientWidth - 160) + 'px';
  tip.textContent = `${km(r.d[i])} · ${fmt(r.z[i], 1)} m`;
  if (!hoverMk) { const el = document.createElement('div'); el.className = 'pt s'; el.style.background = '#fff'; hoverMk = new maplibregl.Marker({ element: el }).setLngLat(r.pts[i]).addTo(map); } else hoverMk.setLngLat(r.pts[i]);
});
$('#pCanvas').addEventListener('mouseleave', () => { $('#pTip').hidden = true; if (lastRes) drawProfile(lastRes); if (hoverMk) { hoverMk.remove(); hoverMk = null; } });
window.addEventListener('resize', () => { if (!$('#profile').hidden && lastRes) drawProfile(lastRes); });
$('#pClose').onclick = () => $('#profile').hidden = true;

/* ---------- viewshed ---------- */
function setImage(id, canvas, bbox) {
  const src = map.getSource(id); if (!src) return;
  if (!canvas) { src.updateImage({ url: BLANK, coordinates: [[0, 1], [1, 1], [1, 0], [0, 0]] }); return; }
  src.updateImage({ url: canvas.toDataURL('image/png'), coordinates: [[bbox[0], bbox[3]], [bbox[2], bbox[3]], [bbox[2], bbox[1]], [bbox[0], bbox[1]]] });
}
async function runViewshed() {
  const v = S.viewshed; if (!v) return; busy(true); await tick();
  try {
    const R = Math.min(30, Math.max(0.5, +S.set.vsR)) * 1000, hObs = +S.set.hObs, hTgt = +S.set.hTgt, curv = S.set.curv, Reff = R_EARTH / (1 - K_REF);
    const [lon0, lat0] = v.ll, mLat = 111132, mLon = 111320 * Math.cos(lat0 * RAD);
    const bbox = [lon0 - R / mLon, lat0 - R / mLat, lon0 + R / mLon, lat0 + R / mLat];
    await DEM.ensureBBox(...bbox);
    const cell = Math.max(30, 2 * R / 1400), nPx = Math.ceil(2 * R / cell);
    const cv = document.createElement('canvas'); cv.width = cv.height = nPx; const g = cv.getContext('2d'), img = g.createImageData(nPx, nPx), d = img.data;
    const z0 = DEM.at(lon0, lat0), H0 = z0 + hObs, rays = Math.ceil(2 * Math.PI * (R / cell) * 1.6), step = cell * 0.6;
    let vis = 0, tot = 0; const seen = new Uint8Array(nPx * nPx);   // 1 = hidden, 2 = visible
    for (let k = 0; k < rays; k++) {
      const th = k / rays * 2 * Math.PI, sx = Math.sin(th), cy = Math.cos(th); let maxA = -1e9;
      for (let r = step; r <= R; r += step) {
        const lon = lon0 + r * sx / mLon, lat = lat0 + r * cy / mLat, z = DEM.at(lon, lat); if (!Number.isFinite(z)) continue;
        const drop = curv ? r * r / (2 * Reff) : 0, ze = z - drop, a = (ze + hTgt - H0) / r, show = a >= maxA;
        const col = Math.floor((r * sx + R) / cell), row = Math.floor((R - r * cy) / cell), idx = row * nPx + col;
        if (col >= 0 && col < nPx && row >= 0 && row < nPx) { if (show) seen[idx] = 2; else if (!seen[idx]) seen[idx] = 1; }
        const t = (ze - H0) / r; if (t > maxA) maxA = t;
      }
    }
    for (let i = 0; i < seen.length; i++) {
      if (seen[i] === 2) { d.set([64, 214, 110, 150], i * 4); vis++; tot++; } else if (seen[i] === 1) { d.set([200, 40, 40, 70], i * 4); tot++; }
    }
    g.putImageData(img, 0, 0);
    // fill any one-pixel gaps between rays (dilate visible into empty neighbours)
    const out = g.getImageData(0, 0, nPx, nPx), o = out.data, inpt = new Uint8ClampedArray(o);
    for (let y = 1; y < nPx - 1; y++) for (let x = 1; x < nPx - 1; x++) {
      const i = (y * nPx + x) * 4; if (inpt[i + 3]) continue; if ((x - nPx / 2) ** 2 + (y - nPx / 2) ** 2 > (nPx / 2) ** 2) continue;
      let gN = 0, rN = 0; for (const j of [i - 4, i + 4, i - nPx * 4, i + nPx * 4]) { if (inpt[j + 3] === 150) gN++; else if (inpt[j + 3] === 70) rN++; }
      if (gN >= 2) o.set([64, 214, 110, 150], i); else if (rN >= 2) o.set([200, 40, 40, 70], i);
    }
    g.putImageData(out, 0, 0);
    v.res = { z0, hObs, R, pct: tot ? 100 * vis / tot : 0, areaKm2: vis * cell * cell / 1e6 };
    v.bbox = bbox; setImage('viewshed', cv, bbox); renderItems();
    if (!markers.has('vs')) { const m = mkMarker(v.ll, 'v', 'observer', true, async ll => { v.ll = ll; await runViewshed(); autosave(); }); markers.set('vs', [m]); } else markers.get('vs')[0].setLngLat(v.ll);
  } finally { busy(false); }
}

/* ---------- focus area + analysis ---------- */
async function roiChanged() {
  renderAll(false); const b = roiBBox(); await DEM.ensureBBox(...b);
  const area = S.roi.type === 'circle' ? Math.PI * S.roi.r ** 2 : polyArea(S.roi.pts);
  $('#roiInfo').textContent = `${S.roi.type === 'circle' ? 'Circle r ' + km(S.roi.r) : 'Polygon ' + S.roi.pts.length + ' pts'} · ${fmt(area / 1e6, 2)} km²`;
  autosave();
}
const polyArea = pts => { const c = pts[0], k = [111320 * Math.cos(c[1] * RAD), 110574]; let a = 0; for (let i = 0; i < pts.length; i++) { const p = pts[i], q = pts[(i + 1) % pts.length]; a += (p[0] - c[0]) * k[0] * (q[1] - c[1]) * k[1] - (q[0] - c[0]) * k[0] * (p[1] - c[1]) * k[1]; } return Math.abs(a) / 2; };
async function runAnalysis() {
  if (!S.roi) { toast('Set a focus area first (Area tab).'); return; }
  busy(true); await tick();
  try {
    const bb = roiBBox(); await DEM.ensureBBox(...bb);
    const [w, s, e, n] = bb, mLat = 110574, mLon = 111320 * Math.cos((s + n) / 2 * RAD), wm = (e - w) * mLon, hm = (n - s) * mLat;
    const cs = Math.max(30, Math.sqrt(wm * hm) / 600), nx = Math.max(3, Math.ceil(wm / cs)), ny = Math.max(3, Math.ceil(hm / cs)), ring = roiRing();
    const Z = new Float32Array(nx * ny).fill(NaN), In = new Uint8Array(nx * ny);
    const lonAt = i => w + (i + 0.5) * cs / mLon, latAt = j => n - (j + 0.5) * cs / mLat;
    // include a 1-cell border outside of the ROI for slopes
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const lo = lonAt(i), la = latAt(j); Z[j * nx + i] = DEM.at(lo, la); In[j * nx + i] = pointInPoly(lo, la, ring) ? 1 : 0; }
    const sl = new Float32Array(nx * ny).fill(NaN), asp = new Float32Array(nx * ny).fill(NaN);
    for (let j = 1; j < ny - 1; j++) for (let i = 1; i < nx - 1; i++) {
      const k = j * nx + i, a = Z[k - nx - 1], b = Z[k - nx], c = Z[k - nx + 1], d = Z[k - 1], f = Z[k + 1], g = Z[k + nx - 1], h = Z[k + nx], ii = Z[k + nx + 1];
      const dzdx = ((c + 2 * f + ii) - (a + 2 * d + g)) / (8 * cs), dzdy = ((g + 2 * h + ii) - (a + 2 * b + c)) / (8 * cs);   // y grows south
      sl[k] = Math.atan(Math.hypot(dzdx, dzdy)) / RAD;
      asp[k] = (Math.atan2(-dzdx, dzdy) / RAD + 360) % 360;      // downslope compass direction (facing)
    }
    // stats
    let cnt = 0, mn = 1e9, mx = -1e9, sum = 0, ssum = 0, sc = 0;
    for (let k = 0; k < Z.length; k++) if (In[k] && Number.isFinite(Z[k])) { cnt++; sum += Z[k]; mn = Math.min(mn, Z[k]); mx = Math.max(mx, Z[k]); if (Number.isFinite(sl[k])) { ssum += sl[k]; sc++; } }
    // feature mask
    const num = id => { const v = $(id).value; return v === '' ? null : +v; };
    const eMin = num('#fEmin'), eMax = num('#fEmax'), sMin = num('#fSmin'), sMax = num('#fSmax'), aC = num('#fAsp'), aT = +$('#fAspTol').value || 45;
    const useF = [eMin, eMax, sMin, sMax, aC].some(x => x !== null);
    const cv = document.createElement('canvas'); cv.width = nx; cv.height = ny; const gx = cv.getContext('2d'), im = gx.createImageData(nx, ny); let mc = 0;
    if (useF) for (let k = 0; k < Z.length; k++) {
      if (!In[k] || !Number.isFinite(Z[k])) continue; let ok = true;
      if (eMin !== null && Z[k] < eMin) ok = false; if (eMax !== null && Z[k] > eMax) ok = false;
      if (sMin !== null && !(sl[k] >= sMin)) ok = false; if (sMax !== null && !(sl[k] <= sMax)) ok = false;
      if (aC !== null) { const dd = Math.abs(((asp[k] - aC + 540) % 360) - 180); if (!(dd <= aT) || !(sl[k] > 1)) ok = false; }
      if (ok) { im.data.set([64, 214, 194, 190], k * 4); mc++; }
    }
    gx.putImageData(im, 0, 0);
    const bbox2 = [w, n - ny * cs / mLat, w + nx * cs / mLon, n];
    setImage('mask', useF ? cv : null, bbox2);
    // contours
    const iv = Math.max(1, +$('#cInt').value || 20), feats = [], Zm = (i, j) => Z[j * nx + i];
    const lv0 = Math.ceil(mn / iv) * iv;
    for (let L = lv0; L <= mx; L += iv) {
      const idx = Math.round(L / iv) % 5 === 0 ? 1 : 0, segs = [];
      for (let j = 0; j < ny - 1; j++) for (let i = 0; i < nx - 1; i++) {
        if (!In[j * nx + i]) continue;
        const v = [Zm(i, j), Zm(i + 1, j), Zm(i + 1, j + 1), Zm(i, j + 1)]; if (v.some(x => !Number.isFinite(x))) continue;
        const c = (v[0] > L) * 8 + (v[1] > L) * 4 + (v[2] > L) * 2 + (v[3] > L); if (c === 0 || c === 15) continue;
        const pt = k => { const a = [[0, 0], [1, 0], [1, 1], [0, 1]], A = a[k], B = a[(k + 1) % 4], t = (L - v[k]) / (v[(k + 1) % 4] - v[k]); return [lonAt(i + A[0] + (B[0] - A[0]) * t), latAt(j + A[1] + (B[1] - A[1]) * t)]; };
        const T = { 1: [[3, 2]], 2: [[2, 1]], 3: [[3, 1]], 4: [[0, 1]], 5: [[0, 3], [2, 1]], 6: [[0, 2]], 7: [[0, 3]], 8: [[0, 3]], 9: [[0, 2]], 10: [[0, 1], [2, 3]], 11: [[0, 1]], 12: [[3, 1]], 13: [[2, 1]], 14: [[3, 2]] }[c];
        for (const [a, b] of T) segs.push([pt(a), pt(b)]);
      }
      if (segs.length) feats.push({ type: 'Feature', properties: { idx, level: L }, geometry: { type: 'MultiLineString', coordinates: segs } });
    }
    map.getSource('contours').setData({ type: 'FeatureCollection', features: feats });
    // peaks
    const sep = Math.max(100, +$('#pSep').value || 600), rw = Math.max(2, Math.round(sep / cs)), peaks = [];
    for (let j = rw; j < ny - rw; j++) for (let i = rw; i < nx - rw; i++) {
      const k = j * nx + i; if (!In[k]) continue; const z = Z[k]; if (!Number.isFinite(z)) continue; let top = true;
      for (let dj = -rw; dj <= rw && top; dj++) for (let di = -rw; di <= rw; di++) { const q = Z[(j + dj) * nx + i + di]; if (q > z || (q === z && (dj < 0 || (dj === 0 && di < 0)))) { top = false; break; } }
      if (top) peaks.push({ ll: [lonAt(i), latAt(j)], z });
    }
    peaks.sort((a, b) => b.z - a.z); const top = peaks.slice(0, 25);
    map.getSource('peaks').setData({ type: 'FeatureCollection', features: top.map(p => ({ type: 'Feature', properties: { z: p.z }, geometry: { type: 'Point', coordinates: p.ll } })) });
    $('#anInfo').textContent = `Grid ${nx}×${ny} (cell ${fmt(cs, 0)} m) · ${fmt(cnt * cs * cs / 1e6, 2)} km² analysed\nElevation min ${fmt(mn, 1)} · mean ${fmt(sum / cnt, 1)} · max ${fmt(mx, 1)} m\nMean slope ${fmt(ssum / sc, 1)}°` + (useF ? `\nMatching terrain: ${fmt(mc * cs * cs / 1e6, 3)} km² (${fmt(100 * mc / cnt, 1)}% of area)` : '\n(no filters set – only contours and peaks)') + `\nContours every ${iv} m (${feats.length} levels) · ${peaks.length} peaks`;
    const pb = $('#peaks'); pb.innerHTML = '';
    top.forEach((p, i) => { const el = document.createElement('div'); el.className = 'item'; el.style.cursor = 'pointer'; const itm = toITM(...p.ll); el.innerHTML = `<div class="t"><span>Peak ${i + 1} · ${fmt(p.z, 0)} m</span></div><div class="d">${p.ll[1].toFixed(5)}, ${p.ll[0].toFixed(5)}  ITM ${fmt(itm[0])}/${fmt(itm[1])}</div>`; el.onclick = () => map.easeTo({ center: p.ll, zoom: Math.max(map.getZoom(), 13) }); pb.append(el); });
    analysis = { Z, sl, asp, nx, ny, cs };
  } finally { busy(false); }
}
function clearAnalysis() { analysis = null; setImage('mask', null); map.getSource('contours').setData(E()); map.getSource('peaks').setData(E()); $('#anInfo').textContent = ''; $('#peaks').innerHTML = ''; }

/* ---------- projects ---------- */
const LS = 'terrainlab.projects.v1';
const loadAll = () => { try { return JSON.parse(localStorage.getItem(LS)) || {}; } catch (_) { return {}; } };
const saveAll = o => { try { localStorage.setItem(LS, JSON.stringify(o)); } catch (_) { toast('Could not save in this browser (storage blocked). Use Export file.'); } };
let saveT = 0;
function autosave() { clearTimeout(saveT); saveT = setTimeout(saveNow, 600); }
function slim() { const o = JSON.parse(JSON.stringify(S, (k, v) => (k === 'res' || k === 'itm' || k === 'elev' || k === 'flat' || k === 'surf' || k === 'up' || k === 'down' || k === 'label' && false) ? undefined : v)); if (o.viewshed) delete o.viewshed.res, delete o.viewshed.bbox; return o; }
function saveNow() { readSettings(); const a = loadAll(); if (!PID) PID = uid(); a[PID] = { name: S.name, updated: Date.now(), data: slim() }; saveAll(a); try { localStorage.setItem(LS + '.cur', PID); } catch (_) { } renderProjects(); }
function renderProjects() {
  const a = loadAll(), box = $('#projList'); box.innerHTML = '';
  Object.entries(a).sort((x, y) => y[1].updated - x[1].updated).forEach(([id, p]) => {
    const el = document.createElement('div'); el.className = 'item' + (id === PID ? '' : ' m'); el.innerHTML = `<div class="t"><span>${p.name.replace(/</g, '&lt;')}</span><span></span></div><div class="d">${new Date(p.updated).toLocaleString()} · ${p.data.items.length} items${id === PID ? ' · open' : ''}</div>`;
    const bt = el.querySelector('.t span:last-child');
    const add = (t, fn) => { const b = document.createElement('button'); b.className = 'mini'; b.textContent = t; b.onclick = fn; bt.append(b); };
    if (id !== PID) add('open', () => openProject(id));
    add('✕', () => { if (!confirm(`Delete project "${p.name}"?`)) return; const x = loadAll(); delete x[id]; saveAll(x); if (id === PID) { PID = null; } renderProjects(); });
    box.append(el);
  });
}
async function openProject(id, data) {
  const a = loadAll(); const d = data || (a[id] && a[id].data); if (!d) return;
  for (const k of [...markers.keys()]) clearMarkers(k); clearAnalysis(); setImage('viewshed', null); $('#profile').hidden = true;
  S = Object.assign(DEF(), d); S.set = Object.assign(DEF().set, d.set || {}); PID = id || null; if (!id) PID = null;
  writeSettings(); applySettings(); $('#projName').value = S.name;
  map.jumpTo({ center: S.view.center, zoom: S.view.zoom, pitch: S.view.pitch, bearing: S.view.bearing });
  busy(true); await tick();
  try { for (const it of S.items) await compute(it); if (S.roi) await roiChanged(); renderAll(); if (S.viewshed) await runViewshed(); } finally { busy(false); }
  if (!PID) saveNow(); else { try { localStorage.setItem(LS + '.cur', PID); } catch (_) { } renderProjects(); }
}
const SETIDS = { base: 'base', l3d: 'l3d', ex: 'ex', pitch: 'pitch', lHs: 'lHs', hs: 'hs', lItems: 'lItems', lRoi: 'lRoi', lView: 'lView', lCont: 'lCont', lMask: 'lMask', lPeaks: 'lPeaks', hObs: 'hObs', hTgt: 'hTgt', vsR: 'vsR', curv: 'curv', fEmin: 'fEmin', fEmax: 'fEmax', fSmin: 'fSmin', fSmax: 'fSmax', fAsp: 'fAsp', fAspTol: 'fAspTol', cInt: 'cInt', pSep: 'pSep' };
function readSettings() { for (const [k, id] of Object.entries(SETIDS)) { const el = $('#' + id); S.set[k] = el.type === 'checkbox' ? el.checked : el.value; } S.name = $('#projName').value || 'Untitled'; }
function writeSettings() { for (const [k, id] of Object.entries(SETIDS)) { const el = $('#' + id); if (el.type === 'checkbox') el.checked = !!S.set[k]; else el.value = S.set[k]; } $('#exOut').textContent = (+S.set.ex).toFixed(1) + '×'; $('#pitchOut').textContent = S.set.pitch + '°'; $('#hsOut').textContent = (+S.set.hs).toFixed(2); }

/* ---------- wire up ---------- */
function bind() {
  $$('#tabs button').forEach(b => b.onclick = () => { $$('#tabs button').forEach(x => x.classList.toggle('on', x === b)); $$('.pane').forEach(p => p.classList.toggle('on', p.dataset.pane === b.dataset.tab)); });
  $$('[data-tool]').forEach(b => b.onclick = () => { setTool(b.dataset.tool); b.blur(); document.body.classList.remove('open'); });
  $('#fFinish').onclick = finish; $('#fCancel').onclick = () => setTool(null);
  const hh = () => document.documentElement.style.setProperty('--hh', $('#top').offsetHeight + 'px'); hh(); addEventListener('resize', hh); if (document.fonts) document.fonts.ready.then(hh);
  $('#btnPanel').onclick = () => document.body.classList.toggle('open');
  const lay = () => { readSettings(); $('#exOut').textContent = (+S.set.ex).toFixed(1) + '×'; $('#pitchOut').textContent = S.set.pitch + '°'; $('#hsOut').textContent = (+S.set.hs).toFixed(2); applySettings(); autosave(); };
  ['base', 'l3d', 'ex', 'hs', 'lHs', 'lItems', 'lRoi', 'lView', 'lCont', 'lMask', 'lPeaks'].forEach(id => $('#' + id).addEventListener('input', lay));
  $('#pitch').addEventListener('input', () => { lay(); map.setPitch(+$('#pitch').value); });
  ['hObs', 'hTgt', 'curv'].forEach(id => $('#' + id).addEventListener('change', async () => { readSettings(); autosave(); for (const it of S.items.filter(i => i.type === 'line')) await update(it); if (S.viewshed) await runViewshed(); }));
  $('#vsR').addEventListener('change', async () => { readSettings(); autosave(); if (S.viewshed) await runViewshed(); });
  ['fEmin', 'fEmax', 'fSmin', 'fSmax', 'fAsp', 'fAspTol', 'cInt', 'pSep'].forEach(id => $('#' + id).addEventListener('change', () => { readSettings(); autosave(); }));
  $('#roiR').addEventListener('input', () => { $('#roiROut').textContent = $('#roiR').value; if (S.roi && S.roi.type === 'circle') { S.roi.r = +$('#roiR').value * 1000; roiChanged(); } });
  $('#roiSet').onclick = async () => { const c = parseCoord($('#roiCenter').value); if (!c) return toast('Enter "lat, lon" or "ITM E N".'); S.roi = { type: 'circle', c: [c.lon, c.lat], r: +$('#roiR').value * 1000 }; await roiChanged(); zoomRoi(); };
  $('#roiZoom').onclick = zoomRoi; $('#roiClear').onclick = () => { S.roi = null; $('#roiInfo').textContent = ''; renderAll(false); autosave(); };
  $('#btnAnalyse').onclick = runAnalysis; $('#btnAnClear').onclick = clearAnalysis;
  $('#btnClearItems').onclick = () => { if (!S.items.length && !S.viewshed) return; S.items.forEach(i => clearMarkers(i.id)); clearMarkers('vs'); S.items = []; S.viewshed = null; setImage('viewshed', null); $('#profile').hidden = true; renderAll(false); autosave(); };
  $('#goForm').onsubmit = async ev => { ev.preventDefault(); const c = parseCoord($('#goInput').value); if (!c) return toast('Enter "lat, lon" or "ITM E N".'); map.flyTo({ center: [c.lon, c.lat], zoom: 13 }); const it = { id: uid(), type: 'point', ll: [c.lon, c.lat], label: 'P' + (S.items.filter(i => i.type === 'point').length + 1) }; S.items.push(it); await update(it); renderItems(); };
  $('#btnRefresh').onclick = async () => { busy(true); await tick(); DEM.clear(); for (const sid of ['dem', 'hs']) { const src = map.getSource(sid); if (src && src.reload) src.reload(); } try { for (const it of S.items) await compute(it); renderAll(); if (S.roi) await roiChanged(); if (S.viewshed) await runViewshed(); if (analysis) await runAnalysis(); toast('Refreshed'); } finally { busy(false); } };
  $('#projName').addEventListener('input', () => { S.name = $('#projName').value; autosave(); });
  $('#pNew').onclick = () => { if (PID) saveNow(); openProject(null, DEF()); };
  $('#pSave').onclick = () => { saveNow(); toast('Saved'); };
  $('#pExport').onclick = () => { readSettings(); const b = new Blob([JSON.stringify({ app: 'israel-terrain-lab', ...slim() }, null, 1)], { type: 'application/json' }), a = document.createElement('a'); a.href = URL.createObjectURL(b); a.download = (S.name || 'project').replace(/[^\w\-]+/g, '_') + '.terrain.json'; a.click(); };
  $('#pImport').onclick = () => $('#pFile').click();
  $('#pFile').onchange = async ev => { const f = ev.target.files[0]; if (!f) return; try { const d = JSON.parse(await f.text()); if (!d.items || !d.set) throw 0; await openProject(null, d); } catch (_) { toast('Not a Terrain Lab project file.'); } ev.target.value = ''; };
  window.addEventListener('keydown', ev => { if (ev.target.tagName === 'INPUT' || ev.target.tagName === 'SELECT') return; if (ev.key === 'Escape') setTool(null); if (ev.key === 'Enter') { ev.preventDefault(); finish(); } });
}
function zoomRoi() { if (!S.roi) return; const b = roiBBox(); map.fitBounds([[b[0], b[1]], [b[2], b[3]]], { padding: 60 }); }

/* ---------- start ---------- */
(async function () {
  let cur = null; try { cur = localStorage.getItem(LS + '.cur'); } catch (_) { }
  const a = loadAll(); if (cur && a[cur]) { S = Object.assign(DEF(), a[cur].data); S.set = Object.assign(DEF().set, a[cur].data.set || {}); PID = cur; }
  writeSettings(); $('#projName').value = S.name; bind(); initMap();
  $('#dataNote').textContent = 'Elevation: Copernicus GLO-30 DSM (ESA, TanDEM-X 2011–2015, 30 m, vertical accuracy about ±2–4 m). It is a surface model: tree tops and buildings are included.';
  renderProjects();
  map.on('load', async () => { busy(true); await tick(); try { for (const it of S.items) await compute(it); if (S.roi) await roiChanged(); renderAll(); if (S.viewshed) await runViewshed(); } finally { busy(false); } });
  window.__app = { S: () => S, DEM, lineProfile, parseCoord, hav, runViewshed };
})();

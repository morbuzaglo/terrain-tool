"""Build the Israel elevation tile pyramid from Copernicus GLO-30 (ESA/TanDEM-X, 30 m DSM).

Downloads the 1x1 degree COG tiles from the public AWS bucket (cached in %LOCALAPPDATA%/drill-spl-cache),
resamples them to Web-Mercator XYZ tiles in Terrarium PNG encoding (R*256 + G + B/256 - 32768 = metres)
and writes tiles/{z}/{x}/{y}.png  (z = 7..12).  Usage: python tools/build_tiles.py
"""
import os, math, sys, io, json
import numpy as np, requests, tifffile
from PIL import Image
from scipy.ndimage import map_coordinates

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(os.environ.get("LOCALAPPDATA", HERE), "drill-spl-cache")
OUT = os.path.join(HERE, "tiles")
LAT0, LAT1, LON0, LON1 = 29.35, 33.50, 34.15, 36.00      # Israel + West Bank + Gaza + Golan + Dead Sea margin
ZMAX, ZMIN = 12, 7
BASE = "https://copernicus-dem-30m.s3.amazonaws.com/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM/Copernicus_DSM_COG_10_{ns}{lat:02d}_00_{ew}{lon:03d}_00_DEM.tif"
os.makedirs(CACHE, exist_ok=True)

def get_tile(lat, lon):
    p = os.path.join(CACHE, f"N{lat:02d}E{lon:03d}.tif")
    if not os.path.exists(p):
        r = requests.get(BASE.format(ns="N", lat=lat, ew="E", lon=lon), stream=True, timeout=120)
        if r.status_code != 200:
            print("missing", lat, lon, r.status_code); return None
        with open(p + ".part", "wb") as f:
            for ch in r.iter_content(1 << 20): f.write(ch)
        os.replace(p + ".part", p)
    return tifffile.imread(p).astype(np.float32)

# ---- mosaic ---------------------------------------------------------------------------------------------
lats = list(range(int(math.floor(LAT0)), int(math.floor(LAT1)) + 1))
lons = list(range(int(math.floor(LON0)), int(math.floor(LON1)) + 1))
N = 3600
mos = np.zeros(((lats[-1] - lats[0] + 1) * N, (lons[-1] - lons[0] + 1) * N), np.float32)   # row 0 = northern edge
for la in lats:
    for lo in lons:
        a = get_tile(la, lo)
        r0 = (lats[-1] - la) * N; c0 = (lo - lons[0]) * N
        if a is not None:
            a = a[:N, :N] if a.shape[0] >= N else np.pad(a, ((0, N - a.shape[0]), (0, N - a.shape[1])))
            a = np.where(a < -500, 0, a)
            mos[r0:r0 + N, c0:c0 + N] = a[:N, :N]
        print("tile", la, lo, None if a is None else (a.shape, float(a.min()), float(a.max())), flush=True)
print("mosaic", mos.shape)
TOP = lats[-1] + 1.0     # latitude of mosaic row 0 edge
LEFT = float(lons[0])    # longitude of mosaic column 0 edge
DEG = 1.0 / N            # degrees per pixel (both axes; GLO-30 is 1 arcsec at this latitude)

# ---- XYZ tiles ------------------------------------------------------------------------------------------
def lonlat2tile(lon, lat, z):
    n = 2 ** z; x = (lon + 180) / 360 * n
    y = (1 - math.asinh(math.tan(math.radians(lat))) / math.pi) / 2 * n
    return x, y
def tile_arr(z, x, y):
    n = 2 ** z
    px = (x + (np.arange(256) + 0.5) / 256) / n
    py = (y + (np.arange(256) + 0.5) / 256) / n
    lon = px * 360 - 180
    lat = np.degrees(np.arctan(np.sinh(math.pi * (1 - 2 * py))))
    cc = (lon - LEFT) / DEG - 0.5
    rr = (TOP - lat) / DEG - 0.5
    C, R = np.meshgrid(cc, rr)
    return map_coordinates(mos, [R, C], order=1, mode="nearest").astype(np.float32)

def encode(a):
    v = np.clip(a + 32768.0, 0, 65535.99)
    r = np.floor(v / 256); g = np.floor(v % 256); b = np.floor((v - np.floor(v)) * 16) * 16
    return np.dstack([r, g, b]).astype(np.uint8)

x0, y1 = lonlat2tile(LON0, LAT0, ZMAX); x1, y0 = lonlat2tile(LON1, LAT1, ZMAX)
X0, X1, Y0, Y1 = int(x0), int(x1), int(y0), int(y1)
print("z12 tile range", X0, X1, Y0, Y1, (X1 - X0 + 1) * (Y1 - Y0 + 1), "tiles")
cur = {}
for x in range(X0, X1 + 1):
    for y in range(Y0, Y1 + 1):
        cur[(x, y)] = tile_arr(ZMAX, x, y)
total = 0
def save(z, d):
    global total
    for (x, y), a in d.items():
        p = os.path.join(OUT, str(z), str(x)); os.makedirs(p, exist_ok=True)
        buf = io.BytesIO(); Image.fromarray(encode(a)).save(buf, "PNG", optimize=True)
        with open(os.path.join(p, f"{y}.png"), "wb") as f: f.write(buf.getvalue())
        total += len(buf.getvalue())
save(ZMAX, cur)
for z in range(ZMAX - 1, ZMIN - 1, -1):
    par = {}
    for (x, y) in cur:
        par.setdefault((x // 2, y // 2), None)
    new = {}
    for (px, py) in par:
        big = np.zeros((512, 512), np.float32)
        for dx in (0, 1):
            for dy in (0, 1):
                a = cur.get((px * 2 + dx, py * 2 + dy))
                if a is None: a = np.zeros((256, 256), np.float32)
                big[dy * 256:(dy + 1) * 256, dx * 256:(dx + 1) * 256] = a
        new[(px, py)] = big.reshape(256, 2, 256, 2).mean(axis=(1, 3))
    cur = new
    save(z, cur)
    print("z", z, len(cur), "tiles", flush=True)
meta = dict(source="Copernicus GLO-30 DSM (ESA, TanDEM-X 2011-2015), release 2022_1", encoding="terrarium", zmin=ZMIN, zmax=ZMAX,
            bounds=[LON0, LAT0, LON1, LAT1], z12_range=[X0, X1, Y0, Y1], bytes=total)
json.dump(meta, open(os.path.join(OUT, "meta.json"), "w"), indent=1)
print("done", total / 1e6, "MB")

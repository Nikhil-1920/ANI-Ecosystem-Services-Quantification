"""
ANI Ecosystem Services — 3D Web-GL Topographic Hotspot Visualization
=====================================================================
Drapes the Ecosystem Collapse Index (ECI) over the 30 m SRTM terrain
using Plotly Surface to produce a fully interactive 3-D map where
triple-collapse hotspots glow in magma tones against lush green hills.

Inputs  : data/processed/ANI_SRTM_DEM_30m_clipped.tif
          results/eci_collapse_hotspots.tif
Outputs : presentation/ani_3d_hotspots_flythrough.html

Run with: venv/bin/python src/render_3d_hotspots.py
"""

import warnings
warnings.filterwarnings('ignore')

import numpy as np
import rasterio
import plotly.graph_objects as go
from scipy.ndimage import uniform_filter
from pathlib import Path

# ── Directory Paths ────────────────────────────────────────────────────
SCRIPT_DIR   = Path(__file__).parent
PROC_DIR     = SCRIPT_DIR.parent / 'data' / 'processed'
RES_DIR      = SCRIPT_DIR.parent / 'results'
PRESENT_DIR  = SCRIPT_DIR.parent / 'presentation'

# ── Print Separator ────────────────────────────────────────────────────
SEP = '=' * 60

# ── Rendering Parameters ──────────────────────────────────────────────
EPICENTER_SEARCH_RADIUS = 500    # Uniform filter size for density smoothing (px)
CROP_HALF_BOX           = 2000   # Half-width of the focal bounding box (px → 60 km)
WEBGL_DOWNSAMPLE        = 4      # Decimation factor for Web-GL surface fluidity
SURFACE_TERRAIN_MAX     = 0.40   # Terrain occupies colour range 0.0–0.40
HOTSPOT_COLOUR_MIN      = 0.60   # Hotspot glow begins at 0.60
HOTSPOT_ECI_THRESHOLD   = 0.50   # ECI score above which hotspot colouring activates
Z_EXAGGERATION          = 0.40   # Vertical exaggeration ratio (terrain drama)

# ── Custom Colourscale: terrain transitions → hotspot glow ────────────
SURFACE_COLOURSCALE = [
    # Normal terrain (0.0 → 0.40)
    [0.000, 'rgb(10, 40, 60)'],     # Ocean deep blue
    [0.020, 'rgb(30, 80, 50)'],     # Coastal / mangrove green
    [0.200, 'rgb(46, 125, 50)'],    # Mid-elevation lush green
    [0.400, 'rgb(141, 110, 99)'],   # High-elevation mountain brown
    # Bridge zone (no data in this gap — avoid visible discontinuity)
    [0.401, 'rgb(141, 110, 99)'],
    [0.599, 'rgb(141, 110, 99)'],
    # Ecosystem collapse hotspots (0.60 → 1.00)
    [0.600, 'rgb(255, 150, 0)'],    # Warning glow — orange
    [0.800, 'rgb(255, 50, 0)'],     # High intensity — deep red
    [1.000, 'rgb(255, 0, 100)'],    # Catastrophic — magenta
]


# ══════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ══════════════════════════════════════════════════════════════════════
if __name__ == '__main__':
    print(SEP)
    print("  ANI Ecosystem Services — 3D Web-GL Hotspot Visualization")
    print(SEP)

    # 1. Load elevation and ECI data ───────────────────────────────────
    print("  Loading 30 m DEM & ECI hotspot rasters …")

    dem_path = PROC_DIR / 'ANI_SRTM_DEM_30m_clipped.tif'
    eci_path = RES_DIR  / 'eci_collapse_hotspots.tif'

    if not dem_path.exists() or not eci_path.exists():
        print("  ❌  Required data missing — run the full pipeline first.")
        exit(1)

    with rasterio.open(dem_path) as dem_src:
        dem_data    = dem_src.read(1).astype(float)
        dem_nodata  = dem_src.nodata

    with rasterio.open(eci_path) as eci_src:
        eci_data    = eci_src.read(1).astype(float)
        eci_nodata  = eci_src.nodata

    # Clean invalid pixels
    dem_data = np.where(dem_data == dem_nodata, np.nan, dem_data)
    dem_data = np.where(dem_data < -10, np.nan, dem_data)   # Remove ocean artefacts
    eci_data = np.where(eci_data == eci_nodata, 0.0, eci_data)
    eci_data = np.nan_to_num(eci_data, nan=0.0)

    # 2. Locate highest-density epicentre ──────────────────────────────
    print("  Searching for maximum collapse epicentre …")
    density_smooth   = uniform_filter(eci_data, size=EPICENTER_SEARCH_RADIUS)
    epicentre_r, epicentre_c = np.unravel_index(
        np.argmax(density_smooth, axis=None), density_smooth.shape
    )

    # Crop a geographic window around the epicentre
    r_start = max(0,                 epicentre_r - CROP_HALF_BOX)
    r_end   = min(dem_data.shape[0], epicentre_r + CROP_HALF_BOX)
    c_start = max(0,                 epicentre_c - CROP_HALF_BOX)
    c_end   = min(dem_data.shape[1], epicentre_c + CROP_HALF_BOX)

    dem_cropped = dem_data[r_start:r_end, c_start:c_end]
    eci_cropped = eci_data[r_start:r_end, c_start:c_end]
    print(f"  Epicentre found. Cropped grid: {dem_cropped.shape} (~60 × 60 km focal zone).")

    # 3. Decimate for Web-GL performance ───────────────────────────────
    ds           = WEBGL_DOWNSAMPLE
    dem_surface  = dem_cropped[::ds, ::ds]
    eci_surface  = eci_cropped[::ds, ::ds]

    # Replace NaN ocean with flat near-zero (renders cleanly in plotly)
    dem_surface  = np.nan_to_num(dem_surface, nan=-5.0)
    print(f"  Decimated Web-GL grid: {dem_surface.shape} pixels.")

    # 4. Build hybrid surface colour array ─────────────────────────────
    print("  Draping ECI signature onto 3D elevation mesh …")
    dem_max_val    = np.nanmax(dem_surface) if np.nanmax(dem_surface) > 0 else 1.0
    surface_colour = (dem_surface / dem_max_val) * SURFACE_TERRAIN_MAX

    # Overwrite hotspot pixels with ECI-scaled glow (0.60 → 1.00)
    hotspot_pixels = eci_surface > HOTSPOT_ECI_THRESHOLD
    surface_colour[hotspot_pixels] = (
        HOTSPOT_COLOUR_MIN + (eci_surface[hotspot_pixels] * 0.40)
    )

    # 5. Build Plotly 3-D surface figure ───────────────────────────────
    surface_trace = go.Surface(
        z              = dem_surface,
        surfacecolor   = surface_colour,
        colorscale     = SURFACE_COLOURSCALE,
        cmin           = 0.0,
        cmax           = 1.0,
        showscale      = False,
        lighting       = dict(
            ambient   = 0.4,
            diffuse   = 0.6,
            roughness = 0.5,
            specular  = 0.5,
            fresnel   = 0.2,
        ),
    )

    fig = go.Figure(data=[surface_trace])
    fig.update_layout(
        title=dict(
            text  = '<b>Interactive 3D Ecosystem Collapse Hotspots</b>'
                    '<br>Andaman & Nicobar Terrain Rendering',
            font  = dict(size=20, color='white'),
            x     = 0.5,
            y     = 0.95,
        ),
        paper_bgcolor = '#0d0d1a',
        scene = dict(
            xaxis = dict(showgrid=False, showticklabels=False, title=''),
            yaxis = dict(showgrid=False, showticklabels=False, title=''),
            zaxis = dict(showgrid=False, showticklabels=False, title='',
                         range=[-50, np.nanmax(dem_surface)]),
            aspectratio = dict(x=1, y=1, z=Z_EXAGGERATION),
            bgcolor     = '#0d0d1a',
        ),
        margin = dict(l=0, r=0, b=0, t=60),
    )

    # 6. Export standalone HTML ────────────────────────────────────────
    PRESENT_DIR.mkdir(parents=True, exist_ok=True)
    html_out = PRESENT_DIR / 'ani_3d_hotspots_flythrough.html'
    fig.write_html(str(html_out), include_plotlyjs='cdn')

    print(f"\n{SEP}")
    print(f"  🌍  3D Web-GL Visualization saved → {html_out.name}")
    print(f"      Open in any web browser to explore interactively.")
    print(SEP + "\n")

import numpy as np
import os
import rasterio
from rasterio.enums import Resampling as resampleenum
from rasterio.warp import reproject as rasterreproj

# folder configuration
rawdirpath  = "data/raw"       
procdirpath = "data/processed" 

# match everything to loss layer grid
masterkeyid = "loss"                    
rasterfiles = {
    "loss":      "ani-lossyear.tif",
    "bio":       "ani-biomass.tif",
    "slope":     "ani-slope.tif",
    "tree":      "ani-treecover.tif",
    "landcover": "ani-landcover.tif"
}

print("--- robust data prep (aligned to loss grid via reproject) ---")

# make sure folders exist
os.makedirs(rawdirpath,  exist_ok=True)
os.makedirs(procdirpath, exist_ok=True)

# 1) open master and store exact target grid definition
masterpathv = os.path.join(rawdirpath, rasterfiles[masterkeyid])
if not os.path.exists(masterpathv):
    raise FileNotFoundError(f"missing master raster: {masterpathv}")

with rasterio.open(masterpathv) as masterrast:
    targetrows = masterrast.height
    targetcols = masterrast.width
    targettran = masterrast.transform
    targetcrsc = masterrast.crs

    print(f"target grid (from {masterkeyid}): {targetcols} x {targetrows}")
    print(f"target crs: {targetcrsc}")

    # nicely formatted affine transform
    print("target transform (affine):")
    print("  |{:9.2f} {:9.2f} {:11.2f} |".format(
        targettran.a, targettran.b, targettran.c
    ))
    print("  |{:9.2f} {:9.2f} {:11.2f} |".format(
        targettran.d, targettran.e, targettran.f
    ))
    print("  |{:9.2f} {:9.2f} {:11.2f} |".format(
        0.0, 0.0, 1.0
    ))

    # write config into processed folder
    configpath = os.path.join(procdirpath, "config.txt")
    with open(configpath, "w") as configfile:
        configfile.write(f"{targetcols} {targetrows}")

# 2) process each file into the target grid
for rasterkeyid, rastername in rasterfiles.items():
    rasterpathv = os.path.join(rawdirpath, rastername)

    if not os.path.exists(rasterpathv):
        print(f"missing: {rasterpathv} (please check your data/raw folder)")
        continue

    print(f"processing {rasterkeyid}...")

    # choose resampling method
    # - categorical: nearest
    # - continuous: bilinear
    if rasterkeyid in ["loss", "landcover"]:
        resamplemode = resampleenum.nearest
    elif rasterkeyid in ["tree"]:
        resamplemode = resampleenum.nearest
    else:
        resamplemode = resampleenum.bilinear

    with rasterio.open(rasterpathv) as sourcerast:
        # destination array (target grid)
        targetarray = np.full(
            (targetrows, targetcols),
            np.nan,
            dtype=np.float32
        )

        # reproject and resample into the loss grid exactly
        rasterreproj(
            source=rasterio.band(sourcerast, 1),
            destination=targetarray,
            src_transform=sourcerast.transform,
            src_crs=sourcerast.crs,
            dst_transform=targettran,
            dst_crs=targetcrsc,
            resampling=resamplemode,
            src_nodata=sourcerast.nodata,
            dst_nodata=np.nan
        )

        # convert nans for downstream c++ code
        targetarray = np.nan_to_num(targetarray, nan=0.0).astype(np.float32)
        outbinpath = os.path.join(procdirpath, f"{rasterkeyid}.bin")
        targetarray.tofile(outbinpath)
        print(f"  -> saved {outbinpath} (resampling={resamplemode})")

print("--- data prep complete. ready for ecovuln ---")
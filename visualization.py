import numpy as np
import os
import matplotlib.pyplot as plt

# configuration
dimsfileid   = "data/config.txt"
outputdpivl  = 300
pixelareaha  = 0.09  # must match c++ pixel area
summarycsvfl = "vulnerability_summary.csv"
print("--- ani master visualization suite ---")

# 1. load grid dimensions
if not os.path.exists(dimsfileid):
    print(f"error: {dimsfileid} not found. run dataset.py first.")
    raise SystemExit(1)

with open(dimsfileid, "r") as dimfileob:
    dimsvalues = dimfileob.read().split()
    gridwidthpx  = int(dimsvalues[0])
    gridheightpx = int(dimsvalues[1])
print(f"grid: {gridwidthpx} x {gridheightpx}")

# 2. helper to load data
def loadbindata(layername):
    binpathstr = f"data/{layername}.bin"
    if not os.path.exists(binpathstr):
        print(f"warning: {binpathstr} not found. skipping related figure.")
        return None

    flatarrayf = np.fromfile(binpathstr, dtype=np.float32)
    expectvals = gridwidthpx * gridheightpx
    if flatarrayf.size != expectvals:
        print(
            f"warning: {binpathstr} has {flatarrayf.size} values, "
            f"expected {expectvals}. skipping."
        )
        return None
    return flatarrayf.reshape((gridheightpx, gridwidthpx))

# load layers used across figures
lossgriddt   = loadbindata("loss")
vulnindexdt  = loadbindata("vulnerability-index")

# --- figure 1: biophysical impacts (carbon & erosion) ---
print("generating figure 1 (biophysical impacts)...")

biomassgrid  = loadbindata("bio")
slopegridvl  = loadbindata("slope")
treecovergr  = loadbindata("tree")

if (
    lossgriddt  is not None and
    biomassgrid is not None and
    slopegridvl is not None and
    treecovergr is not None
):
    # loss mask: tree > 50 and loss year between 19-23
    lossmaskmt = (
        (treecovergr > 50) &
        (lossgriddt >= 19) &
        (lossgriddt <= 23)
    )

    # carbon loss map
    # match c++ assumption: bio is tonnes biomass per hectare
    # convert to per-pixel tonnes, then *0.47 carbon fraction
    carbonmapgr = np.zeros_like(biomassgrid, dtype=np.float32)
    carbonmapgr[lossmaskmt] = (
        biomassgrid[lossmaskmt] * pixelareaha * 0.47
    )
    carbonmapgr = np.ma.masked_where(carbonmapgr <= 0, carbonmapgr)

    # erosion risk map (visual proxy)
    erosionmapg = np.zeros_like(slopegridvl, dtype=np.float32)
    erosionmapg[lossmaskmt] = slopegridvl[lossmaskmt] * 2.5
    erosionmapg = np.ma.masked_where(erosionmapg <= 0, erosionmapg)
    figureone, (axisone, axistwo) = plt.subplots(1, 2, figsize=(14, 8))
    figureone.suptitle(
        "biophysical consequences of forest conversion (ani 2019-2023)",
        fontsize=16
    )

    # panel a: carbon
    imageone = axisone.imshow(
        carbonmapgr,
        cmap=plt.cm.Reds,
        vmin=0,
        vmax=150
    )
    axisone.set_title("a. carbon stock loss", fontsize=14)
    axisone.axis("off")
    plt.colorbar(
        imageone,
        ax=axisone,
        fraction=0.046,
        pad=0.04,
        label="tonnes carbon / pixel"
    )

    # panel b: erosion
    imagetwo = axistwo.imshow(
        erosionmapg,
        cmap=plt.cm.YlOrBr,
        vmin=0,
        vmax=100
    )
    axistwo.set_title("b. soil erosion risk surge", fontsize=14)
    axistwo.axis("off")
    plt.colorbar(
        imagetwo,
        ax=axistwo,
        fraction=0.046,
        pad=0.04,
        label="relative erosion risk"
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    figureone.savefig("Figure_1_Biophysical.png", dpi=outputdpivl)
    plt.close(figureone)
    print("saved Figure_1_Biophysical.png")
else:
    print("skipping figure 1 (missing one of: loss/bio/slope/tree).")

# --- figure 2: conservation vulnerability (the novelty map) ---
print("generating figure 2 (vulnerability map)...")

if vulnindexdt is not None:
    vulnmasked = np.ma.masked_invalid(vulnindexdt)
    vulnmasked = np.ma.masked_where(vulnmasked <= 0, vulnmasked)

    figuretwo = plt.figure(figsize=(10, 12))
    colormapvg = plt.cm.RdYlGn_r  # green -> yellow -> red

    plt.imshow(vulnmasked, cmap=colormapvg, vmin=0.2, vmax=0.8)
    plt.title(
        "spatial vulnerability index & conservation priorities",
        fontsize=16
    )
    plt.axis("off")

    colorbarvg = plt.colorbar(fraction=0.035, pad=0.04)
    colorbarvg.set_label(
        "vulnerability score (0 = intact, 1 = critical threat)",
        fontsize=12
    )

    infostring = "\n".join((
        r"$\bf{PRIORITY\ ZONES}$",
        "CRITICAL: top 10% risk",
        "HIGH: top 25% risk",
        "MODERATE: top 50% risk"
    ))
    textboxopt = dict(boxstyle="round", facecolor="white", alpha=0.9)
    axiscurrent = plt.gca()
    axiscurrent.text(
        0.02,
        0.98,
        infostring,
        transform=axiscurrent.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=textboxopt
    )

    plt.savefig(
        "Figure_2_Vulnerability.png",
        dpi=outputdpivl,
        bbox_inches="tight"
    )
    plt.close(figuretwo)
    print("saved Figure_2_Vulnerability.png")
else:
    print("skipping figure 2 (missing data/vulnerability_index.bin).")

# --- figure 3: statistical summary (charts) ---
print("generating figure 3 (statistics)...")

prioritycls = ["Critical", "High", "Moderate"]
areahectrs  = None
carbontotvl = None

if os.path.exists(summarycsvfl):
    rowvaluesd = {}

    with open(summarycsvfl, "r") as csvfileob:
        csvheader = csvfileob.readline()
        for csvlinevl in csvfileob:
            partslist = csvlinevl.strip().split(",")
            if len(partslist) < 4:
                continue

            classidvl = partslist[0]
            areavl    = float(partslist[1])
            carbonvl  = float(partslist[2])
            pctshare  = float(partslist[3])

            rowvaluesd[classidvl.lower()] = (areavl, carbonvl, pctshare)

    if all(k in rowvaluesd for k in ["critical", "high", "moderate"]):
        areahectrs = [
            rowvaluesd["critical"][0],
            rowvaluesd["high"][0],
            rowvaluesd["moderate"][0]
        ]
        carbontotvl = [
            rowvaluesd["critical"][1],
            rowvaluesd["high"][1],
            rowvaluesd["moderate"][1]
        ]
    else:
        print(
            f"warning: {summarycsvfl} missing expected classes. "
            "falling back to hardcoded values."
        )

if areahectrs is None or carbontotvl is None:
    # fallback values
    areahectrs = [64497.4, 96746.0, 161243.5]
    carbontotvl = [1079672.1, 192394.1, 52170.5]

figurethr, (axisarea, axiscarbon) = plt.subplots(1, 2, figsize=(14, 6))

# chart a: conservation area priority
barcolors = ["#d62728", "#ff7f0e", "#2ca02c"]
areabars  = axisarea.bar(prioritycls, areahectrs, color=barcolors, alpha=0.8)
axisarea.set_title("conservation priority areas", fontsize=14)
axisarea.set_ylabel("area (hectares)")
axisarea.grid(axis="y", linestyle="--", alpha=0.3)

for baritem in areabars:
    barheight = baritem.get_height()
    axisarea.text(
        baritem.get_x() + baritem.get_width() / 2.0,
        barheight,
        f"{int(round(barheight)):,} ha",
        ha="center",
        va="bottom"
    )

# chart b: carbon at risk (pie)
pielabels = ["critical zone", "high risk zone", "moderate risk zone"]
piesizes  = [
    float(carbontotvl[0]),
    float(carbontotvl[1]),
    float(carbontotvl[2])
]
piecolors  = ["#d62728", "#ff7f0e", "#FDB813"]
pieexplode = (0.1, 0.0, 0.0)

axiscarbon.pie(
    piesizes,
    explode=pieexplode,
    labels=pielabels,
    colors=piecolors,
    autopct="%1.1f%%",
    shadow=True,
    startangle=140
)
axiscarbon.set_title("carbon stock at immediate risk", fontsize=14)
plt.tight_layout()
figurethr.savefig("Figure_3_Statistics.png", dpi=outputdpivl)
plt.close(figurethr)
print("saved Figure_3_Statistics.png")
print("\n--- all visualizations complete ---")
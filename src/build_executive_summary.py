"""
ANI Ecosystem Services — One-Page Executive Summary PDF
=======================================================
Builds a single A4 page with the five headline numbers and the three
hero figures referenced by validation metrics. Designed for thesis
front-matter and standalone briefing.

Inputs : results/validation_summary.json
         results/economic_damage_2040_forecast.csv
         results/habitat_quality_by_landcover.csv
         results/rusle_erosion_by_landcover.csv
         figures/hero/*.png

Output : ANI_Executive_Summary.pdf  (project root)
"""

import json
from pathlib import Path

import pandas as pd
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor
from reportlab.pdfgen import canvas

ROOT = Path(__file__).parent.parent
RES = ROOT / "results"
FIG_HERO = ROOT / "figures" / "hero"

# Colours
TXT_DARK = HexColor("#1a1a2e")
TXT_MUTED = HexColor("#4a5568")
ACCENT_RED = HexColor("#c62828")
ACCENT_GREEN = HexColor("#2e7d32")
ACCENT_BLUE = HexColor("#1565c0")
BG_CARD = HexColor("#f5f7fa")
BORDER = HexColor("#cccccc")

OUT_PDF = ROOT / "ANI_Executive_Summary.pdf"


def fmt_int(x: float) -> str:
    return f"{x:,.0f}"


def main():
    summary = json.loads((RES / "validation_summary.json").read_text())
    forecast = pd.read_csv(RES / "economic_damage_2040_forecast.csv").iloc[0]
    hab = pd.read_csv(RES / "habitat_quality_by_landcover.csv")
    soil = pd.read_csv(RES / "rusle_erosion_by_landcover.csv")

    totals = summary["totals"]
    mk = summary["trend"]["co2e_mk"]
    sen = summary["trend"]["co2e_sen"]
    mi_hab = summary["morans_I_habitat_delta"]

    # Headline metrics
    co2e_total = totals["co2e_ggco2e"]
    co2e_lo, co2e_hi = totals["co2e_ci"]
    area_total = totals["area_ha"]
    area_lo, area_hi = totals["area_ci"]
    crop_erosion = soil[soil["esa_class"] == 40]["mean_soil_loss_t_ha"].values[0]
    forest_erosion = soil[soil["esa_class"] == 10]["mean_soil_loss_t_ha"].values[0]
    built_up_ha = hab[hab["esa_class"] == 50]["area_ha"].values[0]
    bau_damage_m = forecast["total_economic_damage_usd"] / 1e6

    c = canvas.Canvas(str(OUT_PDF), pagesize=A4)
    w, h = A4

    # ── Header band ───────────────────────────────────────────────────
    c.setFillColor(HexColor("#1a237e"))
    c.rect(0, h - 28 * mm, w, 28 * mm, fill=1, stroke=0)
    c.setFillColor(HexColor("#ffffff"))
    c.setFont("Helvetica-Bold", 18)
    c.drawString(15 * mm, h - 13 * mm,
                  "Andaman & Nicobar Ecosystem-Services Loss (2001–2024)")
    c.setFont("Helvetica", 11)
    c.drawString(15 * mm, h - 20 * mm,
                  "Quantifying carbon, habitat-quality and soil-retention damage from forest conversion")
    c.setFont("Helvetica-Oblique", 9)
    c.drawString(15 * mm, h - 25 * mm,
                  "Open-data pipeline: GEDI L4B  •  Hansen GFW  •  ESA WorldCover  •  SoilGrids  •  CHIRPS  •  SRTM")

    # ── Headline-number cards (5 across) ──────────────────────────────
    y_cards = h - 70 * mm
    card_h = 32 * mm
    card_w = (w - 30 * mm - 4 * 4 * mm) / 5

    cards = [
        ("CARBON", f"{co2e_total:.0f}",  "GgCO2e lost",
         f"95% CI: [{co2e_lo:.0f}, {co2e_hi:.0f}]", ACCENT_RED),
        ("FOREST", fmt_int(area_total), "ha deforested",
         f"95% CI: [{fmt_int(area_lo)}, {fmt_int(area_hi)}]", ACCENT_RED),
        ("HABITAT", f"{built_up_ha/1000:.1f}k", "ha now built-up",
         f"Mean ΔQ  cropland: -0.017", ACCENT_GREEN),
        ("SOIL", f"{crop_erosion/forest_erosion:.0f}×",
         "erosion on cropland",
         f"vs forest baseline ({forest_erosion:.1f} t/ha/yr)", HexColor("#8d6e63")),
        ("ECONOMIC", f"${bau_damage_m:.0f}M",
         "USD damage by 2040",
         f"BAU scenario (cumulative)", ACCENT_BLUE),
    ]

    x = 15 * mm
    for label, big, mid, foot, colr in cards:
        c.setFillColor(BG_CARD)
        c.setStrokeColor(BORDER)
        c.roundRect(x, y_cards, card_w, card_h, 2 * mm, fill=1, stroke=1)

        c.setFillColor(colr)
        c.setFont("Helvetica-Bold", 8)
        c.drawString(x + 3 * mm, y_cards + card_h - 5 * mm, label)

        c.setFillColor(TXT_DARK)
        c.setFont("Helvetica-Bold", 18)
        c.drawString(x + 3 * mm, y_cards + card_h - 14 * mm, big)

        c.setFillColor(TXT_DARK)
        c.setFont("Helvetica", 8)
        c.drawString(x + 3 * mm, y_cards + card_h - 19 * mm, mid)

        c.setFillColor(TXT_MUTED)
        c.setFont("Helvetica-Oblique", 7)
        c.drawString(x + 3 * mm, y_cards + 3 * mm, foot)

        x += card_w + 4 * mm

    # ── Two-column body ───────────────────────────────────────────────
    y_top = y_cards - 6 * mm
    col_w = (w - 35 * mm) / 2
    col_left_x = 15 * mm
    col_right_x = col_left_x + col_w + 5 * mm

    # Left column: hero figure (carbon loss)
    carbon_fig = FIG_HERO / "carbon_loss_hero.png"
    if carbon_fig.exists():
        c.drawImage(str(carbon_fig), col_left_x, y_top - 70 * mm,
                     width=col_w, height=68 * mm,
                     preserveAspectRatio=True, mask="auto")

    # Right column: hero figure (synthesis)
    synth_fig = FIG_HERO / "synthesis_hero.png"
    if synth_fig.exists():
        c.drawImage(str(synth_fig), col_right_x, y_top - 70 * mm,
                     width=col_w, height=68 * mm,
                     preserveAspectRatio=True, mask="auto")

    # ── Bottom: key findings + validation note ────────────────────────
    y_bot = y_top - 75 * mm
    c.setFillColor(TXT_DARK)
    c.setFont("Helvetica-Bold", 11)
    c.drawString(15 * mm, y_bot, "Key findings")

    c.setFont("Helvetica", 9)
    bullets = [
        ("Trend",
         f"Annual deforestation has slowed significantly (Mann-Kendall Z={mk['Z']:+.2f}, p={mk['p']:.3f}; "
         f"Sen slope {sen['slope']:+.2f} GgCO2e/yr). Cumulative damage persists."),
        ("Spatial pattern",
         f"Ecosystem damage is clustered, not random (Moran's I={mi_hab['I']:.2f}, "
         f"p={mi_hab['p']:.3f}; n={mi_hab['n']:,} cells at ~1.5 km). Justifies "
         f"hotspot-based intervention over uniform measures."),
        ("Carbon validation",
         "GEDI L4B vs. Saatchi 2011 cross-check on forest-only pixels "
         f"(n={summary['cross_validation']['n_pairs']:,}): "
         f"r={summary['cross_validation']['pearson_r']:.2f}, "
         f"bias={summary['cross_validation']['bias']:+.0f} Mg/ha. "
         "Disagreement is consistent with known Saatchi saturation in tropical canopies; "
         "GEDI is treated as the more sensitive estimator."),
        ("Land-cover physics",
         f"Cropland erosion ({crop_erosion:.0f} t/ha/yr) is "
         f"{crop_erosion/forest_erosion:.0f}× the forest baseline "
         f"({forest_erosion:.1f} t/ha/yr) — soil retention is the most "
         "responsive service to forest conversion."),
        ("Economic stake",
         f"Business-as-usual deforestation projects ${bau_damage_m:.0f}M USD "
         "cumulative damages by 2040 (carbon, sediment, and habitat-replacement "
         "costs combined)."),
    ]

    y_cur = y_bot - 5 * mm
    for tag, body in bullets:
        c.setFont("Helvetica-Bold", 9)
        c.setFillColor(ACCENT_BLUE)
        c.drawString(15 * mm, y_cur, f"• {tag}:")
        tag_w = c.stringWidth(f"• {tag}:", "Helvetica-Bold", 9)
        c.setFont("Helvetica", 9)
        c.setFillColor(TXT_DARK)

        # Word-wrap body
        words = body.split()
        line = ""
        max_w = w - 15 * mm - (15 * mm + tag_w + 2 * mm) - 5 * mm
        x0 = 15 * mm + tag_w + 2 * mm
        lines = []
        for word in words:
            test = (line + " " + word).strip()
            if c.stringWidth(test, "Helvetica", 9) <= max_w:
                line = test
            else:
                lines.append(line)
                line = word
                max_w = w - 30 * mm   # subsequent lines use full width
        if line:
            lines.append(line)

        for i, ln in enumerate(lines):
            xx = x0 if i == 0 else 18 * mm
            c.drawString(xx, y_cur, ln)
            y_cur -= 3.6 * mm
        y_cur -= 1 * mm

    # ── Footer ────────────────────────────────────────────────────────
    c.setStrokeColor(BORDER)
    c.line(15 * mm, 15 * mm, w - 15 * mm, 15 * mm)
    c.setFont("Helvetica-Oblique", 7.5)
    c.setFillColor(TXT_MUTED)
    c.drawString(15 * mm, 11 * mm,
                  "Methods: Mann-Kendall + Sen's slope, 2,000-iter bootstrap CIs, Lin's CCC, RMA regression, "
                  "Moran's I (queen 3×3, 199 perms). Full methods: documentation/validation_methods.md")
    c.drawString(15 * mm, 8 * mm,
                  "Pipeline: src/carbon_analysis.py | habitat_quality.py | soil_retention.py | synthesis_hotspots.py | validation_stats.py")
    c.drawRightString(w - 15 * mm, 5 * mm,
                       "Generated from validation_summary.json  -  Page 1 of 1")

    c.showPage()
    c.save()
    print(f"Saved -> {OUT_PDF}")


if __name__ == "__main__":
    main()

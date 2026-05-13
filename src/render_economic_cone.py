"""
Standalone re-renderer for forecast_economic_damages_2060.png.

Reads the existing results/economic_scenarios_2024_2060.csv in-place
(no re-running of the morphological simulation) and rewrites the chart
using the shared _render_economic_cone() in predictive_scenarios.py.

Run with: venv/bin/python src/render_economic_cone.py
"""
from pathlib import Path
import pandas as pd

from predictive_scenarios import _render_economic_cone

SCRIPT_DIR = Path(__file__).parent
RES_DIR    = SCRIPT_DIR.parent / 'results'
FIG_DIR    = SCRIPT_DIR.parent / 'figures'

CSV_IN  = RES_DIR / 'economic_scenarios_2024_2060.csv'
PNG_OUT = FIG_DIR / 'predictive' / 'forecast_economic_damages_2060.png'

if __name__ == '__main__':
    df = pd.read_csv(CSV_IN)
    # Normalise either snake_case or Title_Case column variants
    rename = {
        'year': 'Year',
        'scenario': 'Scenario',
        'economic_damage_usd': 'Economic_Damages_USD',
    }
    df = df.rename(columns=rename)
    _render_economic_cone(df, PNG_OUT)
    print(f"✅  Rewrote {PNG_OUT}")

from pathlib import Path
import nbformat as nbf

nb = nbf.v4.new_notebook()
cells = []
def md(text): cells.append(nbf.v4.new_markdown_cell(text))
def code(text): cells.append(nbf.v4.new_code_cell(text.strip()))
md('''# Denali basal ECM: basic AC / DC plots

One figure per section, with **AC on the left** and **DC on the right**, following the COLDEX basic track-map layout used for the Allan Hills cores. Depth increases downward; the horizontal axis is the measured lateral track position (mm).

Uses cleaned data and `ECM.norm_outside()` to match the two outer tracks to the inner-track mean. **No `norm_all()` normalization is applied**, so the colorbars retain current units (A). The default 10 mm end trim and 10 mm median smoothing follow the COLDEX notebook. Each panel uses its own 5th–95th percentile color limits from unflagged, trimmed measurements; colors are therefore not directly comparable between sections or AC/DC.

The key includes surface imperfections (black), tephra/layers (purple crosses), rock/debris (black circles), and **core breaks (blue diagonal hatching)**. Annotation spans come from the trimmed, unsmoothed flags, preserving narrow breaks. Overlapping annotations are retained. Missing DC measurements are labeled explicitly.

Run with the `allan_hills_exp` environment. Figures are saved to `figures/basic_plots/DEN13A/`.''')
code('''from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch, Rectangle
from IPython.display import display

# Supports launching from either the repository root or scripts/plotting.
ROOT = next((p for p in (Path.cwd(), *Path.cwd().parents)
             if (p / "scripts/core_scripts/ECMclass.py").is_file()), None)
if ROOT is None:
    raise FileNotFoundError("Launch this notebook from within denali_basal_ECM.")
sys.path.insert(0, str(ROOT / "scripts/core_scripts"))
from ECMclass import ECM
from core_names import source_core''')
md('## Settings')
code('''CORE = "DEN13A"
FACE = "t"
SECTIONS = None  # All sections; or, for example, ["224abc", "225cd"].
SMOOTH_MM = 10
TRIM_MM = 10
COLOR_PERCENTILES = (5, 95)
SHOW_FIGURES = True
DPI = 200

DATA_DIR = ROOT / "data/clean"
METADATA_PATH = ROOT / "data/metadata.csv"
FIGURE_DIR = ROOT / "figures/basic_plots" / CORE
FIGURE_DIR.mkdir(parents=True, exist_ok=True)
CMAP = plt.get_cmap("Spectral")
meta = pd.read_csv(METADATA_PATH, dtype={"section": str})

FLAG_STYLES = {
    "surface_imperfection": dict(facecolor="black", edgecolor="black", linewidth=0),
    "tephra_or_layer": dict(facecolor="none", edgecolor="purple", hatch="xx", linewidth=0),
    "rock_or_debris": dict(facecolor="none", edgecolor="black", hatch="oo", linewidth=0),
    "core_break": dict(facecolor="none", edgecolor="blue", hatch="///", linewidth=0),
}
FLAG_LABELS = ["Surface imperfection", "Tephra / layer", "Rock / debris", "Core break"]''')
md('## Loading and plotting helpers')
code('''def load_section(section):
    rows = meta.loc[(meta.core == CORE) & (meta.section == section) & (meta.face == FACE)]
    data = {}
    for mode in ("AC", "DC"):
        mode_rows = rows.loc[rows.ACorDC == mode]
        if mode_rows.empty:
            continue
        # Repeat metadata runs share one cleaned file; load that file once.
        filename = DATA_DIR / f"{source_core(CORE)}-{section}-{FACE}-{mode}.csv"
        if not filename.is_file():
            raise FileNotFoundError(f"Missing cleaned measurement: {filename}")
        item = ECM(CORE, section, FACE, mode, DATA_DIR, METADATA_PATH)
        item.rem_ends(TRIM_MM)
        if len(item.depth) == 0 or np.ptp(item.depth) <= 2 * SMOOTH_MM / 1000:
            raise ValueError(f"Section {section} ({mode}) is too short for the trim/smoothing settings")
        if len(item.y_vec) < 3:
            raise ValueError(f"Section {section} ({mode}) needs inner tracks for norm_outside")
        item.norm_outside()
        item.smooth(SMOOTH_MM)
        data[mode] = item
    if not data:
        raise ValueError(f"No measurements for {CORE}-{section}-{FACE}")
    return data


def sample_edges(depth):
    """Cell boundaries for sorted sample centers, clipped to measured coverage."""
    return np.r_[depth[0], (depth[:-1] + depth[1:]) / 2, depth[-1]]


def flag_runs(mask):
    """Return half-open slices of contiguous flagged samples."""
    changes = np.diff(np.r_[False, mask, False].astype(int))
    return zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1))


def color_norm(item):
    valid = (item.button == 0) & np.isfinite(item.meas)
    for flag in FLAG_STYLES:
        valid &= getattr(item, flag) == 0
    values = item.meas[valid]
    if not values.size:
        raise ValueError(f"No unflagged values for color limits: {item.section} {item.ACorDC}")
    lo, hi = np.percentile(values, COLOR_PERCENTILES)
    if lo == hi:
        pad = max(abs(lo) * 0.01, 1e-15)
        lo, hi = lo - pad, hi + pad
    return Normalize(lo, hi, clip=True)


def plot_tracks(ax, item, norm):
    tracks = np.sort(item.y_vec)
    spacing = np.min(np.diff(tracks))
    width = spacing - min(0.2, spacing * 0.05)
    for track in tracks:
        selected = np.flatnonzero(item.y_s == track)
        selected = selected[np.argsort(item.depth_s[selected])]
        depths = item.depth_s[selected]
        x0 = track - width / 2
        ax.pcolormesh([x0, x0 + width], sample_edges(depths),
                      item.meas_s[selected, None], cmap=CMAP, norm=norm,
                      shading="flat", rasterized=True)

        # Use original flags, not interpolated or smoothed flags, for annotations.
        raw = np.flatnonzero(item.y == track)
        raw = raw[np.argsort(item.depth[raw])]
        edges = sample_edges(item.depth[raw])
        for flag, style in FLAG_STYLES.items():
            for start, stop in flag_runs(getattr(item, flag)[raw] > 0):
                ax.add_patch(Rectangle((x0, edges[start]), width,
                                       edges[stop] - edges[start], **style))
    return tracks[0] - spacing / 2, tracks[-1] + spacing / 2


def plot_section(section):
    data = load_section(section)
    dmin = min(item.depth.min() for item in data.values())
    dmax = max(item.depth.max() for item in data.values())
    fig = plt.figure(figsize=(11, 7), layout="constrained")
    grid = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.65], height_ratios=[1, 0.05])
    left = fig.add_subplot(grid[0, 0])
    right = fig.add_subplot(grid[0, 1], sharey=left)
    limits = []
    for col, (mode, ax) in enumerate(zip(("AC", "DC"), (left, right))):
        ax.set_title(mode)
        ax.set_xlabel("Lateral track position (mm)")
        ax.set_ylabel("Depth (m)")
        ax.set_ylim(dmax, dmin)
        ax.ticklabel_format(axis="y", style="plain", useOffset=False)
        cax = fig.add_subplot(grid[1, col])
        if mode not in data:
            ax.text(0.5, 0.5, f"No {mode} data", ha="center", va="center", transform=ax.transAxes)
            cax.set_axis_off()
            continue
        norm = color_norm(data[mode])
        limits.append(plot_tracks(ax, data[mode], norm))
        fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=cax,
                     orientation="horizontal", label="Current (A)", extend="both")
    # Common lateral limits keep AC/DC track positions aligned.
    for ax in (left, right):
        ax.set_xlim(min(x[0] for x in limits), max(x[1] for x in limits))
    key = fig.add_subplot(grid[:, 2])
    key.set_axis_off()
    key.legend([Patch(**style) for style in FLAG_STYLES.values()], FLAG_LABELS,
               loc="center left", frameon=False, title="Key", handleheight=2, labelspacing=1.5)
    fig.suptitle(f"{CORE.upper()} · section {section} · outside tracks normalized")
    output = FIGURE_DIR / f"{CORE}-{section}-{FACE}.png"
    fig.savefig(output, dpi=DPI, bbox_inches="tight")
    if SHOW_FIGURES:
        display(fig)
    plt.close(fig)
    return output''')
md('## Generate figures\nSet `SECTIONS` above to a short list to plot selected sections. The default runs every section in the metadata for the selected core and face.')
code('''available = meta.loc[(meta.core == CORE) & (meta.face == FACE), "section"].drop_duplicates().tolist()
sections = available if SECTIONS is None else [str(section) for section in SECTIONS]
if not sections:
    raise ValueError(f"No sections selected for {CORE}, face {FACE}")
unknown = set(sections) - set(available)
if unknown:
    raise ValueError(f"Sections missing from metadata: {sorted(unknown)}")
outputs = []
for section in sections:
    output = plot_section(section)
    outputs.append(output)
    print(f"Saved {output.relative_to(ROOT)}")
print(f"Finished: {len(outputs)} figures")''')
nb['cells'] = cells
nb['metadata'] = {'kernelspec': {'display_name':'Python (allan_hills_exp)', 'language':'python', 'name':'python3'}, 'language_info': {'name':'python', 'version':'3.10'}}
nbf.validate(nb)
nbf.write(nb, 'scripts/plotting/basic_plots.ipynb')

from pathlib import Path
import nbformat as nbf
from nbclient import NotebookClient
HERE=Path(__file__).resolve().parent
md=nbf.v4.new_markdown_cell;code=nbf.v4.new_code_cell
nb=nbf.v4.new_notebook(cells=[
md('''# Simple manual-range AC dip estimates — current workflow

For every track, take the **raw maximum for a peak or raw minimum for a trough inside the exact inclusive manual depth range**. Fit a straight line through all those depth picks. No smoothing of picks, flag exclusion, trimming, prominence test, matching tolerance, directional expansion, outlier removal, angle cap, R² cutoff, residual cutoff, or stability test is applied. Each window is independent; overlapping windows are not deduplicated.

If multiple raw samples share the extremum, use the shallowest tied sample. Missing/nonfinite measurements cannot supply an extremum and are reported explicitly. A line requires two distinct lateral coordinates; a conventional residual-based interval requires at least three picks. These are mathematical requirements, not quality filters. Every selected window in this dataset has a pick on every track.

The slope is depth in mm per lateral mm; angle = arctan(slope). Ordinary least-squares standard errors and Student-t critical values (n−2 degrees of freedom) give the 95% slope interval, transformed to angle. There is no extra uncertainty floor. `ci_width_deg` is the **full** angular interval width, upper minus lower. The interval describes regression uncertainty conditional on the picked depths; it does not establish correct feature identity.'''),
code('''from pathlib import Path
import sys
import json
from IPython.display import display, Image
ROOT = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p / "data/metadata.csv").exists())
WORK = ROOT / "agentic_workspace/dip_angles"
sys.path.insert(0, str(WORK))
from simple_manual import run_simple
fits, picks = run_simple()
print(f"{len(fits)} fits, {len(picks)} raw track picks")'''),
code('''display(fits[["window_id", "core", "section", "kind", "n_tracks", "depth_m", "angle_deg",
              "ci_low_deg", "ci_high_deg", "ci_width_deg", "r2", "rmse_mm"]].round(4))'''),
md('''## Overall basal-core dip profile

All estimates share one continuous depth axis, covering the full measured basal interval. Uniform black markers show estimates with horizontal 95% error bars; no section divisions or offsets are applied.'''),
code('''display(Image(filename=str(ROOT / "figures/dip_angles/simple_manual/DEN13A-overall-dips.png")))'''),
md('''## Section figures

Every section has a figure in `figures/dip_angles/simple_manual/`. On the left, the map uses median smoothing and outer-track normalization **for display only**, with raw flag overlays. Picks and regression always use unmodified values from the cleaned AC files. Flags do not exclude picks. On the half-width right panel, uniform black markers show the estimates and horizontal bars show their 95% angle intervals. Axis limits are determined by the plotted results, with no dip cap.'''),
code('''SECTION = "217abc"
display(Image(filename=str(ROOT / f"figures/dip_angles/simple_manual/DEN13A-{SECTION}-t-simple.png")))'''),
md('''## Files and validation

`simple_results/fits.csv` contains every window fit. `picks.csv` records the raw source row, current, depth, tie count, flag status, and residual for every track. `input_windows.csv` snapshots the unchanged manual ranges; `input_sha256.json` fingerprints inputs. Superseded threshold-based workflows have been removed.

`validate_simple.py` checks every actual pick against the raw within-window extremum, verifies all tracks were included, tests flagged/boundary picks and ties, and independently checks the regression and interval calculation.'''),
code('''from validate_simple import validate
display(validate())''')])
nb.metadata={'kernelspec':{'display_name':'Python (allan_hills_exp)','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.10'}}
path=HERE/'simple_manual_dips.ipynb'
nbf.write(nb,path)
NotebookClient(nb,timeout=300,kernel_name='python3',resources={'metadata':{'path':str(HERE)}}).execute()
nbf.write(nb,path)

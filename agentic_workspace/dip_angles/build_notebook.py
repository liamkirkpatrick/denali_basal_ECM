from pathlib import Path
import nbformat as nbf

HERE=Path(__file__).resolve().parent
nb=nbf.v4.new_notebook()
md=nbf.v4.new_markdown_cell
code=nbf.v4.new_code_cell
nb.cells=[
md('''# Denali AC ECM: individual-feature dip estimates

Exploratory workflow kept in `agentic_workspace/dip_angles/`. Uses **cleaned AC only**, with no inter-track signal correlation. Every section gets a figure in `figures/dip_angles/dic1/`, including sections with no accepted features.

Fit depth against lateral track position: $z=a+bx$, using **mm for both coordinates**, then $\\theta=\\arctan(b)$. Positive angles mean increasing depth with increasing track position. These are **signed apparent dips along the measured face**, not reconstructed 3D dip angles. Depth increases downward in every plot.

Run with the existing **allan_hills_exp** environment. No new packages are required.'''),
code('''from pathlib import Path
import sys
import os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/denali-mpl")
import pandas as pd
from IPython.display import display, Image, Markdown
ROOT = next(p for p in (Path.cwd(), *Path.cwd().parents)
            if (p / "scripts/core_scripts/ECMclass.py").exists())
WORK = ROOT / "agentic_workspace/dip_angles"
sys.path.insert(0, str(WORK))
from dip_workflow import Settings, run_all
print("Python:", sys.executable)'''),
md('''## Detection and exclusion rules

1. Preserve all measured tracks. Use their common depth coverage, trim 10 mm at each end, and interpolate to a 0.25 mm grid. Both bracketing raw samples must be finite and have **every flag equal to zero**. Reject gaps exceeding 0.5 mm. Take the union of excluded depths across **all tracks**, including tracks not eventually picked. Thus an entirely flagged track excludes its whole section.
2. Process each uninterrupted clean interval separately. Detect peaks and troughs after Gaussian smoothing with $\\sigma=1$ mm. Require prominence greater than both six times the robust first-difference noise estimate and 15% of that track's clean 5th–95th percentile amplitude range. Prominence is local, within 25 mm; feature widths at half-prominence must be 1–15 mm. Broad extrema are intentionally omitted.
3. Keep candidates at least 17.5 mm from any clean-interval boundary, so smoothing and local prominence estimation cannot touch flagged measurements. Require a corresponding extremum within 1 mm at $\\sigma=0.75$ and 1.25 mm as a smoothing-stability check.
4. Form same-kind feature hypotheses from distinct-track pairs, bounded by $|\\theta|\\leq30°$. Match individual extrema within 1.25 mm of the proposed line. Reject ambiguous matches with multiple candidates on any track. Require at least four tracks, with no minimum lateral span or track fraction. Picks must remain in the seed's 30° search cone, allowing one grid step for discretization. All picks must belong to the same globally clean interval.
5. Fit ordinary least squares without removing outliers. Require RMS residual $\\leq0.6$ mm, maximum residual $\\leq1.2$ mm, leave-one-track-out angle changes $\\leq2°$. Require $R^2\\geq0.95$, except for nearly flat features with total depth range $\\leq1$ mm, where $R^2$ is uninformative.
6. Prefer hypotheses with more tracks, then narrower intervals. Never reuse a pick, and suppress same-kind duplicate features within 5 mm of mean depth. Peaks and troughs from the same waveform can still be statistically dependent.

The 2° interval cutoff and lateral-coverage requirements are disabled; intervals and span remain in the diagnostics. The remaining defaults prioritize a reviewable set. An empty section is **not** evidence of no dipping layers.'''),
code('''cfg = Settings()
SECTIONS = None  # All 16 AC sections, or e.g. ["217abc", "219b"].
# A subset run replaces result tables with that subset; existing figure files remain.
summary, fits = run_all(cfg, sections=SECTIONS)
display(summary)'''),
md('''## Fit precision and limitations

The slope standard error uses the larger of residual variance and a pick-error floor. The floor is the larger of 0.25 mm and the largest smoothing-induced candidate shift for that feature. A Student-t interval with $n-2$ degrees of freedom is calculated for the slope, then transformed with arctangent. Tables also report slope, $R^2$, RMS and maximum residuals, lateral span, track count, prominence/noise ratio, and leave-one-out sensitivity.

**These intervals are conditional on the automatically selected picks.** They assume independent track errors and do not account for feature-search selection, mistaken feature identity, correlated measurement errors, instrument alignment, or core orientation. They are not calibrated probabilities that the layer identification is correct. The smoothing check and error floor are practical safeguards, not a measurement-error calibration. No section-level average is inferred from these selected features.

Both point size and color represent interval precision (larger points and smaller colorbar values indicate tighter intervals). Horizontal bars show the conditional 95% interval. Up/down triangles indicate peaks/troughs; numerical labels link each estimate to its picks and fitted line. Gray bands in the right panel show depths excluded across tracks.

The left panel follows `basic_plots`: track strips, Spectral colors, outer-track mean scaling, 10 mm median display smoothing, and raw flag overlays. Detection uses the finer smoothing above on unscaled currents. Display normalization cannot move extrema. The display uses the common trimmed grid and color limits from globally clean data, so it is not pixel-identical to `basic_plots`.'''),
code('''columns = ["section", "feature_id", "kind", "depth_m", "angle_deg", "ci_low_deg", "ci_high_deg",
           "slope_mm_per_mm", "r2", "rmse_mm", "n_tracks", "loo_change_deg"]
if len(fits):
    display(fits[columns].round(4))
else:
    print("No fits passed selection.")'''),
md('''## Review a section

Change `REVIEW_SECTION` to inspect another section. All PNGs are saved regardless of whether they have retained features. Inspect peak identity before using the estimates scientifically; manual picking can be added later without changing the fitting stage.'''),
code('''REVIEW_SECTION = "217abc"
display(Image(filename=str(ROOT / f"figures/dip_angles/dic1/dic1-{REVIEW_SECTION}-t-AC-dip.png")))'''),
md('''## Outputs and validation

`results/fits.csv` contains accepted fits, `picks.csv` contains their individual track picks and residuals, and `candidates.csv` contains all extrema passing detection. `proposals.csv` records sufficiently populated, unambiguous hypotheses and their fit rejection/duplicate reasons; it does not enumerate every failed endpoint pair. `section_summary.csv`, `settings.json`, and `input_sha256.json` record coverage, thresholds, and input fingerprints. Candidate IDs are local to each section.

The synthetic checks recover isolated known-angle peaks, including a flat feature, and reject a feature above 30°, a flagged feature, and a feature intersecting a missing-data gap. Real-data checks verify every pick's full detection neighborhood is clean, no picks are reused, and all 16 figures exist. These checks establish implementation behavior, not a real-data false-discovery rate. Run validation after the default full-data run.'''),
code('''if SECTIONS is None:
    from validate_workflow import validate
    display(validate())''')]
nb.metadata={'kernelspec':{'display_name':'Python (allan_hills_exp)','language':'python','name':'python3'},'language_info':{'name':'python','version':'3.10'}}
nbf.write(nb,HERE/'dip_angles.ipynb')
if __name__=='__main__':
    from nbclient import NotebookClient
    NotebookClient(nb,timeout=300,kernel_name='python3',resources={'metadata':{'path':str(HERE)}}).execute()
    nbf.write(nb,HERE/'dip_angles.ipynb')

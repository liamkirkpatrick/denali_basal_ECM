# Exploratory Denali dip workflow

Use the existing `allan_hills_exp` environment. No additional packages are needed.

## Manual feature-range GUI

From the repository root:

```bash
conda run --no-capture-output -n allan_hills_exp python agentic_workspace/dip_angles/feature_picker.py
```

This opens **DC only** by default, using the button editor's map style: 10 mm median smoothing, Spectral colors, reversed lateral axis, and depth increasing downward. All cleaning flags are shaded black. The nine sections with DC data are available; the seven without DC are listed as skipped. `--mode AC` is available if AC reference picking is wanted later. There is no automatic fallback between modes.

- Hold **P**, then drag vertically over a peak search window.
- Hold **T**, then drag vertically over a trough search window.
- **N / B** or the Next / Previous buttons change section.
- **U** or Undo removes the last saved range in the current section/mode, including after reopening.
- **W** or Save explicitly saves; every selection and undo already autosaves.
- **Q** or Close exits. Saved ranges reappear next time.
- Toolbar zoom/pan works normally, but turn it off before selecting ranges. P is reserved for peak selection, not pan.

The CSV is `agentic_workspace/dip_angles/manual_feature_ranges.csv`, created on the first selection/save. Do not run multiple picker instances against the same CSV. The file contains one row per **search window**, not final per-track extrema:

| Field | Meaning |
|---|---|
| selection_id | Stable unique selection ID |
| core, section, face | Measurement identifiers |
| source_mode | DC or AC reference map used to select the range |
| kind | peak or trough in that reference signal |
| depth_min_m, depth_max_m | Ordered absolute depth bounds, metres |

Selections can overlap flags: these are manual search requests, not accepted dip observations. Subsequent AC fitting must still reject flagged samples. DC feature polarity/position need not match AC, so DC windows require AC review before being used for dip estimates. The angle pipeline does not yet consume these windows; that comes after manual selection.

`button_edits_v3.py` reloads this CSV for each section and overlays the matching windows on both AC and DC panels, labeled with the source mode. They are read-only reference bands; they do not set defect flags. In that editor, **T retains its existing tephra-marking behavior**. To see edits saved after a section was opened, reopen that section.

Optional commands:

```bash
# Review only specific DC sections:
conda run --no-capture-output -n allan_hills_exp python agentic_workspace/dip_angles/feature_picker.py --sections 222a 223e

# Render first section without opening a desktop window:
conda run -n allan_hills_exp python agentic_workspace/dip_angles/feature_picker.py --preview
```

GUI snapshots are saved to `figures/dip_angles/manual_ranges/`. Custom `--csv PATH` files are supported by the picker; the button editor reads the default CSV only.

## Automatic dip calculations

Open `dip_angles.ipynb` using `allan_hills_exp`. The notebook explains the algorithm and regenerates section figures and result tables. Supporting implementation: `dip_workflow.py`.

Current relaxed settings disable the 2° confidence-interval cutoff, minimum lateral span, and minimum track fraction. Four picks are still required. Prominence, flag exclusion, residual, R², smoothing stability, and leave-one-track-out stability criteria remain. Angle intervals are still reported; their color scale adapts by section. These are conditional precision intervals, not confidence in feature identity.

`results/` contains fits, individual picks, detected candidates, eligible/rejected fit hypotheses, section summaries, input hashes, and settings. Generated figures are under `figures/dip_angles/dic1/` and ignored by Git. No source data are modified.

## Validation

```bash
conda run -n allan_hills_exp python agentic_workspace/dip_angles/validate_workflow.py
conda run -n allan_hills_exp python agentic_workspace/dip_angles/validate_picker.py
```

The picker checks simulate actual Matplotlib mouse/key events, exercise autosave/reload/undo, and verify reference overlays in the button editor. They use temporary CSVs, never the user selection file. The UI validation screenshot is named `figures/dip_angles/manual_ranges/validation-example.png`. Native desktop interaction still needs the user's visual run.

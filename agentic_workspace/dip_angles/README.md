# Denali manual-range dip workflow

Use `simple_manual_dips.ipynb` with the existing `allan_hills_exp` environment. This is the retained dip-analysis workflow.

For every track, take the raw maximum (peak) or minimum (trough) inside each exact inclusive manual range. Fit a line using every track's pick. There are no quality filters, flag exclusions, smoothing of picks, trimming, outlier removal, angle limits, or window expansion. Tied extrema use the shallowest sample. Nonfinite or missing measurements are reported. A line needs two distinct track positions; a residual-based confidence interval needs three.

Depth increases downward. The signed apparent dip is `atan(depth slope)`, fitting depth and lateral position in mm. The 95% Student-t slope interval is transformed to angle. It has no additional uncertainty floor. These intervals describe regression uncertainty conditional on the picks.

## Manual range picker

From the repository root:

```bash
conda run --no-capture-output -n allan_hills_exp python agentic_workspace/dip_angles/feature_picker.py
```

The picker loads cleaned AC ECM only. Hold **P** or **T** and drag vertically to select a peak or trough window. **N/B** navigates sections, **U** undoes the last range in the current section, **W** saves, and **Q** closes. Every selection/undo autosaves to `manual_feature_ranges.csv`. Saved ranges reappear on reopening. Disable toolbar pan/zoom before selecting. Do not run multiple pickers against the same CSV.

Optional `--sections 214abcdefg 217abc` limits the displayed sections. `--preview` saves a PNG without opening a desktop window. `--csv PATH` uses another selection file; analysis and the button editor read the default file.

The GUI follows `basic_plots` display scaling and shows raw flag annotations. `button_edits_v3.py` reads saved windows as reference overlays without changing defect flags. In that editor, T retains its tephra-marking meaning.

## Outputs

- `simple_results/fits.csv`: every fit, including slope, angle, 95% interval bounds/full width, R², and RMS residual.
- `simple_results/picks.csv`: each pick's raw source row, current, depth, flag status, tie count, and residual.
- `simple_results/input_windows.csv` and `input_sha256.json`: input snapshot and fingerprints.
- `figures/dip_angles/simple_manual/`: one figure per section and `DEN13A-overall-dips.png`, with all estimates on one continuous depth axis.
- `figures/dip_angles/manual_ranges/`: picker snapshots.

AC maps use median smoothing and outer-track normalization **only for display**. They retain a current colorbar and flag legend. Angle estimates use uniform black markers with horizontal 95% error bars; there is no uncertainty color scale. Peak/trough triangles distinguish feature types.

`simple_manual.py` contains the analysis and plotting functions. `build_simple_notebook.py` regenerates and executes the notebook. Validation:

```bash
conda run -n allan_hills_exp python agentic_workspace/dip_angles/validate_simple.py
conda run -n allan_hills_exp python agentic_workspace/dip_angles/validate_picker.py
```

The first checks every actual pick against the raw extremum, complete track coverage, and independent slope/interval calculations. The second checks GUI interactions, persistence, and button-editor overlays using temporary selection files.

"""Build and execute the reproducible Denali oneline notebook."""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
cells = []
def markdown(text):
    cells.append(dict(cell_type='markdown', metadata={}, source=text.splitlines(True)))
def code(text):
    cells.append(dict(cell_type='code', metadata={}, source=text.splitlines(True), execution_count=None, outputs=[]))

markdown('''# Denali basal ECM: 1 cm oneline averages

Run in `allan_hills_exp`. All cleaned files are read once; every finite sample with all flags equal to zero contributes equally. Flags: button, surface imperfection, tephra/layer, rock/debris, core break. No trimming, smoothing, track removal, or normalization is applied.

Following COLDEX, centers are spaced at 0.01 m with half-open bins [center − 0.005, center + 0.005). AC and DC are separate. Unlike the reference's plotted medians, this figure shows **arithmetic means**, as requested. Missing bins remain blank in Excel and break plotted lines. Standard deviations describe within-bin population scatter (ddof=0), not uncertainty on the mean. Excel retains COLDEX column names, with depth in m and currents in A.
''')
code('''from pathlib import Path
import json
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = next(p for p in (Path.cwd(), *Path.cwd().parents) if (p / 'data/metadata.csv').exists())
WORK = ROOT / 'agentic_workspace/oneline'
FIG = ROOT / 'figures/average_oneline/DEN13A/DEN13A_oneline_1cm.png'
XLSX = ROOT / 'data/processed/DEN13A_oneline_1cm.xlsx'
FIG.parent.mkdir(parents=True, exist_ok=True)
XLSX.parent.mkdir(parents=True, exist_ok=True)
FLAGS = ['button', 'surface_imperfection', 'tephra_or_layer', 'rock_or_debris', 'core_break']
files = sorted((ROOT / 'data/clean').glob('dic1-*.csv'))
assert files, 'No cleaned data found'
frames, sources = [], []
for path in files:
    raw = pd.read_csv(path).rename(columns={'Button': 'button'})
    mode = path.stem.rsplit('-', 1)[1]
    assert mode in ('AC', 'DC')
    valid = raw[FLAGS].eq(0).all(axis=1) & np.isfinite(raw['meas']) & np.isfinite(raw['True_depth(m)'])
    frames.append(raw.assign(mode=mode, valid=valid))
    sources.append([str(path.relative_to(ROOT)), mode, len(raw), int(valid.sum()), int((~valid).sum())])
raw = pd.concat(frames, ignore_index=True)
# Integer bin IDs avoid cumulative floating-point drift; snap only roundoff at boundaries.
scaled = raw['True_depth(m)'].to_numpy() * 100 + 0.5
scaled = np.where(np.abs(scaled - np.rint(scaled)) < 1e-9, np.rint(scaled), scaled)
raw['bin'] = np.floor(scaled).astype(int)
ids = np.arange(raw['bin'].min(), raw['bin'].max() + 1)
oneline = pd.DataFrame({'depth': ids / 100}, index=ids)
for mode in ('AC', 'DC'):
    selected = raw.loc[raw.valid & raw['mode'].eq(mode)]
    grouped = selected.groupby('bin')['meas']
    oneline[mode] = grouped.mean()
    oneline[mode + '_median'] = grouped.median()
    oneline[mode + '_std'] = grouped.std(ddof=0)
    oneline[mode + '_n'] = grouped.size().reindex(ids, fill_value=0).astype(int)
    assert oneline[mode + '_n'].sum() == len(selected)
    # Independent direct-mask checks at first, middle, and last populated bins.
    populated = grouped.size().index.to_numpy()
    for b in populated[[0, len(populated)//2, -1]]:
        vals = selected.loc[selected['bin'].eq(b), 'meas'].to_numpy()
        assert np.isclose(oneline.loc[b, mode], np.mean(vals), rtol=1e-12, atol=1e-20)
oneline['bin_top_m'] = (ids - 0.5) / 100
oneline['bin_bottom_m'] = (ids + 0.5) / 100
oneline = oneline.reset_index(drop=True)
oneline.to_csv(WORK / 'oneline.csv', index=False)
print(oneline[['AC_n', 'DC_n']].sum())
print('Populated bins:', {m: int(oneline[m].notna().sum()) for m in ('AC', 'DC')})
''')
code('''fig, axes = plt.subplots(1, 2, figsize=(8, 11), sharey=True, layout='constrained')
for ax, mode, color in zip(axes, ('AC', 'DC'), ('#176B91', '#B35C28')):
    ax.plot(oneline[mode], oneline.depth, color=color, linewidth=1)
    ax.set_title(mode + ' ECM')
    ax.set_xlabel('Mean current (A)')
    ax.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax.grid(alpha=0.22)
    ax.set_ylim(oneline.bin_bottom_m.max(), oneline.bin_top_m.min())
    ax.spines[['top', 'right']].set_visible(False)
axes[0].set_ylabel('Depth (m)')
fig.suptitle('Denali basal ECM · 1 cm means', fontsize=15)
fig.savefig(FIG, dpi=250, bbox_inches='tight')
plt.show()
plt.close(fig)
print(FIG)
''')
markdown('''## Export Excel
The first sheet follows COLDEX column naming (`depth`, `AC`, `AC_median`, `AC_std`, and corresponding DC columns), with sample counts and explicit bin boundaries added. Methods and per-file inclusion counts are recorded on a second sheet. The exporter uses the bundled spreadsheet runtime; no packages are installed in the conda environment.
''')
code('''methods = [
    ['Depth / current units', 'm / A'],
    ['Statistic', 'Arithmetic mean of all retained samples; equal weight per sample.'],
    ['Bins', 'Width 0.01 m; centers at hundredth meters; [center - 0.005, center + 0.005).'],
    ['Exclusions', 'Nonfinite depth/current or any nonzero/missing button, surface_imperfection, tephra_or_layer, rock_or_debris, core_break.'],
    ['Preprocessing', 'Cleaned files, each loaded once. No end trimming, smoothing, track removal, or normalization.'],
    ['Missing data', 'Blank statistics and n=0; no interpolation across gaps.'],
    ['Scatter', 'AC_std and DC_std are population standard deviations (ddof=0), not standard errors.'],
    ['Reference', 'COLDEX make_oneline.ipynb: 1 cm centered bins, AC/DC panels, mean/median/std columns.'],
    ['Plot', 'Arithmetic means on linear axes, preserving signed currents; depth increases downward.'],
]
payload = {'columns': list(oneline.columns), 'rows': json.loads(oneline.to_json(orient='values', double_precision=15)), 'methods': methods, 'sources': sources}
# to_json decimal precision can round tiny currents; retain native floats via object replacement instead.
payload['rows'] = oneline.astype(object).where(pd.notna(oneline), None).values.tolist()
(WORK / 'workbook_input.json').write_text(json.dumps(payload, allow_nan=False))
node = Path.home() / '.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node'
subprocess.run([str(node), str(WORK / 'export_workbook.mjs')], check=True, cwd=ROOT)
print(XLSX)
''')
nb = dict(cells=cells, metadata={'kernelspec': {'display_name':'Python (allan_hills_exp)', 'language':'python', 'name':'allan_hills_exp'}, 'language_info': {'name':'python','version':'3.10'}}, nbformat=4, nbformat_minor=5)
for i,c in enumerate(cells):
    c['id'] = f'oneline-{i}'
(HERE / 'denali_oneline.ipynb').write_text(json.dumps(nb, indent=1))
if __name__ == '__main__':
    ns = {}
    for c in cells:
        if c['cell_type'] == 'code':
            exec(''.join(c['source']), ns)

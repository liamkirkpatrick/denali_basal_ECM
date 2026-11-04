import ast
import sys
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector as RealSpanSelector
from matplotlib.backend_bases import KeyEvent

root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(root / 'scripts/core_scripts'))
from ECMclass import ECM
work = root / 'agentic_workspace/editor_validation'
meta = pd.read_csv(root / 'data/metadata.csv')
# Verify every metadata entry resolves in the flat raw directory.
for _, row in meta.iterrows():
    obj = ECM(row['core'], row['section'], row['face'], row['ACorDC'], root / 'data/not_cleaned', root / 'data/metadata.csv')
    assert len(obj.core_break) == len(obj.depth)
    assert not obj.core_break.any()
print(f'Loaded all {len(meta)} raw files from flat directory')

# Exercise real callbacks with small synthetic tracks, without editing project data.
depth = np.tile(np.linspace(200, 200.1, 101), 2)
df = pd.DataFrame({'True_depth(m)': depth, 'Y_dimension(mm)': np.repeat([0, 10], 101), 'meas': np.arange(202) + 1, 'button': 0})
df.to_csv(work / 'test-1a-t-AC.csv', index=False)
pd.DataFrame([dict(core='test', section='1a', face='t', ACorDC='AC', time='', Y_left=0, Y_right=10)]).to_csv(work / 'metadata.csv', index=False)
d = ECM('test', '1a', 't', 'AC', work, work / 'metadata.csv')
selectors = []
def selector(*args, **kwargs):
    value = RealSpanSelector(*args, **kwargs)
    selectors.append(value)
    return value

def show():
    fig = plt.gcf()
    def key(k, release=False):
        event = 'key_release_event' if release else 'key_press_event'
        fig.canvas.callbacks.process(event, KeyEvent(event, fig.canvas, key=k))
    key('c')
    selectors[0].onselect(200.02, 200.04)
    assert 'core_break' in fig._suptitle.get_text()
    key('c', True)
    key('d')
    key('c')
    selectors[0].onselect(200.06, 200.08)
    key('c', True)
    key('n')
    selectors[0].onselect(200.01, 200.02)
    key('b')
    assert fig.axes[0].yaxis_inverted()
    key('w')

source = root / 'scripts/core_scripts/button_edits_v3.py'
tree = ast.parse(source.read_text())
ns = dict(np=np, pd=pd, plt=plt, Path=Path, sys=sys, SpanSelector=selector, window=10, entry_lookup={}, path_to_output_data=work / 'clean', project_root=root, cbar_line=matplotlib.colormaps['coolwarm'], cbar_2d=matplotlib.colormaps['Spectral'])
exec(compile(ast.Module(body=[node for node in tree.body if isinstance(node, ast.FunctionDef)], type_ignores=[]), str(source), 'exec'), ns)
plt.show = show
ns['process_ecm']({'ecm': d, 'note': ''}, 1, 1)
saved = pd.read_csv(work / 'clean/test-1a-t-AC.csv')
expected = (saved['Y_dimension(mm)'] == 0) & saved['True_depth(m)'].between(200.06, 200.08)
assert np.array_equal(saved.core_break, expected.astype(int))
assert saved.surface_imperfection.sum() > 0
assert not saved.rock_or_debris.any() and not saved.tephra_or_layer.any()
loaded = ECM('test', '1a', 't', 'AC', work / 'clean', work / 'metadata.csv')
assert np.array_equal(loaded.core_break, saved.core_break)
loaded.smooth(10)
assert loaded.core_break_s.max() == 1
loaded.rem_ends(5)
assert len(loaded.core_break) == len(loaded.depth)
assert len(loaded.core_break_s) == len(loaded.depth_s)
print('PASS: mark, release, undo, track navigation, flat save, reload, smoothing, clipping, depth direction')

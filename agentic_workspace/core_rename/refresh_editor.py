"""Re-render editor snapshots from saved flags without writing measurement data."""
from pathlib import Path
import sys, ast, re, os
os.environ['MPLCONFIGDIR']='/tmp/denali-rename-mpl'
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector
from matplotlib.backend_bases import KeyEvent
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/'scripts/core_scripts'),str(ROOT/'agentic_workspace/dip_angles')]
from ECMclass import ECM
from feature_ranges import read_ranges, section_ranges, overlay_ranges
source=ROOT/'scripts/core_scripts/button_edits_v3.py'
namespace=dict(np=np,pd=pd,plt=plt,Path=Path,sys=sys,re=re,SpanSelector=SpanSelector,
    window=10,entry_lookup={},path_to_output_data=ROOT/'agentic_workspace/core_rename/unused',project_root=ROOT,
    cbar_line=matplotlib.colormaps['coolwarm'],cbar_2d=matplotlib.colormaps['Spectral'],
    read_ranges=read_ranges,section_ranges=section_ranges,overlay_ranges=overlay_ranges)
tree=ast.parse(source.read_text())
exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],type_ignores=[]),str(source),'exec'),namespace)
meta=pd.read_csv(ROOT/'data/metadata.csv')
entries=[]
for row in meta.drop_duplicates(['core','section','face','ACorDC'],keep='last').to_dict('records'):
    item=ECM(row['core'],row['section'],row['face'],row['ACorDC'],ROOT/'data/clean',ROOT/'data/metadata.csv')
    entry={'ecm':item,'note':row.get('note','')}
    entries.append(entry)
    namespace['entry_lookup'][(item.core,namespace['normalize_section_for_lookup'](item.section),item.face,item.ACorDC)]=entry
plt.show=lambda:None
for number,entry in enumerate(entries,1):
    namespace['process_ecm'](entry,number,len(entries))
    fig=plt.gcf(); item=entry['ecm']
    # Restore saved classifications into the GUI's selection state, then display its final track.
    funcs=[ref() for ref in fig.canvas.callbacks.callbacks['key_press_event'].values()]
    handler=next(f for f in funcs if getattr(f,'__name__','')=='on_key_press')
    state=dict(zip(handler.__code__.co_freevars,[c.cell_contents for c in handler.__closure__]))
    for y in state['y_values']:
        idx=np.flatnonzero(item.y==y); z=item.depth[idx]
        for mark in ['surface_imperfection','tephra_or_layer','rock_or_debris','core_break']:
            mask=getattr(item,mark)[idx]!=0
            edges=np.diff(np.r_[False,mask,False].astype(int))
            for a,b in zip(np.flatnonzero(edges==1),np.flatnonzero(edges==-1)):
                state['selections'][y].append({'vmin':float(min(z[a:b])),'vmax':float(max(z[a:b])),'mark':mark})
    for _ in state['y_values'][1:]: handler(KeyEvent('key_press_event',fig.canvas,key='n'))
    state['redraw_current_y_patches']()
    fig.savefig(ROOT/f'figures/button_edits/{item.core}-{item.section}-{item.face}-{item.ACorDC}.png',dpi=150)
    plt.close(fig)
print('Refreshed editor snapshots; no measurement writes.')

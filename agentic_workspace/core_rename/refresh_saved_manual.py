"""Redraw existing manual fits without rerunning or changing fit selection."""
from pathlib import Path
import sys, os, ast, json
os.environ['MPLCONFIGDIR']='/tmp/denali-rename-mpl'
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'agentic_workspace/dip_angles'))
import pandas as pd
from dip_workflow import Settings, load, plot_section
for suffix in ('', '_pad20mm'):
    folder=ROOT/'agentic_workspace/dip_angles'/('manual_results'+suffix)
    settings=json.loads((folder/'settings.json').read_text())
    manual=settings.pop('manual_detection')
    cfg=Settings(**settings)
    candidates=pd.read_csv(folder/'candidates.csv')
    fits=pd.read_csv(folder/'fits.csv')
    figures=ROOT/'figures/dip_angles'/folder.name.replace('manual_results','manual')
    for path in sorted((ROOT/'data/clean').glob('*-AC.csv')):
        section=path.stem.split('-')[1]
        tracks,grid,values,safe=load(path,cfg)
        selected=candidates[candidates.section.astype(str)==section]
        records=fits[fits.section.astype(str)==section].to_dict('records')
        for fit in records: fit['candidate_ids']=ast.literal_eval(fit['candidate_ids'])
        stem='DEN13A-'+path.stem.split('-',1)[1]
        plot_section(stem+f" · manual windows (+{manual.get('padding_mm_each_end',0):g} mm each end)",tracks,grid,values,safe,selected,records,figures/f'{stem}-manual-dip.png',cfg)
    print('Refreshed saved plots:',figures,flush=True)
from feature_picker import FeaturePicker
import matplotlib.pyplot as plt
paths=sorted((ROOT/'data/clean').glob('*-AC.csv'))
gui=FeaturePicker(paths)
for i in range(len(paths)):
    gui.index=i; gui.show_section(); gui.snapshot()
plt.close(gui.fig)
print('Refreshed manual snapshots',flush=True)

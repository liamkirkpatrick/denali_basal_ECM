"""Headless mouse/key interaction and CSV persistence checks; no user data writes."""
import tempfile
from pathlib import Path
import ast
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backend_bases import KeyEvent, MouseEvent
from feature_picker import FeaturePicker
from feature_ranges import ROOT, read_ranges, section_ranges, overlay_ranges


def validate():
    # Keep automated test selections out of ordinary GUI snapshots and user CSVs.
    original_snapshot=FeaturePicker.snapshot
    FeaturePicker.snapshot=lambda self: None
    paths=sorted((ROOT/'data/clean').glob('*-DC.csv'))[:2]
    with tempfile.TemporaryDirectory() as tmp:
        csv=Path(tmp)/'ranges.csv'
        gui=FeaturePicker(paths,csv)
        canvas=gui.fig.canvas
        def key(value,release=False):
            name='key_release_event' if release else 'key_press_event'
            canvas.callbacks.process(name,KeyEvent(name,canvas,key=value))
        def mouse(name,depth):
            x=np.mean(gui.ax.get_xlim())
            xp,yp=gui.ax.transData.transform((x,depth))
            canvas.callbacks.process(name,MouseEvent(name,canvas,xp,yp,button=1))
        def drag(kind,low,high,release_early=False):
            canvas.draw()
            if kind: key(kind)
            mouse('button_press_event',low)
            if release_early: key(kind,True)
            mouse('motion_notify_event',high)
            mouse('button_release_event',high)
            if kind and not release_early: key(kind,True)
        lo,hi=gui.bounds
        a,b=lo+.1*(hi-lo),lo+.2*(hi-lo)
        drag('p',a,b,True)
        saved=read_ranges(csv)
        assert len(saved)==1 and saved.iloc[0].kind=='peak'
        assert abs(saved.iloc[0].depth_min_m-a)<.001
        drag('t',b+.04,a+.04)  # Reverse drag is sorted.
        saved=read_ranges(csv)
        assert len(saved)==2 and saved.iloc[1].kind=='trough'
        assert (saved.depth_min_m<saved.depth_max_m).all()
        drag(None,a,b)
        assert len(read_ranges(csv))==2
        key('n'); assert gui.index==1
        key('b'); assert gui.index==0 and len(gui.current_rows())==2
        key('u'); assert len(read_ranges(csv))==1
        key('w')
        # Real reference overlays preserve both source-mode label and range bounds.
        ranges=section_ranges(read_ranges(csv),gui.core,gui.section,gui.face)
        fig,ax=plt.subplots(); artists=overlay_ranges(ax,ranges)
        assert len(artists)==4 and 'DC peak' in artists[-1].get_text()
        assert np.allclose(artists[1].get_ydata(),ranges.iloc[0].depth_min_m)
        plt.close(fig)
        # Save a clearly named UI validation image, not a user selection file.
        out=ROOT/'figures/dip_angles/manual_ranges/validation-example.png'
        gui.fig.savefig(out,dpi=120)
        plt.close(gui.fig)
        reopened=FeaturePicker(paths,csv)
        assert len(reopened.current_rows())==1
        reopened.undo()
        assert read_ranges(csv).empty
        plt.close(reopened.fig)
        # Execute just process_ecm with a safe namespace to check editor integration.
        import pandas as pd
        import sys
        from feature_picker import ECM
        from matplotlib.widgets import SpanSelector
        source=ROOT/'scripts/core_scripts/button_edits_v3.py'
        tree=ast.parse(source.read_text())
        namespace=dict(np=np,pd=pd,plt=plt,Path=Path,sys=sys,SpanSelector=SpanSelector,
                       window=10,entry_lookup={},path_to_output_data=Path(tmp),project_root=ROOT,
                       cbar_line=matplotlib.colormaps['coolwarm'],cbar_2d=matplotlib.colormaps['Spectral'],
                       read_ranges=lambda: ranges,section_ranges=section_ranges,overlay_ranges=overlay_ranges)
        exec(compile(ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef)],
                                type_ignores=[]),str(source),'exec'),namespace)
        core,section,face,mode=paths[0].stem.split('-')
        item=ECM(core,section,face,mode,ROOT/'data/clean',ROOT/'data/metadata.csv')
        original_show=plt.show
        plt.show=lambda: None
        namespace['process_ecm']({'ecm':item,'note':'GUI integration validation'},1,1)
        editor=plt.gcf()
        assert any('DC peak' in text.get_text() for text in editor.axes[0].texts)
        assert any('DC peak' in text.get_text() for text in editor.axes[3].texts)
        plt.close(editor); plt.show=original_show
    FeaturePicker.snapshot=original_snapshot
    print('PASS: P/T mouse drags, early key release, reverse drag, unmodified drag ignored, '
          'navigation, autosave/reload, undo, empty CSV, and button-editor reference overlays.')

if __name__=='__main__': validate()

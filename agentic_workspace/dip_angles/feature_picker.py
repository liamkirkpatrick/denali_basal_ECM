"""Native Matplotlib GUI: hold p/t and drag vertically to save feature windows.

Run: conda run --no-capture-output -n allan_hills_exp python \
agentic_workspace/dip_angles/feature_picker.py
"""
from pathlib import Path
import argparse
import sys
import uuid
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import SpanSelector, Button
from matplotlib.patches import Patch
from feature_ranges import (ROOT, DEFAULT_CSV, COLORS, read_ranges, write_ranges,
                            section_ranges, overlay_ranges)

sys.path.insert(0, str(ROOT/'scripts/core_scripts'))
from ECMclass import ECM


class FeaturePicker:
    def __init__(self, paths, csv_path=DEFAULT_CSV):
        self.paths = list(paths)
        if not self.paths:
            raise ValueError('No matching cleaned ECM files.')
        self.csv_path = Path(csv_path)
        self.rows = read_ranges(self.csv_path)
        self.index = 0
        self.held = set()
        self.drag_kind = None
        self.artists = []
        self.fig, self.ax = plt.subplots(figsize=(9,9))
        self.fig.subplots_adjust(left=.13,right=.96,bottom=.22,top=.88)
        self.status = self.fig.text(.13,.155,'',fontsize=10)
        self.fig.text(.13,.085,'Hold P + drag: peak   |   Hold T + drag: trough\n'
                      'N/B: next/previous section   U: undo last range   W: save   Q: close\n'
                      'Selections autosave. Release pan/zoom before selecting.',fontsize=10)
        self.buttons = []
        for i,(label,callback) in enumerate([
            ('Previous',lambda e:self.move(-1)), ('Next',lambda e:self.move(1)),
            ('Undo',self.undo), ('Save',self.save), ('Close',lambda e:plt.close(self.fig))]):
            button=Button(self.fig.add_axes([.13+i*.17,.025,.15,.04]),label)
            button.on_clicked(callback); self.buttons.append(button)
        self.fig.canvas.mpl_connect('key_press_event', self.key_press)
        self.fig.canvas.mpl_connect('key_release_event', self.key_release)
        # Capture the modifier at mouse-down, even if released before mouse-up.
        self.fig.canvas.mpl_connect('button_press_event', self.mouse_press)
        self.fig.canvas.mpl_connect('figure_leave_event', self.leave)
        self.selector = None
        self.show_section()

    def current_rows(self):
        return section_ranges(self.rows,self.core,self.section,self.face,self.mode)

    def show_section(self):
        self.held.clear(); self.drag_kind=None
        if self.selector is not None:
            self.selector.disconnect_events()
        self.ax.clear(); self.artists=[]
        self.core,self.section,self.face,self.mode = self.paths[self.index].stem.split('-')
        item=ECM(self.core,self.section,self.face,self.mode,self.paths[self.index].parent,ROOT/'data/metadata.csv')
        item.smooth(10)
        # Same top-down map as button_edits_v3's far-left subplot: Spectral,
        # 10-mm median smoothing, no outside-track scaling, lateral axis reversed.
        table=pd.DataFrame({'z':item.depth_s,'x':item.y_s,'v':item.meas_s}).pivot(index='z',columns='x',values='v')
        self.ax.pcolormesh(table.columns,table.index,table.to_numpy(),shading='auto',cmap='Spectral')
        # Include every raw cleaning flag; preserve narrow flagged intervals.
        bad=np.logical_or.reduce([getattr(item,f)!=0 for f in
            ['button','surface_imperfection','tephra_or_layer','rock_or_debris','core_break']])
        tracks=np.sort(item.y_vec)
        spacing=np.min(np.diff(tracks))
        from matplotlib.patches import Rectangle
        for x in tracks:
            idx=np.flatnonzero(item.y==x); idx=idx[np.argsort(item.depth[idx])]
            z=item.depth[idx]; edge=np.r_[z[0],(z[:-1]+z[1:])/2,z[-1]]
            changes=np.diff(np.r_[False,bad[idx],False].astype(int))
            for a,b in zip(np.flatnonzero(changes==1),np.flatnonzero(changes==-1)):
                self.ax.add_patch(Rectangle((x-spacing/2,edge[a]),spacing,edge[b]-edge[a],
                                           facecolor='black',alpha=.5,zorder=3))
        self.bounds=(float(item.depth.min()),float(item.depth.max()))
        self.ax.set(xlabel='Distance across core (mm)',ylabel='Depth (m)',
                    xlim=(tracks[-1]+spacing/2,tracks[0]-spacing/2),
                    ylim=(self.bounds[1],self.bounds[0]))
        self.ax.ticklabel_format(axis='y',style='plain',useOffset=False)
        self.ax.set_title(f'{self.core} · {self.section} · face {self.face} · {self.mode}\n'
                          f'Section {self.index+1}/{len(self.paths)} · select feature-search windows')
        self.fig.suptitle('ECM manual peak / trough range picker',fontsize=14)
        self.ax.legend(handles=[Patch(facecolor=COLORS[k],alpha=.35,label=k.title()) for k in COLORS]
                       +[Patch(facecolor='black',alpha=.5,label='Any cleaning flag')],
                       loc='upper right',fontsize=8)
        self.selector=SpanSelector(self.ax,self.select,'vertical',useblit=True,
                                   props=dict(facecolor='white',edgecolor='black',alpha=.3),
                                   button=1,minspan=.00025,interactive=False)
        self.redraw_ranges()
        self.snapshot()

    def redraw_ranges(self):
        for artist in self.artists: artist.remove()
        self.artists=overlay_ranges(self.ax,self.current_rows())
        self.message(f'{len(self.current_rows())} saved windows · {self.csv_path.name}')
        self.fig.canvas.draw_idle()

    def message(self,text):
        self.status.set_text(text); self.fig.canvas.draw_idle()

    def mouse_press(self,event):
        self.drag_kind=None
        toolbar=getattr(self.fig.canvas.manager,'toolbar',None)
        if event.inaxes is self.ax and event.button==1 and not (toolbar and toolbar.mode):
            if self.held=={'p'}: self.drag_kind='peak'
            elif self.held=={'t'}: self.drag_kind='trough'

    def select(self,low,high):
        kind=self.drag_kind; self.drag_kind=None
        if kind is None:
            self.selector.clear()
            self.message('Hold P (peak) or T (trough) while starting the drag.')
            return
        low,high=sorted([low,high])
        low=max(low,self.bounds[0]); high=min(high,self.bounds[1])
        if high-low < .00025: return
        row=dict(selection_id=uuid.uuid4().hex,core=self.core,section=self.section,
                 face=self.face,source_mode=self.mode,kind=kind,depth_min_m=low,depth_max_m=high)
        updated=(pd.DataFrame([row]) if self.rows.empty else
                 pd.concat([self.rows,pd.DataFrame([row])],ignore_index=True))
        write_ranges(updated,self.csv_path)
        self.rows=updated
        self.selector.clear(); self.redraw_ranges(); self.snapshot()

    def key_press(self,event):
        key=(event.key or '').lower()
        if key in ('p','t'):
            self.held.add(key)
            self.message(f'Hold {key.upper()} and drag vertically to mark a {"peak" if key=="p" else "trough"}.')
        elif key=='n': self.move(1)
        elif key=='b': self.move(-1)
        elif key in ('u','backspace'): self.undo()
        elif key=='w': self.save()
        elif key=='q': plt.close(self.fig)
        elif key=='escape': self.held.clear(); self.drag_kind=None; self.selector.clear()

    def key_release(self,event):
        self.held.discard((event.key or '').lower())

    def leave(self,event):
        self.held.clear()

    def move(self,delta):
        target=self.index+delta
        if 0<=target<len(self.paths):
            self.index=target; self.show_section()
        else: self.message('Already at the first/last available section.')

    def undo(self,event=None):
        selected=self.current_rows()
        if selected.empty:
            self.message('No saved range to undo in this section.'); return
        updated=self.rows.drop(selected.index[-1]).reset_index(drop=True)
        write_ranges(updated,self.csv_path); self.rows=updated
        self.redraw_ranges(); self.snapshot()

    def save(self,event=None):
        write_ranges(self.rows,self.csv_path)
        self.message(f'Saved {len(self.rows)} total windows to {self.csv_path.name}')
        self.snapshot()

    def snapshot(self):
        folder=ROOT/'figures/dip_angles/manual_ranges'; folder.mkdir(parents=True,exist_ok=True)
        self.fig.savefig(folder/f'{self.core}-{self.section}-{self.face}-{self.mode}-ranges.png',dpi=140)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core',default='dic1')
    parser.add_argument('--face',default='t')
    parser.add_argument('--mode',choices=['AC','DC'],default='DC')
    parser.add_argument('--sections',nargs='+',help='Optional exact section IDs, e.g. 222a 223e')
    parser.add_argument('--csv',type=Path,default=DEFAULT_CSV)
    parser.add_argument('--preview',action='store_true',help='Save first-section PNG without opening a window')
    args=parser.parse_args()
    if args.preview: plt.switch_backend('Agg')
    elif 'agg'==matplotlib.get_backend().lower():
        raise RuntimeError('No interactive backend. Run from a local desktop terminal using allan_hills_exp.')
    paths=sorted((ROOT/'data/clean').glob(f'{args.core}-*-{args.face}-{args.mode}.csv'))
    if args.sections:
        paths=[p for p in paths if p.stem.split('-')[1] in args.sections]
        absent=set(args.sections)-{p.stem.split('-')[1] for p in paths}
        if absent: print(f'No {args.mode} data for requested sections: {sorted(absent)}')
    else:
        all_sections={p.stem.split('-')[1] for p in (ROOT/'data/clean').glob(f'{args.core}-*-{args.face}-*.csv')}
        absent=all_sections-{p.stem.split('-')[1] for p in paths}
        print(f'Sections skipped (no {args.mode} data): {sorted(absent)}')
    # Matplotlib normally binds P to pan. Reserve all GUI shortcuts.
    reserved={'p','t','n','b','u','w','q','backspace'}
    for name in [key for key in matplotlib.rcParams if key.startswith('keymap.')]:
        matplotlib.rcParams[name]=[key for key in matplotlib.rcParams[name] if key not in reserved]
    gui=FeaturePicker(paths,args.csv)
    print(f'Loaded {len(paths)} {args.mode} sections. CSV: {args.csv}')
    print('Hold P/T and drag vertically. Autosave after each selection; U undoes; N/B changes section.')
    if args.preview: plt.close(gui.fig)
    else: plt.show()

if __name__=='__main__': main()

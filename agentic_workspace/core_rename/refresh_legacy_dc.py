"""Refresh the existing historical DC range snapshot without enabling DC picking."""
from pathlib import Path
import os, sys
os.environ['MPLBACKEND']='Agg'
os.environ['MPLCONFIGDIR']='/tmp/denali-rename-mpl'
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'agentic_workspace/dip_angles'))
from feature_picker import FeaturePicker
import matplotlib.pyplot as plt
snapshot=ROOT/'figures/dip_angles/manual_ranges/DEN13A-222a-t-DC-ranges.png'
if snapshot.exists():
    gui=FeaturePicker([ROOT/'data/clean/dic1-222a-t-AC.csv'])
    gui.paths=[ROOT/'data/clean/dic1-222a-t-DC.csv']
    gui.show_section(); gui.snapshot(); plt.close(gui.fig)

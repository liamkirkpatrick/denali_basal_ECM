"""Shared CSV contract for manually selected ECM feature-search windows.

These are reference windows, never cleaning/defect flags or fitted peak picks.
"""
from pathlib import Path
import os
import tempfile
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CSV = ROOT / 'agentic_workspace/dip_angles/manual_feature_ranges.csv'
COLUMNS = ['selection_id', 'core', 'section', 'face', 'source_mode', 'kind',
           'depth_min_m', 'depth_max_m']
COLORS = {'peak': '#e76f00', 'trough': '#0087a3'}


def read_ranges(path=DEFAULT_CSV):
    path = Path(path)
    if not path.exists():
        return pd.DataFrame(columns=COLUMNS)
    frame = pd.read_csv(path, dtype={c: str for c in COLUMNS[:6]})
    if set(COLUMNS) - set(frame.columns):
        raise ValueError(f'{path}: missing columns {set(COLUMNS)-set(frame.columns)}')
    frame = frame[COLUMNS].copy()
    for name in ['depth_min_m', 'depth_max_m']:
        frame[name] = pd.to_numeric(frame[name], errors='raise')
    if frame[COLUMNS[:6]].isna().any().any() or frame.selection_id.duplicated().any():
        raise ValueError(f'{path}: blank identifiers or duplicate selection IDs')
    if not frame.kind.isin(COLORS).all() or not frame.source_mode.isin(['AC','DC']).all():
        raise ValueError(f'{path}: invalid feature kind or source mode')
    if not np.isfinite(frame[['depth_min_m','depth_max_m']].to_numpy()).all():
        raise ValueError(f'{path}: non-finite depth')
    if not (frame.depth_min_m < frame.depth_max_m).all():
        raise ValueError(f'{path}: reversed or empty ranges')
    return frame


def write_ranges(frame, path=DEFAULT_CSV):
    """Atomically replace the CSV, so an interrupted save cannot truncate it."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, suffix='.csv', delete=False) as tmp:
        name = Path(tmp.name)
        frame.to_csv(tmp, columns=COLUMNS, index=False, float_format='%.9f')
    try:
        read_ranges(name)  # Validate before replacing the user's saved selections.
        os.replace(name, path)
    finally:
        name.unlink(missing_ok=True)


def section_ranges(frame, core, section, face, source_mode=None):
    mask = (frame.core==str(core)) & (frame.section==str(section)) & (frame.face==str(face))
    if source_mode is not None:
        mask &= frame.source_mode==source_mode
    return frame.loc[mask]


def overlay_ranges(ax, frame, labels=True):
    artists = []
    for row in frame.itertuples():
        color = COLORS[row.kind]
        artists.append(ax.axhspan(row.depth_min_m,row.depth_max_m,facecolor=color,
                                  edgecolor=color,alpha=.18,zorder=5))
        for depth in [row.depth_min_m,row.depth_max_m]:
            artists.append(ax.axhline(depth,color=color,lw=.9,ls='--',zorder=6))
        if labels:
            artists.append(ax.text(.01,(row.depth_min_m+row.depth_max_m)/2,
                                   f'{row.source_mode} {row.kind} {row.selection_id[:6]}',
                                   transform=ax.get_yaxis_transform(),color=color,fontsize=8,
                                   va='center',zorder=7,bbox=dict(facecolor='white',alpha=.8,edgecolor='none')))
    return artists

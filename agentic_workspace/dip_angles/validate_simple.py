"""Check direct raw extrema, flags/boundaries, and analytical OLS."""
from types import SimpleNamespace
import json
import numpy as np
import pandas as pd
from scipy.stats import linregress,t
from simple_manual import ROOT,WORK,pick_extrema,fit_line
from core_names import source_core


def validate():
    window=SimpleNamespace(depth_min_m=1.,depth_max_m=1.02,kind='peak')
    raw=pd.DataFrame({'True_depth(m)':[1,1.01,1.02]*3,'Y_dimension(mm)':np.repeat([0,5,10],3),
                      'meas':[10,1,0,0,1,10,10,10,0],'button':1})
    picks=pick_extrema(raw,window)
    assert picks.depth_m.tolist()==[1.,1.02,1.]  # Boundaries, flags, and shallowest tie retained.
    assert picks.flagged.all()
    fit=fit_line(picks); reference=linregress(picks.x_mm,picks.depth_m*1000)
    assert np.isclose(fit['slope_mm_per_mm'],reference.slope)
    expected=np.degrees(np.arctan([reference.slope-t.ppf(.975,1)*reference.stderr,
                                  reference.slope+t.ppf(.975,1)*reference.stderr]))
    assert np.allclose([fit['ci_low_deg'],fit['ci_high_deg']],expected)
    missing=pick_extrema(raw,SimpleNamespace(depth_min_m=2,depth_max_m=3,kind='peak'))
    assert fit_line(missing)['status']=='line_not_identifiable'
    windows=pd.read_csv(WORK/'input_windows.csv').set_index('selection_id')
    picks=pd.read_csv(WORK/'picks.csv');fits=pd.read_csv(WORK/'fits.csv')
    for (core,section,face),group in picks.groupby(['core','section','face']):
        raw=pd.read_csv(ROOT/f'data/clean/{source_core(core)}-{section}-{face}-AC.csv')
        for p in group.itertuples():
            w=windows.loc[p.selection_id]
            track=raw[(raw['Y_dimension(mm)']==p.x_mm)&raw['True_depth(m)'].between(w.depth_min_m,w.depth_max_m)]
            expected=track.meas.max() if w.kind=='peak' else track.meas.min()
            assert p.status=='picked'
            assert np.isclose(p.current_A,expected,rtol=1e-12,atol=0)
        for sid,pp in group.groupby('selection_id'):
            assert len(pp)==raw['Y_dimension(mm)'].nunique()
            ff=fits[fits.selection_id==sid].iloc[0]
            reference=linregress(pp.x_mm,pp.depth_m*1000)
            assert np.isclose(ff.slope_mm_per_mm,reference.slope,atol=1e-8)
    assert len(fits)==len(windows) and fits.status.eq('fitted').all()
    report=dict(windows=len(windows),fits=len(fits),picks=len(picks),
                every_track_included=True,raw_extrema_verified=True,flagged_and_boundary_picks_retained=True,
                slope_and_interval_verified=True)
    (WORK/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':print(validate())

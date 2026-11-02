"""Scientific checks: known inclinations, flagged features, gaps, and real-data invariants."""
import json
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd
from dip_workflow import ROOT, FLAGS, Settings, load, detect, match


def synthetic(angle, flag=False, gap=False):
    rng = np.random.default_rng(42)
    z = np.arange(0,250,.25)
    records = []
    for x in np.arange(0,50,5):
        center = 125+np.tan(np.radians(angle))*(x-22.5)
        v = np.exp(-.5*((z-center)/3)**2)+rng.normal(0,.003,len(z))
        part = pd.DataFrame({'True_depth(m)':z/1000,'Y_dimension(mm)':x,'meas':v})
        for f in FLAGS: part[f] = 0
        if flag and x==20: part.loc[(z>95)&(z<155),'core_break'] = 1
        if gap and x==20: part = part.loc[~((z>95)&(z<155))]
        records.append(part)
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp)/'synthetic.csv'; pd.concat(records).to_csv(path,index=False)
        cfg = Settings(); tracks,grid,v,safe = load(path,cfg)
        candidates = detect(tracks,grid,v,safe,cfg)
        fits,_ = match(candidates,tracks,cfg)
    return fits


def validate():
    report = {}
    for angle in [-20,0,15,29]:
        fits = synthetic(angle)
        assert len(fits)==1, (angle, fits)
        assert abs(fits[0]['angle_deg']-angle)<.5, (angle,fits)
        assert fits[0]['ci_low_deg'] <= angle <= fits[0]['ci_high_deg']
        report[f'recovered_{angle}_deg'] = fits[0]['angle_deg']
    assert not synthetic(35)
    assert not synthetic(15,flag=True)
    assert not synthetic(15,gap=True)
    report.update(reject_above_30=True, reject_flagged_feature=True, reject_gap=True)
    cfg=Settings()
    summary = pd.read_csv(ROOT/'agentic_workspace/dip_angles/results/section_summary.csv')
    fits = pd.read_csv(ROOT/'agentic_workspace/dip_angles/results/fits.csv')
    picks = pd.read_csv(ROOT/'agentic_workspace/dip_angles/results/picks.csv')
    for row in summary.itertuples():
        path = ROOT/f'data/clean/dic1-{row.section}-t-AC.csv'
        tracks,grid,values,safe=load(path,cfg)
        selected=picks[picks.section==row.section]
        for p in selected.itertuples():
            k=np.argmin(np.abs(grid-p.depth_mm))
            margin=int((cfg.prominence_window_mm/2+5*cfg.sigma_mm)/cfg.step_mm)
            assert safe[max(0,k-margin):k+margin+1].all()
        assert (ROOT/f'figures/dip_angles/dic1/dic1-{row.section}-t-AC-dip.png').exists()
    assert not picks.duplicated(['section','candidate_id']).any()
    assert (fits.angle_deg.abs()<=cfg.max_dip_deg).all()
    assert (fits.rmse_mm<=cfg.max_rmse_mm).all()
    if cfg.max_ci_half_deg is not None:
        assert (fits.ci_half_deg<=cfg.max_ci_half_deg).all()
    report.update(real_picks_clear_of_flags=True, no_reused_picks=True,section_figures=len(summary),accepted_fits=len(fits))
    out=ROOT/'agentic_workspace/dip_angles/results/validation.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__': print(validate())

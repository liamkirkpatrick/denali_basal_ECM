"""Conservative, feature-based AC ECM apparent dips. All distances internal: mm."""
from pathlib import Path
from dataclasses import dataclass, asdict
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import find_peaks
from scipy.stats import t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle, Patch

ROOT = Path(__file__).resolve().parents[2]
FLAGS = ['button', 'surface_imperfection', 'tephra_or_layer', 'rock_or_debris', 'core_break']
STYLES = dict(surface_imperfection=dict(facecolor='black', edgecolor='black'),
              tephra_or_layer=dict(facecolor='none', edgecolor='purple', hatch='xx'),
              rock_or_debris=dict(facecolor='none', edgecolor='black', hatch='oo'),
              core_break=dict(facecolor='none', edgecolor='blue', hatch='///'))

@dataclass
class Settings:
    step_mm: float = .25
    trim_mm: float = 10
    sigma_mm: float = 1
    prominence_window_mm: float = 25
    prominence_fraction: float = .15
    noise_multiplier: float = 6
    max_dip_deg: float = 30
    match_mm: float = 1.25
    min_tracks: int = 4
    min_coverage: float = 0
    min_span_fraction: float = 0
    max_rmse_mm: float = .6
    max_residual_mm: float = 1.2
    min_r2: float = .95
    flat_range_mm: float = 1
    max_ci_half_deg: float | None = None
    max_loo_change_deg: float = 2


def runs(mask):
    changes = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)))


def edges(z):
    return np.r_[z[0], (z[:-1] + z[1:]) / 2, z[-1]]


def load(path, cfg):
    frame = pd.read_csv(path).rename(columns={'Button': 'button'})
    missing = set(FLAGS) - set(frame.columns)
    if missing:
        raise ValueError(f'{path}: missing flag columns {missing}')
    tracks = []
    for x, part in frame.groupby('Y_dimension(mm)', sort=True):
        part = part.sort_values('True_depth(m)')
        z = part['True_depth(m)'].to_numpy() * 1000
        if np.any(np.diff(z) <= 0):
            raise ValueError(f'{path}: duplicate/non-increasing depths on track {x}')
        tracks.append(dict(x=x, z=z, v=part.meas.to_numpy(),
                           bad=part[FLAGS].ne(0).any(axis=1).to_numpy() | ~np.isfinite(part.meas.to_numpy()),
                           flags={f: part[f].ne(0).to_numpy() for f in FLAGS}))
    lo = max(tr['z'][0] for tr in tracks) + cfg.trim_mm
    hi = min(tr['z'][-1] for tr in tracks) - cfg.trim_mm
    grid = np.arange(lo, hi, cfg.step_mm)
    safe = np.ones(len(grid), dtype=bool)
    signals = []
    for tr in tracks:
        z, v = tr['z'], tr['v']
        # Both bracketing raw samples must be clean; never bridge a missing interval.
        j = np.clip(np.searchsorted(z, grid), 1, len(z)-1)
        safe &= ~(tr['bad'][j] | tr['bad'][j-1]) & (z[j]-z[j-1] <= 2*cfg.step_mm)
        signals.append(np.interp(grid, z, v))
    return tracks, grid, np.array(signals), safe


def detect(tracks, grid, values, safe, cfg):
    rows = []
    step = cfg.step_mm
    margin = cfg.prominence_window_mm / 2 + 4 * 1.25 * cfg.sigma_mm
    wlen = int(cfg.prominence_window_mm / step) | 1
    for ti, tr in enumerate(tracks):
        valid_values = values[ti, safe]
        if len(valid_values) < 3:
            continue
        amplitude = np.percentile(valid_values, 95) - np.percentile(valid_values, 5)
        for run_id, (a,b) in enumerate(runs(safe)):
            if (b-a)*step < 2*margin:
                continue
            raw = values[ti,a:b]
            noise = np.median(np.abs(np.diff(raw)-np.median(np.diff(raw)))) / (.67449*np.sqrt(2))
            threshold = max(cfg.prominence_fraction*amplitude, cfg.noise_multiplier*noise)
            if threshold <= 0:
                continue
            sm = gaussian_filter1d(raw, cfg.sigma_mm/step)
            variants = [gaussian_filter1d(raw, cfg.sigma_mm*s/step) for s in (.75,1.25)]
            for kind, sign in [('peak',1),('trough',-1)]:
                peaks, props = find_peaks(sign*sm, prominence=threshold, wlen=wlen,
                                         width=(1/step,15/step), distance=3/step)
                alt_peaks = [find_peaks(sign*v, prominence=threshold*.8, wlen=wlen)[0] for v in variants]
                for k, p in enumerate(peaks):
                    if p*step < margin or (b-a-1-p)*step < margin:
                        continue
                    shifts = [np.min(np.abs(q-p))*step if len(q) else np.inf for q in alt_peaks]
                    if max(shifts) > 1:
                        continue
                    rows.append(dict(candidate_id=len(rows), track=ti, x_mm=tr['x'],
                                     depth_mm=grid[a+p], kind=kind, run_id=run_id,
                                     prominence=props['prominences'][k],
                                     prominence_snr=props['prominences'][k]/max(noise,1e-30),
                                     width_mm=props['widths'][k]*step,
                                     smoothing_shift_mm=max(shifts)))
    return pd.DataFrame(rows, columns=['candidate_id','track','x_mm','depth_mm','kind','run_id',
                                     'prominence','prominence_snr','width_mm','smoothing_shift_mm'])


def fit_points(picks, cfg):
    x, z = picks.x_mm.to_numpy(), picks.depth_mm.to_numpy()
    xc = x - x.mean()
    slope = np.dot(xc, z-z.mean()) / np.dot(xc,xc)
    intercept = z.mean() - slope*x.mean()
    residual = z - (intercept+slope*x)
    sse = np.dot(residual,residual)
    sst = np.sum((z-z.mean())**2)
    # A conservative per-pick floor reflects grid spacing and smoothing sensitivity.
    floor = max(cfg.step_mm, picks.smoothing_shift_mm.max())
    se = np.sqrt(max(sse/(len(x)-2), floor**2)/np.dot(xc,xc))
    half = t.ppf(.975,len(x)-2)*se
    low, high = np.degrees(np.arctan([slope-half,slope+half]))
    angle = np.degrees(np.arctan(slope))
    loo = []
    for i in range(len(x)):
        keep = np.arange(len(x)) != i
        loo.append(np.degrees(np.arctan(np.polyfit(x[keep]-x[keep].mean(), z[keep]-z[keep].mean(),1)[0])))
    return dict(slope_mm_per_mm=slope, intercept_mm=intercept, depth_m=z.mean()/1000,
                angle_deg=angle, ci_low_deg=low, ci_high_deg=high,
                ci_half_deg=(high-low)/2, slope_se=se,
                r2=1-sse/sst if sst > 1e-12 else 1., rmse_mm=np.sqrt(sse/len(x)),
                max_residual_mm=np.max(np.abs(residual)), n_tracks=len(x),
                span_mm=np.ptp(x), depth_range_mm=np.ptp(z),
                loo_change_deg=max(np.abs(np.array(loo)-angle)),
                min_prominence_snr=picks.prominence_snr.min(),
                uncertainty_floor_mm=floor)


def match(candidates, tracks, cfg):
    total_span = tracks[-1]['x']-tracks[0]['x']
    required = max(cfg.min_tracks, int(np.ceil(cfg.min_coverage*len(tracks))))
    slope_limit = np.tan(np.radians(cfg.max_dip_deg))
    seen, proposals = set(), []
    for (kind, run_id), group in candidates.groupby(['kind','run_id']):
        items = list(group.itertuples())
        for left in items:
            for right in items:
                dx = right.x_mm-left.x_mm
                if dx <= 0 or dx < cfg.min_span_fraction*total_span:
                    continue
                slope = (right.depth_mm-left.depth_mm)/dx
                if abs(slope) > slope_limit:
                    continue
                predicted = left.depth_mm+slope*(group.x_mm-left.x_mm)
                nearby = group.loc[np.abs(group.depth_mm-predicted) <= cfg.match_mm]
                # Do not choose between competing extrema in a single track.
                if nearby.track.duplicated().any() or len(nearby) < required:
                    continue
                key = tuple(sorted(nearby.candidate_id))
                if key in seen:
                    continue
                seen.add(key)
                fit = fit_points(nearby, cfg)
                # Every pick is within the seed's physical 30-degree search cone (one grid tolerance).
                physical = np.all(np.abs(nearby.depth_mm-left.depth_mm) <=
                                  slope_limit*np.abs(nearby.x_mm-left.x_mm)+cfg.step_mm)
                reasons = []
                for failed, reason in [
                    (not physical,'search_cone'),
                    (abs(fit['angle_deg']) > cfg.max_dip_deg,'angle'),
                    (fit['rmse_mm'] > cfg.max_rmse_mm,'rmse'),
                    (fit['max_residual_mm'] > cfg.max_residual_mm,'residual'),
                    (fit['r2'] < cfg.min_r2 and fit['depth_range_mm'] > cfg.flat_range_mm,'r2'),
                    (cfg.max_ci_half_deg is not None and fit['ci_half_deg'] > cfg.max_ci_half_deg,'interval'),
                    (fit['loo_change_deg'] > cfg.max_loo_change_deg,'leave_one_out')]:
                    if failed: reasons.append(reason)
                fit.update(kind=kind, run_id=run_id, candidate_ids=key,
                           status='rejected' if reasons else 'eligible', reason=';'.join(reasons))
                proposals.append(fit)
    # Prefer coverage, then precision. A pick may contribute to only one accepted fit.
    proposals.sort(key=lambda f: (-f['n_tracks'], f['ci_half_deg'], f['rmse_mm']))
    used, accepted = set(), []
    for fit in proposals:
        if fit['status'] != 'eligible':
            continue
        overlap = used.intersection(fit['candidate_ids'])
        close = any(a['kind']==fit['kind'] and abs(a['depth_m']-fit['depth_m']) < .005 for a in accepted)
        if overlap or close:
            fit.update(status='duplicate', reason='shared_picks_or_nearby_feature')
            continue
        fit.update(status='accepted', feature_id=len(accepted)+1)
        used.update(fit['candidate_ids'])
        accepted.append(fit)
    return accepted, proposals


def plot_section(name, tracks, grid, values, safe, candidates, fits, output, cfg):
    fig, (ax, right) = plt.subplots(1,2,figsize=(10,9),sharey=True,
                                   gridspec_kw={'width_ratios':[2,1]},layout='constrained')
    xs = np.array([tr['x'] for tr in tracks]); spacing = np.min(np.diff(xs))
    # Same 10-mm median display and outer-track scaling as basic_plots.
    means = np.array([np.mean(tr['v'][(tr['z']>=grid[0])&(tr['z']<=grid[-1])]) for tr in tracks])
    scales = np.ones(len(tracks)); scales[[0,-1]] = np.mean(means[1:-1])/means[[0,-1]]
    shown = values*scales[:,None]
    valid = shown[:,safe].ravel()
    if not len(valid): valid = shown[np.isfinite(shown)]
    lo,hi = np.percentile(valid,[5,95]); norm = Normalize(lo,hi if hi>lo else lo+1e-15)
    for ti,tr in enumerate(tracks):
        color = median_filter(shown[ti],size=41,mode='nearest')
        ax.pcolormesh([tr['x']-spacing*.48,tr['x']+spacing*.48],edges(grid/1000),
                      color[:,None],cmap='Spectral',norm=norm,shading='flat',rasterized=True)
        e = edges(tr['z']/1000)
        for flag, style in STYLES.items():
            for a,b in runs(tr['flags'][flag]):
                ax.add_patch(Rectangle((tr['x']-spacing*.48,e[a]),spacing*.96,e[b]-e[a],
                                       linewidth=0,**style))
        other = tr['flags']['button'] & ~np.logical_or.reduce([tr['flags'][f] for f in STYLES])
        for a,b in runs(other):
            ax.add_patch(Rectangle((tr['x']-spacing*.48,e[a]),spacing*.96,e[b]-e[a],color='gray'))
    # Gray right-panel bands mark globally excluded depths.
    for a,b in runs(~safe):
        right.axhspan(grid[a]/1000,grid[b-1]/1000,color='.9',zorder=0)
    quality_max = max(2., max((f['ci_half_deg'] for f in fits), default=2.))
    quality_norm = Normalize(0,quality_max)
    cmap = plt.get_cmap('viridis_r')
    for fit in fits:
        picks = candidates[candidates.candidate_id.isin(fit['candidate_ids'])]
        marker = '^' if fit['kind']=='peak' else 'v'
        q = fit['ci_half_deg']; color = cmap(quality_norm(q))
        size = 35+90*(1-q/quality_max)
        xx = np.array([picks.x_mm.min(),picks.x_mm.max()])
        yy = (fit['intercept_mm']+fit['slope_mm_per_mm']*xx)/1000
        ax.plot(xx,yy,color='white',lw=3)
        ax.plot(xx,yy,color=color,lw=1.5)
        ax.scatter(picks.x_mm,picks.depth_mm/1000,s=32,marker=marker,c=[color],edgecolors='white',linewidths=.6,zorder=4)
        label = str(fit['feature_id'])
        ax.text(xx[-1]+1,yy[-1],label,fontsize=8,va='center')
        right.errorbar(fit['angle_deg'],fit['depth_m'],
                       xerr=[[fit['angle_deg']-fit['ci_low_deg']],[fit['ci_high_deg']-fit['angle_deg']]],
                       color=color,capsize=3,lw=1.3)
        right.scatter(fit['angle_deg'],fit['depth_m'],s=size,c=[color],marker=marker,edgecolors='black',linewidths=.5,zorder=4)
        right.annotate(label,(fit['angle_deg'],fit['depth_m']),xytext=(7,5),textcoords='offset points',fontsize=8)
    ax.set(xlabel='Lateral track position (mm)',ylabel='Depth (m)',title='AC ECM · individual picks and fits',
           xlim=(xs[0]-spacing/2,xs[-1]+spacing*1.1),ylim=(grid[-1]/1000,grid[0]/1000))
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    right.set(xlabel='Signed apparent dip (°)',title=f'{len(fits)} accepted features',xlim=(-32,32))
    right.axvline(0,color='.5',lw=.7,ls=':'); right.grid(axis='x',alpha=.2)
    if not fits:
        right.text(.5,.5,'No features passed\nstrict selection',ha='center',va='center',transform=right.transAxes)
    fig.colorbar(ScalarMappable(norm=norm,cmap='Spectral'),ax=ax,orientation='horizontal',fraction=.045,pad=.06,label='AC current (A); display smoothing 10 mm',extend='both')
    fig.colorbar(ScalarMappable(norm=quality_norm,cmap=cmap),ax=right,orientation='horizontal',fraction=.045,pad=.06,label='95% interval half-width (°)\nSmaller = stronger fit; larger point')
    handles = [Patch(**v,label=k.replace('_',' ')) for k,v in STYLES.items()]
    handles += [Patch(facecolor='gray',label='Other button flag'), Patch(facecolor='.9',label='Excluded depth (right)')]
    handles += [plt.Line2D([],[],marker='^',color='k',ls='',label='Peak'),plt.Line2D([],[],marker='v',color='k',ls='',label='Trough')]
    fig.legend(handles=handles,loc='outside lower center',ncol=3,fontsize=8)
    fig.suptitle(f'{name} · feature-based apparent dip\n95% intervals conditional on selected picks',fontsize=12)
    fig.savefig(output,dpi=180)
    plt.close(fig)


def run_all(cfg=None, sections=None):
    cfg = cfg or Settings()
    work = ROOT/'agentic_workspace/dip_angles/results'; work.mkdir(parents=True,exist_ok=True)
    figures = ROOT/'figures/dip_angles/dic1'; figures.mkdir(parents=True,exist_ok=True)
    all_candidates, all_fits, all_proposals, all_picks, summaries = [],[],[],[],[]
    for path in sorted((ROOT/'data/clean').glob('*-AC.csv')):
        name = path.stem; section = name.split('-')[1]
        if sections is not None and section not in sections: continue
        tracks, grid, values, safe = load(path,cfg)
        candidates = detect(tracks,grid,values,safe,cfg)
        fits, proposals = match(candidates,tracks,cfg)
        plot_section(name,tracks,grid,values,safe,candidates,fits,figures/f'{name}-dip.png',cfg)
        all_candidates.extend(candidates.assign(section=section).to_dict('records'))
        for fit in fits:
            all_fits.append(dict(section=section,**fit))
            picks = candidates[candidates.candidate_id.isin(fit['candidate_ids'])].copy()
            picks['residual_mm'] = picks.depth_mm-(fit['intercept_mm']+fit['slope_mm_per_mm']*picks.x_mm)
            all_picks.extend(picks.assign(section=section,feature_id=fit['feature_id']).to_dict('records'))
        all_proposals.extend(dict(section=section,**p) for p in proposals)
        summary = dict(section=section,tracks=len(tracks),clean_fraction=safe.mean(),
                       candidates=len(candidates),proposals=len(proposals),accepted=len(fits))
        summaries.append(summary); print(summary,flush=True)
    for name, records in [('candidates',all_candidates),('fits',all_fits),('proposals',all_proposals),('picks',all_picks),('section_summary',summaries)]:
        pd.DataFrame(records).to_csv(work/f'{name}.csv',index=False)
    (work/'settings.json').write_text(json.dumps(asdict(cfg),indent=2)+'\n')
    inputs = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted((ROOT/'data/clean').glob('*-AC.csv'))}
    (work/'input_sha256.json').write_text(json.dumps(inputs,indent=2)+'\n')
    return pd.DataFrame(summaries), pd.DataFrame(all_fits)

if __name__ == '__main__':
    run_all()

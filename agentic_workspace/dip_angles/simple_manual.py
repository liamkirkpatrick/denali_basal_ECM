"""Raw min/max picks in exact manual windows; unfiltered OLS fits."""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.stats import t
from scipy.ndimage import median_filter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Patch,Rectangle
from feature_ranges import ROOT,DEFAULT_CSV,read_ranges
sys.path.insert(0,str(ROOT/'scripts/core_scripts'))
from core_names import canonical_core

WORK=ROOT/'agentic_workspace/dip_angles/simple_results'
FIGURES=ROOT/'figures/dip_angles/simple_manual'
STYLES={
 'surface_imperfection':dict(facecolor='black',edgecolor='black'),
 'tephra_or_layer':dict(facecolor='none',edgecolor='purple',hatch='xx'),
 'rock_or_debris':dict(facecolor='none',edgecolor='black',hatch='oo'),
 'core_break':dict(facecolor='none',edgecolor='blue',hatch='///')}


def pick_extrema(frame, window):
    records=[]
    for x,track in frame.groupby('Y_dimension(mm)',sort=True):
        inside=track.loc[track['True_depth(m)'].between(window.depth_min_m,window.depth_max_m)]
        finite=inside.loc[np.isfinite(inside.meas)&np.isfinite(inside['True_depth(m)'])]
        finite=finite.sort_values('True_depth(m)',kind='stable')
        if finite.empty:
            records.append(dict(x_mm=x,depth_m=np.nan,current_A=np.nan,
                                status='no_finite_measurements_in_window'))
            continue
        value=finite.meas.max() if window.kind=='peak' else finite.meas.min()
        tied=finite.loc[finite.meas==value]
        row=tied.iloc[0]  # Deterministic tie: shallowest measured extremum.
        flags=[f for f in ['button','Button',*STYLES] if f in frame]
        records.append(dict(x_mm=x,depth_m=row['True_depth(m)'],current_A=row.meas,
                            raw_row_index=int(tied.index[0]),tied_extrema=len(tied),
                            flagged=bool(row[flags].ne(0).any()),status='picked'))
    return pd.DataFrame(records)


def fit_line(picks):
    good=picks[picks.status=='picked']; x=good.x_mm.to_numpy(); z=good.depth_m.to_numpy()*1000
    result=dict(n_tracks=len(x),depth_m=np.mean(z)/1000 if len(z) else np.nan,
                slope_mm_per_mm=np.nan,intercept_mm=np.nan,angle_deg=np.nan,
                ci_low_deg=np.nan,ci_high_deg=np.nan,ci_width_deg=np.nan,
                r2=np.nan,rmse_mm=np.nan,status='line_not_identifiable')
    if len(x)<2 or np.ptp(x)==0: return result
    xc=x-x.mean(); zc=z-z.mean(); sxx=np.dot(xc,xc)
    slope=np.dot(xc,zc)/sxx; intercept=z.mean()-slope*x.mean()
    residual=zc-slope*xc; sse=np.dot(residual,residual); sst=np.dot(zc,zc)
    result.update(slope_mm_per_mm=slope,intercept_mm=intercept,
                  angle_deg=np.degrees(np.arctan(slope)),rmse_mm=np.sqrt(sse/len(x)),
                  r2=1-sse/sst if sst>0 else np.nan,status='fitted')
    if len(x)>2:
        se=np.sqrt(sse/(len(x)-2)/sxx)
        half=t.ppf(.975,len(x)-2)*se
        low,high=np.degrees(np.arctan([slope-half,slope+half]))
        result.update(ci_low_deg=low,ci_high_deg=high,ci_width_deg=high-low)
    return result


def edges(z):
    return np.r_[z[0],(z[:-1]+z[1:])/2,z[-1]]


def plot_section(name,frame,picks,fits,path):
    fig,(ax,right)=plt.subplots(1,2,figsize=(11,9),sharey=True,
                                gridspec_kw={'width_ratios':[2,1]},layout='constrained')
    xs=np.sort(frame['Y_dimension(mm)'].unique()); spacing=min(np.diff(xs)) if len(xs)>1 else 5
    groups=[g.sort_values('True_depth(m)') for _,g in frame.groupby('Y_dimension(mm)',sort=True)]
    scales=np.ones(len(xs))
    if len(xs)>2:
        means=np.array([g.meas.mean() for g in groups])
        for i in [0,len(xs)-1]:
            if np.isfinite(means[i]) and means[i]!=0: scales[i]=np.mean(means[1:-1])/means[i]
    valid=[]
    for scale,g in zip(scales,groups):
        flags=[f for f in ['button','Button',*STYLES] if f in g]
        valid.extend(g.loc[~g[flags].ne(0).any(axis=1),'meas'].to_numpy()*scale)
    valid=np.asarray(valid);valid=valid[np.isfinite(valid)]
    if not len(valid): valid=frame.meas[np.isfinite(frame.meas)].to_numpy()
    lo,hi=np.percentile(valid,[5,95]); norm=Normalize(lo,hi if hi>lo else lo+1e-15,clip=True)
    for x,g,scale in zip(xs,groups,scales):
        z=g['True_depth(m)'].to_numpy();v=g.meas.to_numpy()*scale
        # Display only: 10-mm median smoothing; picks use raw values above.
        step=np.median(np.diff(z)); size=max(1,int(round(.01/step))|1)
        ax.pcolormesh([x-spacing*.48,x+spacing*.48],edges(z),median_filter(v,size=size)[:,None],
                      norm=norm,cmap='Spectral',shading='flat',rasterized=True)
        masks={f:g[f].ne(0).to_numpy() if f in g else np.zeros(len(g),bool) for f in STYLES}
        button='Button' if 'Button' in g else 'button'
        masks['other_button']=g[button].ne(0).to_numpy()&~np.logical_or.reduce(list(masks.values()))
        for flag,mask in masks.items():
            change=np.diff(np.r_[False,mask,False].astype(int)); e=edges(z)
            style=STYLES.get(flag,dict(facecolor='gray',edgecolor='gray'))
            for a,b in zip(np.flatnonzero(change==1),np.flatnonzero(change==-1)):
                ax.add_patch(Rectangle((x-spacing*.48,e[a]),spacing*.96,e[b]-e[a],linewidth=0,**style))
    for f in fits.itertuples():
        pp=picks[(picks.selection_id==f.selection_id)&(picks.status=='picked')]
        marker='^' if f.kind=='peak' else 'v'
        color='black'
        ax.scatter(pp.x_mm,pp.depth_m,s=26,marker=marker,c=[color],edgecolors='white',linewidths=.5,zorder=5)
        if f.status!='fitted': continue
        xx=np.array([pp.x_mm.min(),pp.x_mm.max()]); yy=(f.intercept_mm+f.slope_mm_per_mm*xx)/1000
        ax.plot(xx,yy,color='white',lw=3);ax.plot(xx,yy,color=color,lw=1.4)
        ax.text(xx[-1]+1,yy[-1],str(f.window_id),fontsize=7,va='center')
        size=45
        if np.isfinite(f.ci_width_deg):
            right.errorbar(f.angle_deg,f.depth_m,xerr=[[f.angle_deg-f.ci_low_deg],[f.ci_high_deg-f.angle_deg]],color=color,capsize=3)
        right.scatter(f.angle_deg,f.depth_m,s=size,marker=marker,c=[color],edgecolors='black',linewidths=.5,zorder=4)
        right.annotate(str(f.window_id),(f.angle_deg,f.depth_m),xytext=(5,5),textcoords='offset points',fontsize=7)
    ax.set(xlabel='Lateral track position (mm)',ylabel='Depth (m)',title='AC ECM',
           ylim=(frame['True_depth(m)'].max(),frame['True_depth(m)'].min()),xlim=(xs[0]-spacing/2,xs[-1]+spacing*1.3))
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    right.set(xlabel='Signed apparent dip (°)',title='Apparent dip estimates')
    right.axvline(0,color='.5',ls=':',lw=.7); right.margins(x=.15)
    if not len(fits):right.text(.5,.5,'No manual windows',transform=right.transAxes,ha='center')
    fig.colorbar(ScalarMappable(norm=norm,cmap='Spectral'),ax=ax,orientation='horizontal',fraction=.045,pad=.06,label='AC current (A) · 10 mm display smoothing',extend='both')
    handles=[Patch(**style,label=flag.replace('_',' ')) for flag,style in STYLES.items()]
    handles+=[Patch(facecolor='gray',label='Other button flag')]
    handles += [plt.Line2D([],[],marker=m,color='k',ls='',label=k) for k,m in [('Peak','^'),('Trough','v')]]
    fig.legend(handles=handles,loc='outside lower center',ncol=4,fontsize=8)
    fig.suptitle(name,fontsize=12)
    fig.savefig(path,dpi=180);plt.close(fig)


def plot_overall(fits, depth_bounds, core):
    """One continuous-depth panel for every estimate, without section grouping."""
    fig,ax=plt.subplots(figsize=(7,11),layout='constrained')
    for f in fits.sort_values('depth_m').itertuples():
        if f.status!='fitted': continue
        marker='^' if f.kind=='peak' else 'v'
        has_interval=np.isfinite(f.ci_width_deg)
        color='black'
        size=45
        if has_interval:
            ax.errorbar(f.angle_deg,f.depth_m,
                        xerr=[[f.angle_deg-f.ci_low_deg],[f.ci_high_deg-f.angle_deg]],
                        color=color,capsize=3,lw=1,alpha=.85,zorder=2)
        ax.scatter(f.angle_deg,f.depth_m,s=size,marker=marker,c=[color],
                   edgecolors='black',linewidths=.5,zorder=3)
    ax.axvline(0,color='.5',ls=':',lw=.8,zorder=0)
    ax.grid(axis='x',alpha=.15)
    pad=.01*(depth_bounds[1]-depth_bounds[0])
    ax.set(ylim=(depth_bounds[1]+pad,depth_bounds[0]-pad),
           xlabel='Signed apparent dip (°)',ylabel='Depth (m)',
           title=f'{core} · Apparent dip estimates')
    ax.ticklabel_format(axis='y',style='plain',useOffset=False)
    ax.margins(x=.08)
    ax.legend(handles=[plt.Line2D([],[],marker=m,color='k',ls='',label=k)
                       for k,m in [('Peak','^'),('Trough','v')]],loc='lower right',frameon=False)
    output=FIGURES/f'{core}-overall-dips.png'
    fig.savefig(output,dpi=200);plt.close(fig)
    return output


def run_simple():
    WORK.mkdir(exist_ok=True);FIGURES.mkdir(parents=True,exist_ok=True)
    windows=read_ranges().copy();windows['core']=windows.core.map(canonical_core)
    windows['window_id']=np.arange(1,len(windows)+1)
    if not windows.source_mode.eq('AC').all():raise ValueError('This analysis uses AC selections only.')
    all_picks=[];all_fits=[];sections=[];seen=set();depth_bounds={}
    hashes={str(DEFAULT_CSV.relative_to(ROOT)):hashlib.sha256(DEFAULT_CSV.read_bytes()).hexdigest()}
    for path in sorted((ROOT/'data/clean').glob('*-AC.csv')):
        source,section,face,_=path.stem.split('-');core=canonical_core(source)
        frame=pd.read_csv(path)
        bounds=(frame['True_depth(m)'].min(),frame['True_depth(m)'].max())
        previous=depth_bounds.get(core,bounds)
        depth_bounds[core]=(min(previous[0],bounds[0]),max(previous[1],bounds[1]))
        selected=windows[(windows.core==core)&(windows.section==section)&(windows.face==face)]
        local_picks=[];local_fits=[]
        for w in selected.itertuples():
            identity=dict(core=core,section=section,face=face,window_id=w.window_id,selection_id=w.selection_id,kind=w.kind)
            picks=pick_extrema(frame,w);fit=fit_line(picks)
            picks=picks.assign(**identity)
            if fit['status']=='fitted':
                picks['residual_mm']=picks.depth_m*1000-(fit['intercept_mm']+fit['slope_mm_per_mm']*picks.x_mm)
            local_picks.extend(picks.to_dict('records'));local_fits.append(dict(**identity,**fit))
            seen.add(w.selection_id)
        pp=pd.DataFrame(local_picks);ff=pd.DataFrame(local_fits)
        plot_section(f'{core}-{section}-{face}',frame,pp,ff,FIGURES/f'{core}-{section}-{face}-simple.png')
        all_picks.extend(local_picks);all_fits.extend(local_fits)
        sections.append(dict(core=core,section=section,windows=len(selected),tracks=frame['Y_dimension(mm)'].nunique()))
        hashes[str(path.relative_to(ROOT))]=hashlib.sha256(path.read_bytes()).hexdigest()
        print(f'{section}: {len(selected)} windows',flush=True)
    if seen!=set(windows.selection_id):raise ValueError('Some manual windows have no source file.')
    picks=pd.DataFrame(all_picks);fits=pd.DataFrame(all_fits)
    picks.to_csv(WORK/'picks.csv',index=False);fits.to_csv(WORK/'fits.csv',index=False)
    pd.DataFrame(sections).to_csv(WORK/'section_summary.csv',index=False)
    windows.to_csv(WORK/'input_windows.csv',index=False)
    (WORK/'input_sha256.json').write_text(json.dumps(hashes,indent=2)+'\n')
    for core,group in fits.groupby('core'):
        plot_overall(group,depth_bounds[core],core)
    return fits,picks

if __name__=='__main__':run_simple()

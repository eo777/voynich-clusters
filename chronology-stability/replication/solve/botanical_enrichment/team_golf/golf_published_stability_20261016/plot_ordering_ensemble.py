"""Render retained-observation ensemble frequencies without smoothing modes."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgb

HERE = Path(__file__).resolve().parent
TEAM = HERE.parent
GROUPS = ['1a','1b','1.5a','2a','3a','3b','3c','4a','4c','5a','6','']
LABELS = {'3a':'Herbal 3a','3b':'Herbal 3b','3c':'Stars 3c','4a':'Bio 4a','4c':'Stars 4c','':'Unassigned'}
INK = '#193244'
MUTED = '#647884'
SATURATION = .20


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_panel(cohort, sampler):
    datasets=[]; sources={}
    for penalty in [.01,.05]:
        filename=f'ensemble_frequency_{cohort}_{sampler}_l{int(penalty*1000):03d}.json'
        path=HERE/filename
        value=json.loads(path.read_text())
        assert value['cohort']==cohort and value['sampler']==sampler and value['penalty']==penalty
        assert value['scopes']['all_starts']['fits']==800
        assert value['scopes']['all_starts']['mask_draws']==100
        datasets.append(value); sources[filename]=digest(path)
    a,b=datasets
    assert a['meta']==b['meta'] and a['histogram_edges']==b['histogram_edges']
    rows=[]
    for meta,ra,rb in zip(a['meta'],a['scopes']['all_starts']['rows'],b['scopes']['all_starts']['rows']):
        assert meta['id']==ra['id']==rb['id']
        assert ra['retained_fits']==rb['retained_fits']==8*ra['retained_mask_draws']==8*rb['retained_mask_draws']
        counts=np.asarray(ra['histogram'])+np.asarray(rb['histogram'])
        denominator=ra['retained_fits']+rb['retained_fits']
        assert int(counts.sum())==denominator==16*ra['retained_mask_draws']
        assert denominator>0
        rows.append(dict(id=meta['id'],group=meta['group'],bifolium=meta['bifolium'],
            counts=counts.tolist(),denominator=denominator,mask_draws_retaining_row=ra['retained_mask_draws'],
            frequency=(counts/denominator).tolist()))
    return dict(cohort=cohort,sampler=sampler,rows=rows,source_sha256=sources,
        histogram_edges=a['histogram_edges'],fits=1600,independent_masks=100,
        group_counts={g:sum(r['group']==g for r in rows) for g in GROUPS})


def rgba(frequency, color):
    """One common declared frequency scale, not per-row peak normalization."""
    rgb=np.asarray(to_rgb(color))
    alpha=np.sqrt(np.clip(np.asarray(frequency)/SATURATION,0,1))
    out=np.ones((len(alpha),4));out[:,:3]=rgb;out[:,3]=alpha
    return out[None,:,:]


def draw(mode='native'):
    source_paths=[HERE/'ensemble_index.json',HERE/'ENSEMBLE_DESIGN.md',Path(__file__)]
    index=json.loads(source_paths[0].read_text())
    assert index['status']=='complete'
    panels=[load_panel('paragraph207','point' if mode=='native' else 'sheet'),load_panel('pool50','point')]
    name='ordering_ensemble' if mode=='native' else 'ordering_ensemble_sheet_control'
    palette_path=TEAM/'nonlinear_sequence_20261010/plot_data.json'
    palette=json.loads(palette_path.read_text())['palette'];palette['']='#7d8990'
    source_paths.append(palette_path)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,
        'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':'#cad5db',
        'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    counts=panels[0]['group_counts']
    heights={g:max(.58,counts[g]/22) for g in GROUPS}
    bounds={};cursor=0.
    for g in GROUPS:
        bounds[g]=(cursor,cursor+heights[g]);cursor+=heights[g]+.11
    fig,axes=plt.subplots(1,2,figsize=(12,9.2))
    fig.subplots_adjust(left=.12,right=.96,bottom=.17,top=.84,wspace=.23)
    strips=0
    for col,(ax,panel) in enumerate(zip(axes,panels)):
        for g in GROUPS:
            lo,hi=bounds[g];color=palette.get(g,'#7d8990')
            rows=[r for r in panel['rows'] if r['group']==g]
            ax.hlines(hi+.05,-.015,1.015,color='#e2e8eb',lw=.45)
            step=(hi-lo)/max(len(rows),1)
            for i,row in enumerate(rows):
                center=lo+(i+.5)*step
                # The fixed 100 bins are shown directly; no KDE or filled interval.
                ax.imshow(rgba(row['frequency'],color),extent=[0,1,center+.37*step,center-.37*step],
                          aspect='auto',interpolation='nearest',zorder=3)
                strips+=1
            if col==0:
                ax.text(-.035,(lo+hi)/2,LABELS.get(g,g),ha='right',va='center',color=color,fontsize=10.3)
            ax.text(1.025,(lo+hi)/2,str(len(rows)) if rows else '—',ha='left',va='center',color=MUTED,fontsize=8.5)
        title='Individual pages · 207' if col==0 else 'Pooled bifolia · 50'
        ax.set_title(title,loc='left',fontsize=13,pad=13,fontweight='bold')
        ax.text(1.025,1.013,'n',transform=ax.transAxes,fontsize=8.5,color=MUTED,ha='left')
        ax.set_xlim(-.012,1.012);ax.set_ylim(cursor-.09,-.1);ax.set_yticks([])
        ax.spines['left'].set_visible(False);ax.set_xticks(np.arange(0,1.01,.2))
        ax.grid(axis='x',color='#dce4e8',lw=.6,alpha=.65);ax.set_axisbelow(True)
        ax.set_xlabel('Position along each fitted curve',labelpad=7)
        ax.text(0,-.055,'earlier',transform=ax.transAxes,ha='left',va='top',fontsize=8.5,color=MUTED)
        ax.text(1,-.055,'later',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color=MUTED)
    title='Golf ordering across resampled curves' if mode=='native' else 'Golf ordering · physical-sheet sampling control'
    fig.suptitle(title,x=.12,y=.976,ha='left',fontsize=18,fontweight='bold')
    fig.text(.12,.933,'100 random 90% subsamples · 8 starts per penalty · λ = 0.01 / 0.05',fontsize=10.8,color=MUTED)
    fig.text(.12,.900,'One strip per observation; darker bins contain more retained fits.',fontsize=10,color=MUTED)
    # The neutral legend displays the same square-root opacity as every Golf hue.
    lax=fig.add_axes([.39,.087,.30,.018])
    freq=np.linspace(0,SATURATION,401)
    lax.imshow(rgba(freq,'#263f50'),extent=[0,SATURATION,0,1],aspect='auto',interpolation='nearest')
    lax.set_yticks([]);lax.set_xticks([0,.01,.05,.1,.2],labels=['0','1','5','10','20+'])
    lax.tick_params(axis='x',length=2,labelsize=8)
    for spine in lax.spines.values():spine.set_visible(False)
    fig.text(.54,.048,'Retained fits per 0.01-wide position bin (%)',ha='center',fontsize=9,color=MUTED)
    fig.text(.12,.016,'Sensitivity frequencies, not historical probabilities. Rows are grouped by Golf label, not ordered by the fitted chronology.',fontsize=8.6,color=MUTED)
    for ext in ['png','svg']:
        fig.savefig(HERE/f'{name}.{ext}',dpi=160,facecolor='white')
    plt.close(fig)
    assert strips==257
    sources={str(p.relative_to(TEAM)):digest(p) for p in source_paths}
    for panel in panels:sources.update(panel['source_sha256'])
    data=dict(status='complete',source_sha256=sources,panels=panels,group_order=GROUPS,
        primary_scope='retained observations, all eight starts per penalty; equal penalty/start weights',
        sampler='native points in both cohorts' if mode=='native' else 'whole physical sheets in both cohorts',
        color_scale=dict(frequency_cap=SATURATION,opacity='sqrt(min(frequency/.20,1)); common to every row/hue'),
        caveats=['100 masks are shared across starts/penalties; 1600 fits are not1600 independent resamples.',
            'Start permutations are selected separately by cohort and penalty; differences are overall model sensitivity, not isolated lambda effects.',
            'Fixed full-data positive-MDS geometry and initial anchors; distances and MDS are not re-estimated.',
            'Native page and pooled Golf labels differ for mixed-label bifolia.'])
    (HERE/f'{name}_data.json').write_text(json.dumps(data,separators=(',',':'))+'\n')
    provenance=dict(status='complete',source_sha256=sources,observations=strips,
        figure=f'{name}.png',histogram_bins=100,frequency_scope=data['primary_scope'],
        scale=data['color_scale'],sampler=data['sampler'],new_fits=0,
        outputs={f'{name}.{ext}':dict(bytes=(HERE/f'{name}.{ext}').stat().st_size,
            sha256=digest(HERE/f'{name}.{ext}')) for ext in ['png','svg']})
    (HERE/f'{name}_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(dict(figure=name,strips=strips,panel_sizes=[len(p['rows']) for p in panels],
        total_retained_positions=[sum(r['denominator'] for r in p['rows']) for p in panels]),indent=2))


if __name__=='__main__':
    mode=sys.argv[1] if len(sys.argv)>1 else 'native'
    assert mode in ['native','sheet_control']
    draw(mode)

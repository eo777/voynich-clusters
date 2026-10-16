"""RF-only, equal-weight combined, and side-by-side retained-fit strips."""
from pathlib import Path
import hashlib
import importlib.util
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
TEAM=HERE.parent
OLD=TEAM/'golf_published_stability_20261016'
s=importlib.util.spec_from_file_location('original_strip_plot',OLD/'plot_ordering_ensemble.py')
P=importlib.util.module_from_spec(s);s.loader.exec_module(P)

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def panel(root,cohort,sampler):
    P.HERE=root
    x=P.load_panel(cohort,sampler)
    x['source_sha256']={root.name+'/'+k:v for k,v in x['source_sha256'].items()}
    return x

def combine(a,b):
    assert a['histogram_edges']==b['histogram_edges']
    rows=[]
    for r,s in zip(a['rows'],b['rows']):
        for k in ['id','group','bifolium','denominator','mask_draws_retaining_row']:assert r[k]==s[k],(k,r['id'])
        counts=np.array(r['counts'])+np.array(s['counts']);denom=r['denominator']+s['denominator']
        rows.append(dict(r,counts=counts.tolist(),denominator=denom,frequency=(counts/denom).tolist()))
    return dict(a,rows=rows,fits=a['fits']+b['fits'],representation_weights={'ZL_Zattera':.5,'RF1b_basic_EVA':.5},
        source_sha256={**a['source_sha256'],**b['source_sha256']})

def draw(view,sheet=False):
    sam='sheet' if sheet else 'point'
    z=[panel(OLD,'paragraph207',sam),panel(OLD,'pool50','point')]
    r=[panel(HERE,'paragraph207',sam),panel(HERE,'pool50','point')]
    if view=='rf1b':
        panels=r;titles=['Individual pages · 207','Pooled bifolia · 50'];title='RF1b-er · basic EVA ordering sensitivity'
    elif view=='combined':
        panels=[combine(a,b) for a,b in zip(z,r)];titles=['Individual pages · 207','Pooled bifolia · 50'];title='ZL + RF1b-er · combined ordering sensitivity'
    else:
        panels=[z[0],r[0],z[1],r[1]];titles=['Pages · ZL / Zattera','Pages · RF / basic EVA','Bifolia · ZL / Zattera','Bifolia · RF / basic EVA'];title='Ordering sensitivity · ZL and RF1b-er'
    palette=json.loads((TEAM/'nonlinear_sequence_20261010/plot_data.json').read_text())['palette'];palette['']='#7d8990'
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':P.INK,'axes.labelcolor':P.INK,
        'xtick.color':P.MUTED,'ytick.color':P.MUTED,'axes.edgecolor':'#cad5db','axes.spines.top':False,
        'axes.spines.right':False,'svg.fonttype':'none'})
    bounds={};cursor=0
    for g in P.GROUPS:
        h=max(.58,z[0]['group_counts'][g]/22);bounds[g]=(cursor,cursor+h);cursor+=h+.11
    wide=len(panels)==4
    fig,axes=plt.subplots(1,len(panels),figsize=(17.5 if wide else 12,9.2))
    left=.085 if wide else .12
    fig.subplots_adjust(left=left,right=.97,bottom=.17,top=.84,wspace=.22)
    strips=0
    for col,(ax,d,title_col) in enumerate(zip(axes,panels,titles)):
        for g in P.GROUPS:
            lo,hi=bounds[g];color=palette.get(g,'#7d8990');rows=[r for r in d['rows'] if r['group']==g]
            ax.hlines(hi+.05,-.015,1.015,color='#e2e8eb',lw=.45);step=(hi-lo)/max(len(rows),1)
            for i,row in enumerate(rows):
                y=lo+(i+.5)*step
                ax.imshow(P.rgba(row['frequency'],color),extent=[0,1,y+.37*step,y-.37*step],aspect='auto',interpolation='nearest',zorder=3)
                strips+=1
            if col==0:ax.text(-.035,(lo+hi)/2,P.LABELS.get(g,g),ha='right',va='center',color=color,fontsize=10.3)
            ax.text(1.025,(lo+hi)/2,str(len(rows)) if rows else '—',ha='left',va='center',color=P.MUTED,fontsize=8)
        ax.set_title(title_col,loc='left',fontsize=11.5 if wide else 13,pad=13,fontweight='bold')
        ax.set_xlim(-.012,1.012);ax.set_ylim(cursor-.09,-.1);ax.set_yticks([]);ax.spines['left'].set_visible(False)
        ax.set_xticks(np.arange(0,1.01,.2));ax.grid(axis='x',color='#dce4e8',lw=.6,alpha=.65);ax.set_axisbelow(True)
        ax.set_xlabel('Position along fitted curve',labelpad=7)
        ax.text(0,-.055,'earlier',transform=ax.transAxes,ha='left',va='top',fontsize=8.5,color=P.MUTED)
        ax.text(1,-.055,'later',transform=ax.transAxes,ha='right',va='top',fontsize=8.5,color=P.MUTED)
    fig.suptitle(title,x=left,y=.976,ha='left',fontsize=18,fontweight='bold')
    note='Whole-bifolium sampling' if sheet else 'Page / bifolium sampling'
    fig.text(left,.933,f'{note} · 100 subsamples · 8 starts · λ = 0.01 / 0.05',fontsize=10.8,color=P.MUTED)
    desc='50% ZL / Zattera + 50% RF / basic EVA · matched rows and sampling masks.' if view=='combined' else 'One strip per observation; darker bins contain more retained fits.'
    fig.text(left,.900,desc,fontsize=10,color=P.MUTED)
    lax=fig.add_axes([.39,.087,.30,.018]);freq=np.linspace(0,P.SATURATION,401)
    lax.imshow(P.rgba(freq,'#263f50'),extent=[0,P.SATURATION,0,1],aspect='auto',interpolation='nearest')
    lax.set_yticks([]);lax.set_xticks([0,.01,.05,.1,.2],labels=['0','1','5','10','20+']);lax.tick_params(axis='x',length=2,labelsize=8)
    for spine in lax.spines.values():spine.set_visible(False)
    fig.text(.54,.048,'Retained fits per 0.01-wide position bin (%)',ha='center',fontsize=9,color=P.MUTED)
    footer='Sensitivity frequencies, not historical probabilities.'
    if view!='rf1b':footer+=' Both transcription and glyph grouping differ.'
    fig.text(left,.016,footer,fontsize=8.6,color=P.MUTED)
    name='ordering_'+view+('_sheet_control' if sheet else '')
    for ext in ['png','svg']:fig.savefig(HERE/f'{name}.{ext}',dpi=160,facecolor='white')
    plt.close(fig)
    sources={k:v for d in panels for k,v in d['source_sha256'].items()}
    sources['plot_comparison.py']=sha(Path(__file__));sources['original_plot_code']=sha(OLD/'plot_ordering_ensemble.py')
    output=dict(status='complete',figure=name,panels=panels,strips=strips,source_sha256=sources,
        scope='retained fits, all starts and penalties',representation_weight='one half each' if view=='combined' else 'separate',
        equal_weight_is_not_a_posterior=True,shared_masks_not_independent_replications=True,
        frequency_cap=P.SATURATION,outputs={ext:sha(HERE/f'{name}.{ext}') for ext in ['png','svg']})
    (HERE/f'{name}_data.json').write_text(json.dumps(output,separators=(',',':'),allow_nan=False),encoding='utf-8')
    print(json.dumps({'figure':name,'strips':strips}),flush=True)

if __name__=='__main__':
    for sheet in [False,True]:
        for view in ['rf1b','combined','comparison']:draw(view,sheet)

"""Original-Golf scEGRET continuity, preserving every original bin number."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Patch
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent
SOURCE=HERE/'scegret_summary.json'
INK='#193244'; MUTED='#647884'
GROUPS=['1a','1b','1.5a','2a','3a','3b','3c','4a','4c','5a','6','']
CHANNEL={'word':'#226eaa','line':'#d56b22'}


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    data=json.loads(SOURCE.read_text());assert data['status']=='complete'
    palette_path=HERE.parent/'nonlinear_sequence_20261010/plot_data.json'
    palette=json.loads(palette_path.read_text())['palette'];palette['']='#7d8990'
    # Middle predeclared resolutions; never chosen after inspecting peaks.
    selected=[]
    for cohort,steps in [('paragraph207',12),('pool50',8)]:
        pair={channel:next(r for r in data['original_golf'] if r['cohort']==cohort
            and r['steps']==steps and r['matrix']==channel+'_per_million') for channel in ['word','line']}
        assert all(r['status']=='complete' for r in pair.values())
        assert pair['word']['occupancy']==pair['line']['occupancy']
        selected.append((cohort,steps,pair))
    maximum=max(v for _,_,pair in selected for v in pair['word']['full_pair_changes'])
    norm=Normalize(vmin=0,vmax=maximum)
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'text.color':INK,
        'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,'axes.edgecolor':'#c4d0d7',
        'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig=plt.figure(figsize=(12,8.5))
    gs=fig.add_gridspec(3,2,left=.08,right=.89,bottom=.23,top=.85,
        height_ratios=[3.1,.95,1.55],hspace=.26,wspace=.22)
    summaries=[]
    for col,(cohort,steps,pair) in enumerate(selected):
        top=fig.add_subplot(gs[0,col]);occupancy_ax=fig.add_subplot(gs[1,col],sharex=top)
        adjacent_ax=fig.add_subplot(gs[2,col],sharex=top)
        primary=pair['word'];retained=set(primary['retained_bins'])
        values={tuple(p):v for p,v in zip(primary['full_pairs'],primary['full_pair_changes'])}
        patches=[];colors=[];missing=0
        for a in range(1,steps+1):
            for b in range(a+1,steps+1):
                cx,cy=(a+b)/2,(b-a)/2
                vertices=[(cx-.5,cy),(cx,cy+.5),(cx+.5,cy),(cx,cy-.5)]
                if (a,b) in values:
                    patches.append(Polygon(vertices));colors.append(values[(a,b)])
                else:
                    top.add_patch(Polygon(vertices,facecolor='#eef1f3',edgecolor='#acbac3',
                        hatch='///',linewidth=.3));missing+=1
        collection=PatchCollection(patches,cmap='YlOrRd',norm=norm,edgecolor='white',linewidth=.45)
        collection.set_array(np.asarray(colors));top.add_collection(collection)
        top.set_xlim(.45,steps+.55);top.set_ylim(-.1,steps/2+.45)
        top.set_aspect('equal',adjustable='box');top.set_yticks([])
        top.set_xticks(range(1,steps+1));top.tick_params(axis='x',length=0,labelbottom=False)
        top.spines[['left','bottom']].set_visible(False)
        title='Individual pages · 207' if cohort=='paragraph207' else 'Pooled bifolia · 50'
        top.set_title(title,loc='left',fontsize=13,fontweight='bold',pad=10)
        bottom=np.zeros(steps)
        for group in GROUPS:
            counts=np.array([c['golf_counts'].get(group,0) for c in primary['bin_composition']])
            occupancy_ax.bar(range(1,steps+1),counts,bottom=bottom,width=.8,
                color=palette.get(group,'#7d8990'),edgecolor='white',linewidth=.3)
            bottom+=counts
        assert bottom.tolist()==primary['occupancy']
        for i,count in enumerate(bottom,1):
            occupancy_ax.text(i,count+.02*max(bottom),str(int(count)),ha='center',va='bottom',fontsize=7.8,color=INK)
            if i not in retained:
                occupancy_ax.bar(i,max(count,.7),width=.8,facecolor='none',edgecolor='#72848e',hatch='///',linewidth=.6)
        occupancy_ax.set_ylim(0,max(bottom)*1.32);occupancy_ax.set_yticks([])
        occupancy_ax.tick_params(axis='x',length=0,labelbottom=False)
        occupancy_ax.spines[['left','bottom']].set_visible(False)
        occupancy_ax.set_ylabel('Bin count',fontsize=9,labelpad=7)
        for channel in ['word','line']:
            for entry in pair[channel]['adjacent']:
                a,b=entry['bins'];value=entry['value'];gap=entry['spans_excluded_bins']
                adjacent_ax.plot([a,b],[value,value],color=CHANNEL[channel],lw=1.25,
                    linestyle='--' if gap else '-',alpha=.8)
                adjacent_ax.scatter([(a+b)/2],[value],s=17,marker='D' if gap else 'o',
                    facecolors='white' if gap else CHANNEL[channel],edgecolors=CHANNEL[channel],linewidth=.9,zorder=4)
        for i in range(1,steps+1):
            if i not in retained:adjacent_ax.axvspan(i-.45,i+.45,color='#d7dfe4',alpha=.35,lw=0)
        adjacent_ax.set_ylim(bottom=0)
        adjacent_ax.set_xticks(range(1,steps+1),labels=[str(i)+('*' if i in primary['H_focal_bins'] else '')
            for i in range(1,steps+1)])
        adjacent_ax.tick_params(axis='both',labelsize=8,length=3)
        adjacent_ax.grid(axis='y',color='#dce4e8',lw=.5)
        adjacent_ax.set_ylabel('Relative change',fontsize=9)
        adjacent_ax.set_xlabel('Original pseudotime bin',fontsize=10)
        summaries.append(dict(cohort=cohort,steps=steps,occupancy=primary['occupancy'],
            retained_bins=primary['retained_bins'],missing_pair_tiles=missing,
            H_focal_bins=primary['H_focal_bins'],H_bypass_comparisons=primary['H_bypass_comparisons'],
            focal_pair_bin_resolution=primary['focal_pair_bin_resolution'],
            pyramid_source='word_per_million',feature_count=primary['feature_count'],
            largest_steps={channel:pair[channel]['maximum_step'] for channel in ['word','line']},
            source_summary={channel:pair[channel]['source_retained_adjacent_summary'] for channel in ['word','line']}))
    fig.suptitle('Original Golf ordering · scEGRET continuity',x=.08,y=.976,ha='left',fontsize=18,fontweight='bold')
    fig.text(.08,.932,'All-bin-pair pyramid above; adjacent retained-bin changes below',fontsize=10.8,color=MUTED)
    fig.text(.08,.902,'Per-million transition profiles · middle tested resolution: 12 / 8 bins',fontsize=10,color=MUTED)
    cax=fig.add_axes([.92,.505,.018,.27])
    cb=fig.colorbar(collection,cax=cax);cb.ax.tick_params(labelsize=8,length=2)
    cb.set_label('Relative word-transition change',fontsize=9,labelpad=7)
    handles=[Line2D([],[],color=CHANNEL[c],marker='o',lw=1.3,ms=4,label=c.title()+' transitions') for c in ['word','line']]
    handles+=[Patch(facecolor='#eef1f3',edgecolor='#acbac3',hatch='///',label='Bin with ≤3 observations'),
              Line2D([],[],color='#647884',ls='--',marker='D',markerfacecolor='white',lw=1.2,ms=4,label='Comparison spans excluded bins')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.142),ncol=2,frameon=False,
        fontsize=9,columnspacing=2.4,handlelength=2.2)
    golf=[Patch(facecolor=palette.get(g,'#7d8990'),label=g or 'unassigned') for g in GROUPS]
    fig.legend(handles=golf,loc='lower center',bbox_to_anchor=(.5,.052),ncol=12,frameon=False,
        fontsize=9,columnspacing=1.05,handlelength=1.2,title='Golf composition of each original bin',title_fontsize=9)
    fig.text(.5,.023,'* Bin containing f58r/v (pages) or H:58–65 (pooled). Feature contrasts do not establish missing manuscript sections.',ha='center',fontsize=8.7,color=MUTED)
    for extension in ['png','svg']:fig.savefig(HERE/f'scegret_original_golf.{extension}',dpi=160,facecolor='white')
    plt.close(fig)
    provenance=dict(status='complete',source_sha256=sha(SOURCE),rendering_source_sha256=sha(__file__),
        palette_sha256=sha(palette_path),panels=summaries,
        selection='Original Golf spline only; middle predeclared resolution12pages/8pools; wordperMillionpyramid, word/lineperMillionadjacent scores.',
        bins='Every original bin ID remains on the axis; omitted bins and pair tiles marked, no factor-code compression.',
        new_fits=0,new_scegret_runs=0,outputs={f'scegret_original_golf.{ext}':dict(
            bytes=(HERE/f'scegret_original_golf.{ext}').stat().st_size,sha256=sha(HERE/f'scegret_original_golf.{ext}')) for ext in ['png','svg']})
    (HERE/'scegret_original_golf_provenance.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(json.dumps(summaries,indent=2))


if __name__=='__main__':main()

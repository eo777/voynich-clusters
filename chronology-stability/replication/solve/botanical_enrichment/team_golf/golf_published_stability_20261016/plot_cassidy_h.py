"""Render corrected Cassidy continuation; no fitting or statistical resampling."""
from pathlib import Path
import hashlib
import json
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle

HERE = Path(__file__).resolve().parent
SOURCE = HERE/'cassidy_h_results.json'
INK = '#193244'


def edges(values):
    a = np.asarray(values)
    assert (np.diff(a)>0).all()
    return np.r_[a[0]-(a[1]-a[0])/2, (a[1:]+a[:-1])/2, a[-1]+(a[-1]-a[-2])/2]


def sign(value):
    return 0 if value is None or abs(value)<1e-9 else (1 if value>0 else -1)


def main():
    data = json.loads(SOURCE.read_text())
    assert data['status']=='complete' and data['design']['penalty']==.01
    xy = np.asarray(data['reference_coords3']); hi = data['h_index']
    gaps = [c['combined_herbal_preference']['gap_percent'] for s in data['settings'] for c in s['cells']]
    limit = math.ceil(max(abs(v) for v in gaps if v is not None)*10)/10
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':11, 'axes.titlesize':12,
        'axes.labelsize':11, 'axes.edgecolor':'#758793', 'xtick.color':'#49616f',
        'ytick.color':'#49616f', 'text.color':INK, 'axes.labelcolor':INK, 'svg.fonttype':'none'})
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 5.9), sharex=True)
    fig.subplots_adjust(left=.085, right=.965, top=.81, bottom=.34, wspace=.18)
    summaries = []
    any_disagreement = False
    for ax, row in zip(axes, data['settings']):
        nx, ny = len(row['xs']), len(row['ys'])
        z = np.array([c['combined_herbal_preference']['gap_percent'] for c in row['cells']], dtype=float).reshape(ny,nx)
        disagreement = np.array([len(set(sign(p['herbal_preference']['gap_percent']) for p in c['paths']))>1
                                 for c in row['cells']]).reshape(ny,nx)
        any_disagreement |= bool(disagreement.any())
        im = ax.pcolormesh(edges(row['xs']), edges(row['ys']), z, cmap='RdBu_r', norm=norm, shading='flat')
        ax.contour(row['xs'], row['ys'], z, levels=[0], colors=INK, linewidths=.9)
        if disagreement.any():
            xe, ye = edges(row['xs']), edges(row['ys'])
            for iy, ix in np.argwhere(disagreement):
                ax.add_patch(Rectangle((xe[ix],ye[iy]),xe[ix+1]-xe[ix],ye[iy+1]-ye[iy],
                    facecolor='none',edgecolor='#526573',hatch='////',linewidth=0,alpha=.6))
        plane = row['plane']
        for i, meta in enumerate(data['meta']):
            if i != hi:
                ax.scatter(*xy[i, plane], s=23, c=data['palette'].get(meta['group'],'#8b949b'),
                           edgecolors='white', linewidths=.7, zorder=3)
        ax.scatter(*xy[hi,plane], s=150, marker='*', facecolor='white', edgecolor=INK, linewidth=1.2, zorder=5)
        ax.set_xlim(edges(row['xs'])[[0,-1]]); ax.set_ylim(edges(row['ys'])[[0,-1]])
        ax.set_xticks([-.2,0,.2]); ax.set_yticks([-.1,0,.1,.2] if plane[1]==1 else [-.1,0,.1])
        ax.set_xlabel('MDS 1'); ax.set_ylabel('MDS '+str(plane[1]+1))
        ax.set_aspect('equal', adjustable='box'); ax.tick_params(length=3)
        old = np.array([c['original_herbal_preference']['gap_percent'] for c in row['cells']])
        native = row['cells'][row['baseline_cell_index']]
        summaries.append(dict(plane=plane, three_a_first_cells=int(np.sum(z>1e-9)),
            three_b_first_cells=int(np.sum(z < -1e-9)), tied_cells=int(np.sum(abs(z)<=1e-9)),
            unavailable_cells=int(np.sum(~np.isfinite(z))), path_preference_disagreement_cells=int(disagreement.sum()),
            original_map_preference_changes=int(np.sum(np.sign(old)!=np.sign(z.ravel()))),
            maximum_gap_change_from_original_percent=float(np.max(abs(old-z.ravel()))),
            native_H_gap_percent=native['combined_herbal_preference']['gap_percent'],
            any_seed_path_order_disagreement_cells=sum(c['path_order_changes']>0 for c in row['cells']),
            any_seed_return_order_change_cells=sum(c['return_order_changes']>0 for c in row['cells']),
            maximum_path_node_rms=max(c['maximum_path_node_rms'] for c in row['cells']),
            maximum_return_node_rms=max(c['maximum_return_node_rms'] for c in row['cells']),
            terminal_4c_first_cells=sum(c['combined_terminal_preference']['gap_percent'] is not None
                and c['combined_terminal_preference']['gap_percent']<0 for c in row['cells'])))
    fig.suptitle('H:58–65 and the Herbal order', y=.972, fontsize=21, fontweight='bold')
    fig.text(.5,.909,'Cassidy derivative prediction + correction · λ = 0.01',ha='center',fontsize=11,color='#526b7b')
    fig.text(.5,.867,'Fixed-coordinate slices · final 1a-before-1b routes',ha='center',fontsize=10,color='#526b7b')
    cax = fig.add_axes([.28,.224,.44,.027])
    cb = fig.colorbar(im,cax=cax,orientation='horizontal')
    cb.ax.tick_params(labelsize=9,length=2); cb.outline.set_edgecolor('#b5c4cd')
    fig.text(.26,.238,'3b before 3a',ha='right',va='center',fontsize=11,color='#2166ac',fontweight='bold')
    fig.text(.74,.238,'3a before 3b',ha='left',va='center',fontsize=11,color='#b2182b',fontweight='bold')
    fig.text(.5,.163,'Relative fit advantage (%) · line: approximate equal score',ha='center',fontsize=9.5,color='#526b7b')
    groups = [g for g in ['1a','1b','1.5a','2a','3a','3b','3c','4a','4c','5a','6','']
              if any(m['group']==g and i!=hi for i,m in enumerate(data['meta']))]
    handles = [Line2D([],[],marker='o',linestyle='',markersize=5.8,
        markerfacecolor=data['palette'].get(g,'#8b949b'),markeredgecolor='white',label=g or 'unassigned') for g in groups]
    handles.append(Line2D([],[],marker='*',linestyle='',markersize=11,
        markerfacecolor='white',markeredgecolor=INK,label='observed H'))
    if any_disagreement:
        handles.append(Patch(facecolor='none',edgecolor='#536776',hatch='////',label='path-dependent preference'))
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.005),ncol=7,
               frameon=False,columnspacing=1.3,handletextpad=.35,fontsize=9.5,title='Golf groups and markers',title_fontsize=9.5)
    for extension in ['png','svg']:
        fig.savefig(HERE/f'cassidy_H_regions.{extension}',dpi=160,facecolor='white')
    plt.close(fig)
    result = dict(status='complete',source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        rendering_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        panels=summaries, color_limit_percent=limit,
        grid='Actual midpoint cell boundaries; displayed zero contour is interpolation of corrected grid scores.',
        score='100*(best3b-first objective-best3a-first objective)/min(best directional objectives); minimum over three paths and ten seed families, final1a-first outcomes only.',
        outputs={name:dict(bytes=(HERE/name).stat().st_size,sha256=hashlib.sha256((HERE/name).read_bytes()).hexdigest())
                 for name in ['cassidy_H_regions.png','cassidy_H_regions.svg']})
    (HERE/'cassidy_H_regions_provenance.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result['panels'],indent=2))


if __name__=='__main__':
    main()

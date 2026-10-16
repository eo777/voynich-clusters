"""Rebuild matched transition geometry; verify the implementation against ZL."""
from pathlib import Path
from collections import Counter
from functools import lru_cache
import csv
import hashlib
import importlib.util
import json
import sys
import numpy as np
from scipy.special import rel_entr
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
TEAM=HERE.parent
ROOT=TEAM.parents[2]
OLD=TEAM/'golf_published_stability_20261016'
DATA=ROOT/'apps/clustermatrix/research-data'
sys.path.insert(0,str(ROOT))
from voynich_toolkit.source_lattice import SourceLattice
from voynich_toolkit.glyph_units import units

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,obj): (HERE/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def module(name,path):
    s=importlib.util.spec_from_file_location(name,path); m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

@lru_cache(None)
def glyphs(word,mode):
    return tuple('other' if u.rare else u.text for u in units(word,mode))

def count_source(path,meta,mode,verify=False):
    doc=SourceLattice(path,verify=verify)
    assert doc.roundtrip()==doc.raw
    selected={m['id'] for m in meta}
    counts={k:{p:Counter() for p in selected} for k in ['word','line']}
    words=Counter(); loci=Counter(); excluded=Counter();all_pages=set();chars=set()
    for line in doc.stream('all_spaces','preferred_first'):
        if line.locus.category!='para': continue
        p=line.locus.page; all_pages.add(p)
        if p not in selected: continue
        ts=line.tokens;good=[not t.has_unreadable and '?' not in t.text for t in ts]
        loci[p]+=1;excluded[p]+=len(ts)-sum(good);words[p]+=sum(good)
        for t,ok in zip(ts,good):
            if not ok: continue
            chars.update(t.text)
            seq=glyphs(t.text,mode)
            if seq: counts['word'][p].update(zip(('^',)+seq,seq+('$',)))
        start=None
        for i in range(len(ts)+1):
            ok=i<len(ts) and good[i]
            if ok and start is None:start=i
            if not ok and start is not None:
                seq=glyphs(''.join(t.text for t in ts[start:i]),mode)
                if seq:
                    counts['line'][p].update(zip(seq,seq[1:]))
                    if start==0:counts['line'][p][('^',seq[0])]+=1
                    if i==len(ts):counts['line'][p][(seq[-1],'$')]+=1
                start=None
    assert all(words[p]>0 for p in selected),sorted(selected-set(words))
    return counts,dict(readable_words=dict(words),loci=dict(loci),excluded_words=dict(excluded),
        missing_pages=sorted(selected-all_pages),extra_paragraph_pages=sorted(all_pages-selected),
        observed_characters=sorted(chars),source_lines=len(doc.raw.splitlines()),source_sha256=sha(path))

def dense(counts,meta,alphabet):
    ix={g:i for i,g in enumerate(alphabet)}
    a=np.zeros((len(meta),len(ix),len(ix)),dtype=np.int32)
    for i,m in enumerate(meta):
        for (g,h),n in counts[m['id']].items():a[i,ix[g],ix[h]]=n
    return a

def distances(x):
    x=np.asarray(x,dtype=float);tot=x.sum(2);share=tot/np.maximum(tot.sum(1)[:,None],1)
    p=x/np.maximum(tot[:,:,None],1);d=np.zeros((len(x),len(x)))
    for i in range(len(x)):
        middle=(p[i]+p)/2
        js=(rel_entr(p[i],middle)+rel_entr(p,middle)).sum(2)/(2*np.log(2))
        w=(share[i]+share)/2*((tot[i]>=3)&(tot>=3))
        d[i]=(js*w).sum(1)
    np.fill_diagonal(d,0)
    assert np.isfinite(d).all() and np.min(d)>=0 and np.max(abs(d-d.T))<1e-14
    return d

def main():
    page=json.loads((OLD/'ensemble_paragraph207_point_l010.json').read_text())['design']
    pool=json.loads((OLD/'ensemble_pool50_point_l010.json').read_text())['design']
    meta=page['meta'];pm=pool['meta'];blocks=np.array(page['blocks'])
    common=module('rf_mds',TEAM/'whole_manuscript_topology_20261015/common.py')
    repro=module('rf_reference_js',DATA/'reproduce.py')
    def pooled(x):return np.array([x[blocks==m['bifolium']].sum(0) for m in pm])
    # Verify all original ZL counts, including unreadable barriers and line ends.
    zl,zl_info=count_source(ROOT/'ZL_ivtff_raw.txt',meta,'zattera',verify=True)
    alphabet=[r['unit'] for r in csv.DictReader((DATA/'alphabet.csv').open(encoding='utf-8'))]
    audit={}; frozen={}
    ids=[m['profile_id'] for m in meta+pm]
    for mode in ['word','line']:
        a=dense(zl[mode],meta,alphabet);x=np.concatenate([a,pooled(a)])
        expected=np.zeros_like(x);lookup={uid:i for i,uid in enumerate(ids)};ai={g:i for i,g in enumerate(alphabet)};seen=0
        with (DATA/f'transitions-{mode}-counts.csv').open(encoding='utf-8',newline='') as f:
            for r in csv.DictReader(f):
                if r['profile_id'] in lookup:
                    expected[lookup[r['profile_id']],ai[r['from_unit']]]=[int(r['to_'+g]) for g in alphabet];seen+=1
        assert seen==len(ids)*len(alphabet)
        np.testing.assert_array_equal(x,expected)
        ds=distances(x);names,old=repro.square(f'pairwise/shared-{mode}-js.csv');ix=[names.index(uid) for uid in ids]
        error=float(abs(ds-old[np.ix_(ix,ix)]).max());assert error<1e-12
        audit[mode]=dict(counts_exact=True,distance_max_error=error)
        frozen[mode]=ds
    # RF units are plain EVA letters, not composite Zattera units.
    rf,info=count_source(HERE/'RF1b-er.txt',meta,'raw')
    assert not info['missing_pages'] and not info['extra_paragraph_pages']
    assert all(len(c)==1 and c in 'abcdefghijklmnopqrstuvwxyz' for c in info['observed_characters'])
    alpha=sorted({'^','$'}|{g for m in rf.values() for c in m.values() for pair in c for g in pair})
    arrays={};rf_ds={}
    for mode in ['word','line']:
        a=dense(rf[mode],meta,alpha);p=pooled(a);arrays['page_'+mode+'_counts']=a;arrays['pool_'+mode+'_counts']=p
        d=distances(np.concatenate([a,p]));check=repro.weighted_js(np.concatenate([a,p]))[0]
        np.testing.assert_allclose(d,check,atol=1e-12,rtol=0)
        arrays['page_'+mode+'_distance']=d[:207,:207];arrays['pool_'+mode+'_distance']=d[207:,207:];rf_ds[mode]=d
        audit['rf_'+mode]=dict(independent_js_error=float(abs(d-check).max()),transition_events_pages=int(a.sum()),transition_events_pools=int(p.sum()))
    mds={}
    with threadpool_limits(1):
        for prefix in ['page','pool']:
            d=(arrays[prefix+'_word_distance']+arrays[prefix+'_line_distance'])/2
            x,eig=common.embed(d);arrays[prefix+'_distance']=d;arrays[prefix+'_coords']=x;arrays[prefix+'_eigenvalues']=eig
            mds[prefix]=dict(positive_dimensions=x.shape[1],negative_eigenvalue_mass=float(-eig[eig<0].sum()/abs(eig).sum()))
    np.savez_compressed(HERE/'inputs.npz',**arrays)
    assert (HERE/'inputs.npz').stat().st_size<10_000_000
    for m in meta:m['zl_readable_words']=m['words'];m['words']=info['readable_words'][m['id']]
    for m in pm:m['zl_readable_words']=m['words'];m['words']=sum(p['words'] for p in meta if p['bifolium']==m['bifolium'])
    sources=[Path(__file__),HERE/'DESIGN.md',HERE/'RF1b-er.txt',ROOT/'ZL_ivtff_raw.txt',
        ROOT/'voynich_toolkit/source_lattice.py',ROOT/'voynich_toolkit/glyph_units.py',ROOT/'voynich_toolkit/eva.py',
        DATA/'reproduce.py',DATA/'profiles.csv',DATA/'alphabet.csv',
        OLD/'ensemble_paragraph207_point_l010.json',OLD/'ensemble_pool50_point_l010.json']
    sources += [DATA/f'transitions-{m}-counts.csv' for m in ['word','line']]
    sources += [DATA/f'pairwise/shared-{m}-js.csv' for m in ['word','line']]
    save('inputs.json',dict(status='complete',representation='RF1b-er / individual basic EVA characters',
        source_url='https://www.voynich.nu/data/RF1b-er.txt',alphabet=alpha,meta=meta,pool_meta=pm,
        rf_census=info,zl_census=zl_info,verification=audit,geometry=mds,
        data_sha256=sha(HERE/'inputs.npz'),source_sha256={p.relative_to(ROOT).as_posix():sha(p) for p in sources}))
    print(json.dumps(dict(status='passed',pages=len(meta),pools=len(pm),units=len(alpha),geometry=mds,checks=audit),indent=2),flush=True)

if __name__=='__main__':main()

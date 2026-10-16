"""Build the local report, forum download bundle and artifact-link checks."""
from pathlib import Path
from html import escape
import hashlib
import json
import re
import zipfile

from markdown_it import MarkdownIt
from PIL import Image

HERE = Path(__file__).resolve().parent
STYLE = '''
:root{color-scheme:light;--ink:#193246;--muted:#536b7b;--line:#dce4e9;--paper:#fff;--wash:#f4f7f9}
*{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink);font:17px/1.6 system-ui,Segoe UI,sans-serif}
main{max-width:1180px;margin:auto;padding:36px 24px 70px}h1,h2,h3{line-height:1.2;letter-spacing:-.02em}h1{font-size:36px;margin-bottom:12px}h2{font-size:25px;margin-top:42px}h3{font-size:20px}p{max-width:92ch}a{color:#075c86;text-underline-offset:3px}nav{display:flex;gap:18px;flex-wrap:wrap;margin:24px 0}.muted{color:var(--muted)}.card{background:var(--paper);border:1px solid var(--line);border-radius:12px;padding:20px;margin:22px 0}img{display:block;width:100%;height:auto}.links{display:flex;gap:20px;flex-wrap:wrap;margin:12px 0 0}.caption{font-size:15px;color:var(--muted)}table{border-collapse:collapse;max-width:100%;font-size:15px;display:block;overflow:auto}th,td{padding:9px 13px;border-bottom:1px solid var(--line);text-align:left}th{background:#eaf0f4}code{font-size:.9em;background:#eaf0f4;padding:2px 4px;overflow-wrap:anywhere}pre{white-space:pre-wrap;background:#eaf0f4;padding:14px}details{border-top:1px solid var(--line);padding:14px 0}summary{cursor:pointer;font-weight:650}footer{margin-top:40px;color:var(--muted);font-size:14px}.report{max-width:1000px}blockquote{border-left:3px solid #89acbd;padding-left:16px;margin-left:0}
@media(max-width:600px){main{padding:20px 12px}h1{font-size:28px}.card{padding:10px}body{font-size:16px}}
'''


def doc(title, body, report=False):
    return '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>'+escape(title)+'</title><style>'+STYLE+'</style></head><body><main'+(' class="report"' if report else '')+'>'+body+'</main></body></html>\n'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    required = ['ensemble_summary.json','cassidy_h_summary.json','scegret_results.json']
    for name in required:
        assert json.loads((HERE/name).read_text(encoding='utf-8'))['status'] == 'complete', name
    md = (HERE/'RESULTS.md').read_text(encoding='utf-8')
    assert 'Pending completion' not in md and '**Status: analysis running.**' not in md
    report = MarkdownIt('commonmark',{'html':False}).enable('table').render(md)
    (HERE/'report.html').write_text(doc('Golf chronology: methods and results',
        '<nav><a href="index.html">Figures and downloads</a><a href="RESULTS.md">Report source</a></nav>'+report,True),encoding='utf-8')

    figures = [
        ('ordering_ensemble','Ordering sensitivity',
         'Each horizontal strip is one observed text unit or pooled bifolium. Color intensity records how often its fitted position lands in that interval. Starts and penalties have equal design weight regardless of fit quality; the report separately compares the best recovered scores. Separate dark areas preserve competing positions. These are empirical sensitivity frequencies, not historical probabilities.'),
        ('cassidy_H_regions','H position and the Herbal ordering',
         'Red favors 3a before 3b; blue favors 3b before 3a. The star is observed H:58–65. The line approximates equal recovered scores. Hatching marks places where the preference depends on the continuation route used to move H. Only the displayed coordinates change; lambda is 0.01.'),
        ('scegret_original_golf','scEGRET continuity of the original Golf route',
         'The original scEGRET R calculation measures changes in transition features between position bins. Larger values mark sharper changes. Original bin numbers, occupancy and Golf composition remain visible, including bins too sparse for the calculation. A jump does not uniquely identify missing bifolia.')]
    body = '<p class="muted">Golf chronology · 16 October 2026</p><h1>Ordering uncertainty and the H frontier</h1><p>Published principal-graph ensembles, derivative-based continuation and scEGRET continuity screening, with explicit text-data adaptations and numerical audits.</p><nav><a href="report.html">Methods and findings</a><a href="forum_figures.zip">Download figures and captions</a></nav>'
    bundle = []  # Forum bundle contains standalone graphics and their captions.
    for stem,title,caption in figures:
        for ext in ['png','svg']:
            assert (HERE/f'{stem}.{ext}').exists(), stem
            bundle.append(f'{stem}.{ext}')
        body += f'<section id="{stem}" class="card"><h2 style="margin-top:0">{escape(title)}</h2><a href="{stem}.png"><img src="{stem}.png" alt="{escape(title)}"></a><p class="caption">{escape(caption)}</p><p class="links"><a href="{stem}.png">PNG</a><a href="{stem}.svg">SVG</a></p></section>'
    body += '<details><summary>Physical-sheet sampling control</summary><p>Paragraphs on one physical sheet are included or omitted together. This changes the sampling unit and is kept separate from the principal point-resampling figure.</p><a href="ordering_ensemble_sheet_control.png"><img src="ordering_ensemble_sheet_control.png" alt="Whole-sheet sampling control"></a><p class="links"><a href="ordering_ensemble_sheet_control.png">PNG</a><a href="ordering_ensemble_sheet_control.svg">SVG</a></p></details>'
    body += '<details><summary>Calculations and audits</summary><p><a href="ensemble_summary.json">Ensemble outcomes</a> · <a href="ensemble_findings.json">Joint order comparisons</a> · <a href="cassidy_h_summary.json">H continuation</a> · <a href="scegret_summary.json">scEGRET summary</a></p><p><a href="ensemble_audit.json">Ensemble accounting audit</a> · <a href="ensemble_independent_audit.json">Independent ensemble replays</a> · <a href="cassidy_independent_audit.json">Independent continuation audit</a> · <a href="scegret_independent_audit.json">Independent scEGRET audit</a></p><p><a href="DESIGN.md">Study design</a> · <a href="ENSEMBLE_DESIGN.md">Ensemble details</a> · <a href="CASSIDY_DESIGN.md">Continuation details</a> · <a href="REPRODUCE.md">Reproduction instructions</a></p></details>'
    body += '<footer>The geometry, feature definitions and finite set of starting routes are held fixed as described in the report. These figures quantify sensitivity within that design; they do not date the manuscript or prove a unique chronology.</footer>'
    (HERE/'index.html').write_text(doc('Golf chronology: published stability methods',body),encoding='utf-8')
    captions = '\n\n'.join(title+'\n'+caption+'\n\n'+(HERE/f'{stem}_caption.txt').read_text(encoding='utf-8').strip() for stem,title,caption in figures)
    captions += '\n\nPhysical-sheet sampling control\n'+(HERE/'ordering_ensemble_sheet_control_caption.txt').read_text(encoding='utf-8').strip()
    (HERE/'forum_captions.txt').write_text(captions+'\n\nMethods: Albergante et al. (2020), doi:10.3390/e22030296; Cassidy (2023), doi:10.1007/s11538-023-01200-0; Zheng et al. (2023), doi:10.1093/bib/bbad356. Text-specific inputs, conditioning and controls are described in the accompanying local study report.\n',encoding='utf-8')
    bundle += ['forum_captions.txt','ordering_ensemble_sheet_control.png','ordering_ensemble_sheet_control.svg']
    with zipfile.ZipFile(HERE/'forum_figures.zip','w',zipfile.ZIP_DEFLATED) as z:
        for name in bundle:
            z.write(HERE/name,name)

    links=[]
    for filename in ['index.html','report.html']:
        for target in re.findall(r'(?:href|src)="([^"]+)"',(HERE/filename).read_text(encoding='utf-8')):
            if re.match(r'https?://|#',target):
                continue
            target=target.split('#')[0]
            assert (HERE/target).is_file(), (filename,target)
            links.append(dict(source=filename,target=target))
    images=[]
    for path in HERE.glob('*.png'):
        with Image.open(path) as im:
            im.verify()
        with Image.open(path) as im:
            images.append(dict(file=path.name,width=im.width,height=im.height,sha256=sha(path)))
    check=dict(status='passed',local_links_checked=len(links),links=links,images=images,
               bundle_members=bundle,scope='Programmatic links and PNG integrity; visual review recorded separately.')
    (HERE/'presentation_verification.json').write_text(json.dumps(check,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status='complete',links=len(links),images=len(images),zip_bytes=(HERE/'forum_figures.zip').stat().st_size)))


if __name__=='__main__':
    main()

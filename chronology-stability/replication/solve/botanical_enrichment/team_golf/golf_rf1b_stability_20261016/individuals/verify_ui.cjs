// Exercise the shipped script, filters, pairs and CSV without refitting.
const fs=require('node:fs'),vm=require('node:vm'),assert=require('node:assert/strict');
const path=require('node:path'),root=__dirname;
const html=fs.readFileSync(process.argv[2]||path.join(root,'index.html'),'utf8');
const payload=html.match(/<script id="dataset" type="application\/json">([\s\S]*?)<\/script>/)[1];
const script=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const defaults={sampler:'point',penalty:'both',group:'all',search:''};
for(const id of ['cohort','rep','scope','sort']){
  const select=html.match(new RegExp(`<select id="${id}">([\\s\\S]*?)<\\/select>`))[1];
  defaults[id]=select.match(/<option value="([^"]+)" selected>/)[1];
}
assert.equal(defaults.cohort,'pool50');assert.equal(defaults.rep,'both');
assert.equal(defaults.scope,'best');assert.equal(defaults.sort,'consensus');
const elements={},downloads=[];
const element=id=>elements[id]??=( {value:defaults[id]??'',innerHTML:'',textContent:id==='dataset'?payload:'',disabled:false,querySelectorAll(){return []},click(){downloads.push(this.download)}} );
const context=vm.createContext({document:{getElementById:element,createElement:()=>element('download')},Blob,URL,setTimeout:fn=>fn(),console});
vm.runInContext(script,context);
const run=s=>vm.runInContext(s,context);
assert.equal(run('state.visible.length'),50);
assert.equal(elements.group.value,'all');
assert.match(elements.selectedTitle.textContent,/I:67-68/);
assert.equal(run('state.settings.length'),4);
assert(run('state.visible.every((r,i,a)=>i===0||r.score>=a[i-1].score)'));
run("$('example3b').onclick()");assert.equal(run('state.visible.length'),15);
assert.match(elements.selectedTitle.textContent,/f86v4/);
assert(run('state.rows.every(r=>Math.abs(r.h.reduce((a,b)=>a+b,0)-1)<1e-12)'));
run("$('example2a').onclick()");assert.equal(run('state.visible.length'),5);
assert.match(elements.selectedTitle.textContent,/I:67-68/);
run("$('sort').value='consensus';refresh()");
assert.equal(run('state.visible.at(-1).m.id'),'I:67-68');
run("$('rep').value='both';$('scope').value='best';refresh()");
assert.equal(run('state.settings.length'),4);
assert(run('state.rows.every(r=>Math.abs(r.h.reduce((a,b)=>a+b,0)-1)<1e-12)'));
assert(run('(()=>{const p=paired(0,1),q=paired(1,0);return Math.abs(p.before-q.after)<1e-12&&Math.abs(p.before+p.tie+p.after-1)<1e-12})()'));
run("$('csv').onclick()");assert.deepEqual(downloads,['individual-ordering-sensitivity.csv']);
let combinations=0;
for(const cohort of ['paragraph207','pool50'])for(const rep of ['RF','ZL','both'])for(const scope of ['all','best'])for(const penalty of ['both','0.01','0.05'])for(const sampler of ['point','sheet']){
  run(`$('cohort').value='${cohort}';$('rep').value='${rep}';$('scope').value='${scope}';$('penalty').value='${penalty}';$('sampler').value='${sampler}';$('group').value='all';refresh()`);
  assert.equal(run('state.visible.length'),cohort==='paragraph207'?207:50);
  assert(run('state.visible.every(r=>Number.isFinite(r.score)&&r.score>=0&&r.score<=1)'));
  combinations++;
}
run("$('search').value='not-a-page';refresh()");assert.equal(run('state.visible.length'),0);
const result={status:'passed',control_combinations:combinations,checks:['normalization','focal selection','paired reversal symmetry','consensus sorting','CSV action','empty search'],browser_visual_check:'separate'};
if(!process.argv[2])fs.writeFileSync(path.join(root,'ui_verification.json'),JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify(result));

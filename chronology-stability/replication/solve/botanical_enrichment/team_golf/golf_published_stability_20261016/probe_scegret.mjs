// Diagnose runtime stages using one original-function case, outside the ensemble.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { WebR } from './_r_runtime/node_modules/webr/dist/webr.mjs';
const here=path.dirname(fileURLToPath(import.meta.url));
const start=Date.now();
const log=(...x)=>console.log(((Date.now()-start)/1000).toFixed(1),...x);
const webR=new WebR({interactive:false});
await webR.init();log('R initialized');
await webR.installPackages(['matrixStats','matrixTests','reshape2','dplyr','jsonlite']);log('packages installed');
for(const name of ['methods','matrixStats','matrixTests','reshape2','dplyr','jsonlite']){
  await webR.evalRVoid(`suppressPackageStartupMessages(library(${name}))`);log('library',name);
}
for(const name of ['egret_flap.R','egret_fly.R']){
  await webR.FS.writeFile('/home/web_user/'+name,new Uint8Array(fs.readFileSync(path.join(here,'scegret_source','R',name))));
  await webR.evalRVoid(`source('/home/web_user/${name}')`);log('source',name);
}
const c=JSON.parse(fs.readFileSync(path.join(here,'scegret_inputs.json'),'utf8')).cohorts[0];
await webR.FS.writeFile('/home/web_user/counts.csv',new Uint8Array(fs.readFileSync(path.join(here,c.matrices.word_per_million))));
await webR.FS.writeFile('/home/web_user/pt.csv',new TextEncoder().encode('t\n'+c.candidates[0].positions.join('\n')+'\n'));
await webR.evalRVoid("cc <- read.csv('/home/web_user/counts.csv',row.names=1,check.names=FALSE); pt <- read.csv('/home/web_user/pt.csv')$t; mm <- data.frame(row.names=colnames(cc))");log('input loaded');
await webR.evalRVoid('obj <- egret.flap(cc,mm,pt)');log('flap complete');
await webR.evalRVoid('obj <- suppressWarnings(egret.fly(obj,steps=10,minCellCount=3,min.detect=.1))');log('fly complete');
log('result',await webR.evalRString('jsonlite::toJSON(obj@H_Summary,auto_unbox=TRUE,digits=16)'));
await webR.close();

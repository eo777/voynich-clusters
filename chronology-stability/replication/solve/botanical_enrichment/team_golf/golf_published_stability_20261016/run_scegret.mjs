// Execute pinned, unmodified scEGRET R functions under the local WebR interpreter.
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import { WebR } from './_r_runtime/node_modules/webr/dist/webr.mjs';

const here=path.dirname(fileURLToPath(import.meta.url));
const segment=process.argv[2] || 'all';
if(!['all','pages_remaining','pools','page_line_per_million','page_line_raw','pool_lines'].includes(segment))throw Error('Unknown segment');
const webR=new WebR({interactive:false});
console.log('Starting WebR for original scEGRET functions');
await webR.init();
console.log('R ready',await webR.evalRString('R.version.string'));
await webR.installPackages(['matrixStats','matrixTests','reshape2','dplyr','jsonlite']);
console.log('scEGRET computation dependencies installed');
for(const name of ['egret_flap.R','egret_fly.R']){
  await webR.FS.writeFile('/home/web_user/'+name,new Uint8Array(fs.readFileSync(path.join(here,'scegret_source','R',name))));
}
await webR.evalRVoid(`
suppressPackageStartupMessages({library(methods);library(matrixStats);library(matrixTests);library(reshape2);library(dplyr);library(jsonlite)})
source('/home/web_user/egret_flap.R')
source('/home/web_user/egret_fly.R')
`);
const inputs=JSON.parse(fs.readFileSync(path.join(here,'scegret_inputs.json'),'utf8'));
const results={status:'running',method:'Unmodified egret.flap and egret.fly from pinned scEGRET source, executed in R/WebR',
  runtime:{webr:'0.6.0',R:await webR.evalRString('R.version.string')},design:inputs.design,cohorts:[]};
const target=path.join(here,segment==='all'?'scegret_results.json':`scegret_job_${segment}.json`);
const save=()=>fs.writeFileSync(target,JSON.stringify(results,null,2)+'\n');
let count=0;
for(const cohort of inputs.cohorts){
  if(segment==='pages_remaining' && cohort.cohort!=='paragraph207')continue;
  if(segment==='pools' && cohort.cohort!=='pool50')continue;
  if(segment.startsWith('page_line') && cohort.cohort!=='paragraph207')continue;
  if(segment==='pool_lines' && cohort.cohort!=='pool50')continue;
  const cc={cohort:cohort.cohort,meta:cohort.meta,runs:[]};results.cohorts.push(cc);
  for(const [matrixId,filename] of Object.entries(cohort.matrices)){
    if(segment==='pages_remaining' && matrixId==='word_per_million')continue;
    if(segment==='page_line_per_million' && matrixId!=='line_per_million')continue;
    if(segment==='page_line_raw' && matrixId!=='line_raw')continue;
    if(segment==='pool_lines' && !matrixId.startsWith('line_'))continue;
    await webR.FS.writeFile('/home/web_user/counts.csv',new Uint8Array(fs.readFileSync(path.join(here,filename))));
    await webR.evalRVoid("input.counts <- read.csv('/home/web_user/counts.csv',row.names=1,check.names=FALSE)");
    for(const candidate of cohort.candidates){
      await webR.FS.writeFile('/home/web_user/pseudotime.csv',new TextEncoder().encode('t\n'+candidate.positions.join('\n')+'\n'));
      for(const steps of inputs.design.steps[cohort.cohort]){
        const code=`
          pt <- read.csv('/home/web_user/pseudotime.csv')$t
          meta <- data.frame(row.names=colnames(input.counts))
          rawbin <- cut(pt,breaks=${steps})
          occupancy <- tabulate(as.integer(rawbin),nbins=${steps})
          if(sum(occupancy>3)<2){
            answer <- list(status='insufficient_occupied_bins',occupancy=occupancy)
          }else{
            obj <- egret.flap(input.counts,meta,pt)
            obj <- suppressWarnings(egret.fly(obj,steps=${steps},minCellCount=3,min.detect=.1))
            kept <- as.integer(names(obj@CellBins))
            pairs <- t(as.matrix(obj@AllCombns))
            contrib <- abs(as.matrix(obj@stepChange))/obj@GeneInfo$Mean
            top <- lapply(seq_len(ncol(contrib)),function(j){
              ix <- head(order(contrib[,j],decreasing=TRUE),10)
              list(features=rownames(obj@GeneInfo)[ix],contribution=contrib[ix,j]/nrow(obj@FilteredCounts))
            })
            answer <- list(status='complete',occupancy=occupancy,retained_bins=kept,
              feature_count=nrow(obj@FilteredCounts),features=rownames(obj@FilteredCounts),
              assignment=as.integer(as.character(obj@Metadata$TimeBin)),
              factor_codes=as.integer(obj@Metadata$TimeBin),
              adjacent=obj@H_Summary,adjacent_bin_pairs=cbind(head(kept,-1),tail(kept,-1)),
              full_pairs=pairs,full_enthalpy=obj@H_Summaryfull$StepEnthalpy/nrow(obj@FilteredCounts),
              top_adjacent_features=top,bin_observations=obj@CellBins)
          }
          jsonlite::toJSON(answer,auto_unbox=TRUE,na='null',digits=16)
        `;
        try{
          const answer=JSON.parse(await webR.evalRString(code));
          cc.runs.push({candidate:candidate.id,matrix:matrixId,steps,...answer});
        }catch(error){
          cc.runs.push({candidate:candidate.id,matrix:matrixId,steps,status:'error',error:String(error)});
          save();throw error;
        }
        count++;
        save();
        if(count%9===0)console.log('Progress',segment,cohort.cohort,matrixId,count);
      }
    }
    save();console.log('Completed',cohort.cohort,matrixId,'total runs',count);
  }
}
results.status='complete';results.run_count=count;
results.runtime.packages=JSON.parse(await webR.evalRString("jsonlite::toJSON(as.data.frame(installed.packages()[,c('Package','Version')]),auto_unbox=TRUE)"));
save();console.log('Complete scEGRET runs',count);
await webR.close();

#!/usr/bin/env python3
"""Verify release bytes, optionally replay sealed computation in a new directory."""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def verify():
    manifest=json.loads((ROOT/'MANIFEST_SHA256.json').read_text())
    for name,digest in manifest['sha256'].items():
        if sha(ROOT/name)!=digest:raise RuntimeError('Changed release bytes: '+name)
    forecasts=json.loads((ROOT/'FORECASTS.json').read_text())
    for name,digest in forecasts['source_sha256'].items():
        if sha(ROOT/name)!=digest:raise RuntimeError('Changed preregistered source: '+name)
    receipt=json.loads((ROOT/'PREREG_RECEIPT.json').read_text())
    public=json.loads((ROOT/'PUBLIC_PREREG_RECEIPT.json').read_text())
    if not sha(ROOT/'FORECASTS.json')==receipt['forecast_sha256']==public['forecast_sha256']:
        raise RuntimeError('Forecast receipts disagree')
    execution=json.loads((ROOT/'EXECUTION.json').read_text())
    for name,digest in execution['seed_files'].items():
        if sha(ROOT/'runs'/name)!=digest:raise RuntimeError('Changed seed checkpoint: '+name)
    metrics=json.loads((ROOT/'METRICS.json').read_text())
    expected='GO_ON_INDEPENDENT_REVIEW_OF_THIS_PILOT' if all(metrics['gates'].values()) else 'NO_GO_ON_FINITE_DOSE_PREDICTIVE_TRANSFER'
    if expected!=metrics['scientific_decision']:raise RuntimeError('Decision contradicts gates')
    return {'status':'BYTE_AND_RECEIPT_VERIFIED','files':len(manifest['sha256']),
            'forecast_sha256':sha(ROOT/'FORECASTS.json'),'decision':expected}


def replay():
    import numpy,scipy
    if numpy.__version__!='2.4.3' or scipy.__version__!='1.17.0':
        raise RuntimeError('Use requirements.txt pinned numpy/scipy')
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    dest=ROOT/'reproduction_runs'/stamp;dest.mkdir(parents=True,exist_ok=False)
    forecasts=json.loads((ROOT/'FORECASTS.json').read_text())
    for name in list(forecasts['source_sha256'])+['FORECASTS.json','PREREG_RECEIPT.json']:
        shutil.copy2(ROOT/name,dest/name)
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MPLBACKEND='Agg')
    for name,args in [('run_experiment.py',['run']),('analyse.py',[])]:
        with (dest/(name+'.stdout.txt')).open('w') as stdout,(dest/(name+'.stderr.txt')).open('w') as stderr:
            subprocess.run([sys.executable,str(dest/name),*args],check=True,env=env,stdout=stdout,stderr=stderr)
    new=json.loads((dest/'METRICS.json').read_text())
    old=json.loads((ROOT/'METRICS.json').read_text())
    if new['gates']!=old['gates'] or new['scientific_decision']!=old['scientific_decision']:
        raise RuntimeError('Replay changes decision; preserve output: '+str(dest))
    return {'status':'SEALED_REPLAY_COMPLETED','directory':str(dest),
            'same_gates_and_decision':True,'release_files_untouched':True}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--rerun',action='store_true');args=parser.parse_args()
    print(json.dumps(verify(),indent=2))
    if args.rerun:print(json.dumps(replay(),indent=2))

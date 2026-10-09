#!/usr/bin/env python3
"""Prepare/seal forecasts first, then execute the single predeclared comparison."""
import argparse
import csv
import datetime as dt
import hashlib
import json
import platform
from pathlib import Path
import sys
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from model import MLP, FixedFeatures, field, inputs

ROOT = Path(__file__).resolve().parent
SOURCES = ['config.json','DESIGN_RU.md','model.py','run_experiment.py','analyse.py']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)


def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def environment():
    return {'python':sys.version,'executable':sys.executable,'numpy':np.__version__,
            'scipy':scipy.__version__,'platform':platform.platform(),
            'learning_rule':'Euclidean mean-square gradient flow; no optimizer memory',
            'advisor_recommended_model':'GPT-6.1 Sol high',
            'actual_codex_model_effort':'not exposed by available runtime metadata; no override claimed'}


def future(model, theta, X, target, config, variation=False, refined=False):
    n = model.npar
    def rhs(t, state):
        if not variation: return field(model,state,X,target)
        F,DF = field(model,state[:n],X,target,True)
        return np.r_[F,(DF@state[n:].reshape(n,n)).ravel()]
    def event(t, state):
        e = model.data(state[:n],X)[0]-target
        return .5*(e@e-config['future_delta']**2)
    event.terminal = True; event.direction = -1
    initial = np.r_[theta,np.eye(n).ravel()] if variation else theta
    prefix = 'refinement_' if refined else ''
    sol = solve_ivp(rhs,(0,config['future_max_time']),initial,events=event,method='DOP853',
                    rtol=config[prefix+'rtol'],atol=config[prefix+'atol'])
    if not sol.success or not len(sol.t_events[0]):
        return {'status':'NOT_REACHED','last_error':float(np.linalg.norm(model.data(sol.y[:n,-1],X)[0]-target)),
                'message':sol.message,'nfev':sol.nfev,'horizon':config['future_max_time']}
    end = sol.y_events[0][0]; theta_end = end[:n]
    out = {'status':'REACHED','tau':float(sol.t_events[0][0]),'nfev':sol.nfev}
    if variation:
        q,J,_ = model.data(theta_end,X); grad = J.T@(q-target)
        denominator = float(grad@field(model,theta_end,X,target))
        if denominator >= 0: raise RuntimeError('Nontransverse first hitting')
        out.update({'event_denominator':denominator,
                    'g_tau':(-grad@end[n:].reshape(n,n)/denominator).tolist()})
    return out


def derivative_checks(model, theta, X, target):
    n = model.npar; eps = 1e-28
    _,J,_ = model.data(theta,X)
    jac_cs = []; df_cs = []
    for d in np.eye(n):
        z = theta.astype(complex)+1j*eps*d
        jac_cs.append(model.data(z,X)[0].imag/eps)
        df_cs.append(field(model,z,X,target).imag/eps)
    _,DF = field(model,theta,X,target,True)
    jac_err = float(np.max(np.abs(J-np.array(jac_cs).T)))
    df_err = float(np.max(np.abs(DF-np.array(df_cs).T)))
    if jac_err>1e-10 or df_err>1e-10: raise RuntimeError('Analytic derivative validation failed')
    return {'Jacobian_complex_step_max_abs':jac_err,'field_derivative_complex_step_max_abs':df_err}


def prepare(config):
    if (ROOT/'runs').exists(): raise RuntimeError('Histories already exist; cannot reseal forecasts')
    if (ROOT/'FORECASTS.json').exists(): raise RuntimeError('Forecasts already sealed; do not overwrite')
    train,all_X,dirs = inputs(config); model = MLP(config['width'])
    forecasts = []; started = utc()
    for seed in config['seeds']:
        theta = model.initial(seed); q0,Ja,_ = model.data(theta,all_X)
        _,s,_ = np.linalg.svd(Ja,full_matrices=False)
        if s[-1]<1e-9:
            forecasts.append({'seed':seed,'status':'INITIAL_RANK_FAILURE','theta0':theta.tolist()}); continue
        right_inverse = Ja.T@np.linalg.solve(Ja@Ja.T,np.eye(len(all_X)))
        Pi = np.eye(model.npar)-right_inverse@Ja
        fa,da = field(model,theta,train,q0[:3]+config['history_amplitude']*dirs['u'],True)
        fb,db = field(model,theta,train,q0[:3]+config['history_amplitude']*dirs['v'],True)
        B = db@fa-da@fb; omega = Pi@B
        checks = derivative_checks(model,theta,train,q0[:3]+config['future_amplitude']*dirs['u'])
        row = {'seed':seed,'status':'PREPARED','theta0':theta.tolist(),'q0_all':q0.tolist(),
               'singular_values_all':s.tolist(),'bracket':B.tolist(),'omega_ret':omega.tolist(),
               'tangent_residual':float(np.linalg.norm(Ja@omega)),'checks':checks,'tasks':{}}
        for task,direction in dirs.items():
            target = q0[:3]+config['future_amplitude']*direction
            ref = future(model,theta,train,target,config,variation=True)
            if ref['status']=='REACHED':
                g = np.array(ref['g_tau']); slope = float(g@omega); scale = float(np.linalg.norm(omega))
                fd = {}
                for step in config['derivative_difference_steps']:
                    direction_k = omega/scale if scale else omega
                    plus = future(model,theta+step*direction_k,train,target,config,refined=True)
                    minus = future(model,theta-step*direction_k,train,target,config,refined=True)
                    if plus['status']!='REACHED' or minus['status']!='REACHED': raise RuntimeError('Derivative first-hit failed')
                    val = (plus['tau']-minus['tau'])/(2*step)*scale
                    fd[str(step)] = {'directional_derivative':val,'relative_error':abs(val-slope)/max(abs(slope),1e-9)}
                ref.update({'target':target.tolist(),'directional_derivative':slope,'finite_difference':fd,
                            'predicted_delta':{str(h):float(2*h*h*slope) for h in config['durations']}})
            row['tasks'][task] = ref
        forecasts.append(row)
        print(json.dumps({'phase':'forecast_only','seed':seed,'sigma_min_all':float(s[-1]),
                          'slopes':{k:v.get('directional_derivative') for k,v in row['tasks'].items()}}),flush=True)
    package = {'phase':'PRE_HISTORY_FORECASTS','started_utc':started,'sealed_utc':utc(),
               'environment':environment(),'source_sha256':{name:digest(ROOT/name) for name in SOURCES},
               'config':config,'initial_states':forecasts,'paired_histories_executed':False}
    write(ROOT/'FORECASTS.json',package)
    with (ROOT/'FORECASTS.csv').open('w',newline='') as f:
        writer = csv.writer(f); writer.writerow(['seed','task','h','tau0','predicted_delta','predicted_relative_delta'])
        for seed in forecasts:
            for task,ref in seed.get('tasks',{}).items():
                if ref['status']=='REACHED':
                    for h in config['durations']:
                        value=ref['predicted_delta'][str(h)];writer.writerow([seed['seed'],task,h,ref['tau'],value,value/ref['tau']])
    print('Forecast SHA256 '+digest(ROOT/'FORECASTS.json'),flush=True)


def history(model, theta, train, target0, dirs, h, steps, orientation, config):
    routes = {'plus':[('u',1),('v',1),('u',-1),('v',-1)],
              'minus':[('v',1),('u',1),('v',-1),('u',-1)]}
    theta = theta.copy(); dt = h/steps
    for name,sign in routes[orientation]:
        target = target0+sign*config['history_amplitude']*dirs[name]
        for _ in range(steps):
            f = lambda z:field(model,z,train,target)
            k1=f(theta); k2=f(theta+dt*k1/2);k3=f(theta+dt*k2/2);k4=f(theta+dt*k3)
            theta += dt*(k1+2*k2+2*k3+k4)/6
    return theta


def restore(model, theta, all_X, q0, config):
    if isinstance(model,FixedFeatures):
        q,J,_ = model.data(theta,all_X)
        end = theta-J.T@np.linalg.solve(J@J.T,q-q0)
        return {'status':'RESTORED','theta':end.tolist(),'error':float(np.linalg.norm(model.data(end,all_X)[0]-q0)),
                'time':'exact gradient-flow limit; not charged as finite-time repair','nfev':0}
    def event(t,z):return np.linalg.norm(model.data(z,all_X)[0]-q0)-config['restore_tolerance']
    event.terminal=True;event.direction=-1
    error0 = float(np.linalg.norm(model.data(theta,all_X)[0]-q0))
    if error0<=config['restore_tolerance']:
        return {'status':'RESTORED','theta':theta.tolist(),'error':error0,'time':0.,'nfev':0}
    sol=solve_ivp(lambda t,z:field(model,z,all_X,q0),(0,config['restore_max_time']),theta,
                  events=event,method='DOP853',rtol=config['rtol'],atol=config['atol'])
    end = sol.y[:,-1]; error=float(np.linalg.norm(model.data(end,all_X)[0]-q0))
    okay=sol.success and len(sol.t_events[0]) and error<=config['restore_tolerance']*1.02
    return {'status':'RESTORED' if okay else 'RESTORE_FAILED','theta':end.tolist(),'error':error,
            'time':float(sol.t[-1]),'nfev':sol.nfev,'message':sol.message}


def linear_time(model, theta, train, target, config):
    q,J,_=model.data(theta,train);vals,U=np.linalg.eigh(J@J.T/len(train));e=U.T@(q-target)
    f=lambda t:float(np.sum(e*e*np.exp(-2*vals*t))-config['future_delta']**2)
    if f(config['future_max_time'])>0:return {'status':'NOT_REACHED'}
    return {'status':'REACHED','tau':float(brentq(f,0,config['future_max_time'],xtol=1e-12))}


def run(config):
    # Invalidate prior success before even validating source hashes.
    (ROOT/'RUN_SUCCESS.json').unlink(missing_ok=True)
    seal=json.loads((ROOT/'FORECASTS.json').read_text())
    for name,h in seal['source_sha256'].items():
        if digest(ROOT/name)!=h:raise RuntimeError('Changed sealed source: '+name)
    if config!=seal['config']:raise RuntimeError('Changed predeclared config')
    expected=json.loads((ROOT/'PREREG_RECEIPT.json').read_text())
    if digest(ROOT/'FORECASTS.json')!=expected['forecast_sha256']:raise RuntimeError('Changed forecasts')
    started=utc();wall=time.monotonic();train,all_X,dirs=inputs(config);model=MLP(config['width'])
    for ref in seal['initial_states']:
        seed=ref['seed'];path=ROOT/'runs'/f'seed_{seed}.json'
        if path.exists():
            previous=json.loads(path.read_text())
            if previous.get('forecast_sha256')!=expected['forecast_sha256']:raise RuntimeError('Stale seed checkpoint')
            print('Resume completed seed '+str(seed),flush=True);continue
        out={'seed':seed,'forecast_sha256':expected['forecast_sha256'],'started_utc':utc(),'status':ref['status'],'conditions':[]}
        if ref['status']=='PREPARED':
            theta=np.array(ref['theta0']);q0=np.array(ref['q0_all']);linear=FixedFeatures(model,theta)
            for h in config['durations']:
                cond={'h':h,'nonlinear':{},'linear':{},'tasks':{}}
                for orientation in ('plus','minus'):
                    endpoint=history(model,theta,train,q0[:3],dirs,h,config['history_rk4_steps'],orientation,config)
                    result=restore(model,endpoint,all_X,q0,config);result['history_endpoint']=endpoint.tolist()
                    cond['nonlinear'][orientation]=result
                    lp=history(linear,linear.initial_theta,train,q0[:3],dirs,h,config['history_rk4_steps'],orientation,config)
                    cond['linear'][orientation]=restore(linear,lp,all_X,q0,config)
                for task,forecast in ref['tasks'].items():
                    row={'forecast':forecast,'nonlinear':{},'linear':{}}
                    if forecast['status']=='REACHED':
                        target=np.array(forecast['target'])
                        for orientation in ('plus','minus'):
                            state=cond['nonlinear'][orientation]
                            row['nonlinear'][orientation]=(future(model,np.array(state['theta']),train,target,config)
                              if state['status']=='RESTORED' else {'status':'NOT_RUN_UNMATCHED_CURRENT_RESPONSE'})
                            row['linear'][orientation]=linear_time(linear,np.array(cond['linear'][orientation]['theta']),train,target,config)
                        if all(row['nonlinear'][s]['status']=='REACHED' for s in ('plus','minus')):
                            row['observed_delta']=row['nonlinear']['plus']['tau']-row['nonlinear']['minus']['tau']
                            _,Ja,_=model.data(theta,all_X);pinv=Ja.T@np.linalg.solve(Ja@Ja.T,np.eye(len(all_X)))
                            plus=np.array(cond['nonlinear']['plus']['theta']);minus=np.array(cond['nonlinear']['minus']['theta'])
                            mismatch=model.data(plus,all_X)[0]-model.data(minus,all_X)[0]
                            g=np.array(forecast['g_tau']);normal=g@pinv
                            row['mismatch_time_linear_estimate']=float(normal@mismatch)
                            row['mismatch_time_bound']=float(np.linalg.norm(normal)*np.linalg.norm(mismatch))
                    cond['tasks'][task]=row
                if seed in config['refinement_seeds'] and h==config['primary_duration']:
                    refined={};tight=dict(config);tight['rtol']=config['refinement_rtol'];tight['atol']=config['refinement_atol']
                    states={}
                    for orientation in ('plus','minus'):
                        end=history(model,theta,train,q0[:3],dirs,h,config['refinement_rk4_steps'],orientation,config)
                        states[orientation]=restore(model,end,all_X,q0,tight)
                    for task,f in ref['tasks'].items():
                        if f['status']=='REACHED' and all(s['status']=='RESTORED' for s in states.values()):
                            times={s:future(model,np.array(z['theta']),train,np.array(f['target']),config,refined=True) for s,z in states.items()}
                            if all(z['status']=='REACHED' for z in times.values()):
                                refined[task]={'delta':times['plus']['tau']-times['minus']['tau'],'times':times}
                    cond['refinement']=refined
                out['conditions'].append(cond)
                write(ROOT/'PROGRESS.json',{'seed':seed,'h_completed':h,'started_utc':started,'updated_utc':utc()})
                print(json.dumps({'phase':'paired_histories','seed':seed,'h':h,
                    'restore_errors':{s:z['error'] for s,z in cond['nonlinear'].items()},
                    'deltas':{k:z.get('observed_delta') for k,z in cond['tasks'].items()}}),flush=True)
        out['completed_utc']=utc();write(path,out)
    receipt={'status':'EXECUTION_COMPLETED_NOT_SCIENTIFIC_ACCEPTANCE','started_utc':started,'completed_utc':utc(),
             'wall_seconds':time.monotonic()-wall,'environment':environment(),'forecast_sha256':expected['forecast_sha256'],
             'prereg_commit':expected['prereg_commit'],'seed_files':{p.name:digest(p) for p in sorted((ROOT/'runs').glob('seed_*.json'))}}
    write(ROOT/'EXECUTION.json',receipt);write(ROOT/'RUN_SUCCESS.json',receipt)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('phase',choices=['prepare','run']);args=parser.parse_args()
    config=json.loads((ROOT/'config.json').read_text())
    if np.__version__!='2.4.3' or scipy.__version__!='1.17.0':raise RuntimeError('Use pinned numpy/scipy versions from requirements.txt')
    if args.phase=='prepare':prepare(config)
    else:run(config)

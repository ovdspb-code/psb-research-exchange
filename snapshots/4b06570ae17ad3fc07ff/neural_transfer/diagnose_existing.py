#!/usr/bin/env python3
"""Post-measurement decomposition of existing endpoints; never changes forecasts."""
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def diagnose():
    forecasts=json.loads((ROOT/'FORECASTS.json').read_text())
    f_by_seed={s['seed']:s for s in forecasts['initial_states']}
    results=json.loads((ROOT/'RESULTS.json').read_text())
    rows=[]; restoration=[]
    for seed in results['runs']:
        ref=f_by_seed[seed['seed']]
        omega=np.array(ref['omega_ret'])
        for cond in seed['conditions']:
            h=cond['h']; state=cond['nonlinear']
            if not all(s['status']=='RESTORED' for s in state.values()):continue
            actual=np.array(state['plus']['theta'])-np.array(state['minus']['theta'])
            geometric=2*h*h*omega
            restoration.extend({'seed':seed['seed'],'h':h,'orientation':side,
                                'time':s['time'],'error':s['error']} for side,s in state.items())
            for task,run in cond['tasks'].items():
                if 'observed_delta' not in run:continue
                f=run['forecast'];g=np.array(f['g_tau']);tau=f['tau']
                predicted=f['predicted_delta'][str(h)]
                linear_actual=float(g@actual)
                observed=run['observed_delta']
                rows.append({'seed':seed['seed'],'h':h,'task':task,'tau0':tau,
                    'observed_relative':observed/tau,'predicted_relative':predicted/tau,
                    'post_hoc_endpoint_linear_relative':linear_actual/tau,
                    'history_remainder_relative':(predicted-linear_actual)/tau,
                    'event_remainder_relative':(linear_actual-observed)/tau,
                    'parameter_contrast_norm':float(np.linalg.norm(actual)),
                    'geometric_parameter_prediction_error':float(np.linalg.norm(actual-geometric)/np.linalg.norm(geometric))})
    summary={}
    for h in forecasts['config']['durations']:
        r=[x for x in rows if x['h']==h]
        error=np.array([x['predicted_relative']-x['observed_relative'] for x in r])
        hist=np.array([x['history_remainder_relative'] for x in r])
        event=np.array([x['event_remainder_relative'] for x in r])
        effects=np.abs([x['observed_relative'] for x in r])
        times=[x['time'] for x in restoration if x['h']==h]
        summary[str(h)]={
            'median_abs_effect_percent':float(np.median(effects)*100),
            'max_abs_effect_percent':float(np.max(effects)*100),
            'restore_time_min_median_max':np.quantile(times,[0,.5,1]).tolist(),
            'baseline_future_time_min_median_max':np.quantile([x['tau0'] for x in r],[0,.5,1]).tolist(),
            'rmse_sealed_relative_contrast':float(np.sqrt(np.mean(error**2))),
            'rmse_history_remainder_relative':float(np.sqrt(np.mean(hist**2))),
            'rmse_event_remainder_relative':float(np.sqrt(np.mean(event**2))),
            'event_remainder_rmse_over_total':float(np.linalg.norm(event)/np.linalg.norm(error)),
            'geometric_parameter_prediction_error_median':float(np.median([x['geometric_parameter_prediction_error'] for x in r]))}
    out={'status':'POST_HOC_DIAGNOSTIC_NOT_A_FORECAST_OR_NEW_GATE',
         'method':'g_tau dotted with measured restored parameter contrast; no new trajectories or fitted coefficients',
         'decomposition':'predicted-observed=(predicted-linear_actual)+(linear_actual-observed); RMSE terms are not additive',
         'summaries':summary,'rows':rows,'restoration':restoration}
    (ROOT/'DIAGNOSTICS.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':diagnose()

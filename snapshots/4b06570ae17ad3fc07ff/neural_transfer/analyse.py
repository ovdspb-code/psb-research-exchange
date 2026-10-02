#!/usr/bin/env python3
"""Predeclared state-level scoring; failures remain in the denominator."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent


def skill(predicted, observed):
    denom=float(np.mean(observed**2))
    return 1-float(np.mean((predicted-observed)**2))/denom if denom else None


def analyse():
    config=json.loads((ROOT/'config.json').read_text())
    forecasts=json.loads((ROOT/'FORECASTS.json').read_text())
    for name,digest in forecasts['source_sha256'].items():
        if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Changed sealed analysis source: '+name)
    runs=[json.loads((ROOT/'runs'/f'seed_{s}.json').read_text()) for s in config['seeds']]
    flat=[];control_max=0.;restore_max=0.;refinements=[];fdmax=0.
    for seed in runs:
        for cond in seed['conditions']:
            restore_max=max(restore_max,max(z['error'] for z in cond['nonlinear'].values()))
            for task,row in cond['tasks'].items():
                ref=row['forecast'];h=cond['h']
                okay='observed_delta' in row and ref['status']=='REACHED'
                record={'seed':seed['seed'],'h':h,'task':task,'status':'MEASURED' if okay else 'FAILED_OR_CENSORED'}
                if okay:
                    p=ref['predicted_delta'][str(h)];v=row['observed_delta'];tau=ref['tau']
                    mismatch=row['mismatch_time_bound']
                    resolution=max(10*mismatch,1e-7*tau)
                    record.update({'tau0':tau,'predicted_delta':p,'observed_delta':v,
                        'predicted_relative':p/tau,'observed_relative':v/tau,
                        'relative_prediction_error':abs(p-v)/max(abs(v),resolution),
                        'mismatch_time_bound':mismatch,'mismatch_fraction':mismatch/max(abs(v),1e-30),
                        'resolved':abs(v)>resolution,'sign_correct':bool(p*v>0 and abs(v)>resolution),
                        'tau_plus':row['nonlinear']['plus']['tau'],'tau_minus':row['nonlinear']['minus']['tau']})
                    fdmax=max(fdmax,max(z['relative_error'] for z in ref['finite_difference'].values()))
                    if all(row['linear'][s]['status']=='REACHED' for s in ('plus','minus')):
                        ld=row['linear']['plus']['tau']-row['linear']['minus']['tau']
                        lc=abs(ld)/max(row['linear']['plus']['tau'],1e-30)
                        control_max=max(control_max,lc);record['linear_control_relative_difference']=lc
                    else: control_max=float('inf')
                    if task in cond.get('refinement',{}):
                        diff=abs(cond['refinement'][task]['delta']-v)
                        refinements.append({'seed':seed['seed'],'task':task,'absolute_difference':diff,
                                            'fraction':diff/max(abs(v),resolution)})
                flat.append(record)
    summaries={}
    for h in config['durations']:
        rows=[r for r in flat if r['h']==h];valid=[r for r in rows if r['status']=='MEASURED']
        complete=len(valid)==3*len(config['seeds'])
        out={'expected_rows':3*len(config['seeds']),'measured_rows':len(valid),'failed_or_censored_rows':3*len(config['seeds'])-len(valid)}
        if valid:
            out.update({'sign_accuracy':sum(r['sign_correct'] for r in valid)/out['expected_rows'],
                        'resolved_rows':sum(r['resolved'] for r in valid),
                        'median_relative_error':float(np.median([r['relative_prediction_error'] for r in valid])),
                        'max_mismatch_fraction':max(r['mismatch_fraction'] for r in valid),
                        'median_absolute_relative_contrast':float(np.median([abs(r['observed_relative']) for r in valid]))})
        if complete:
            pred=np.array([[next(r['predicted_relative'] for r in rows if r['seed']==s and r['task']==t)
                            for t in ('u','v','z')] for s in config['seeds']])
            obs=np.array([[next(r['observed_relative'] for r in rows if r['seed']==s and r['task']==t)
                            for t in ('u','v','z')] for s in config['seeds']])
            error_seed=np.mean((pred-obs)**2,axis=1);null_seed=np.mean(obs**2,axis=1)
            improve=null_seed-error_seed
            rng=np.random.default_rng(config['bootstrap_seed'])
            samples=rng.integers(0,len(improve),size=(config['bootstrap_repeats'],len(improve)))
            boot=improve[samples].mean(axis=1)
            ps=pred-pred.mean(axis=1,keepdims=True);os=obs-obs.mean(axis=1,keepdims=True)
            out.update({'skill':skill(pred,obs),'bootstrap_mse_improvement_CI95':np.quantile(boot,[.025,.975]).tolist(),
                        'bootstrap_unit':'initial state, all three tasks as one block',
                        'shape_skill':skill(ps,os),
                        'shape_energy_fraction':float(np.sum(os**2)/np.sum(obs**2)) if np.sum(obs**2) else 0.,
                        'seed_mse_improvement':dict(zip(map(str,config['seeds']),improve.tolist()))})
        summaries[str(h)]=out
    primary=summaries[str(config['primary_duration'])]
    def enough(name,threshold):return primary.get(name) is not None and primary[name]>=threshold
    all_restore=(len(flat)==len(config['seeds'])*3*len(config['durations']) and
                 all(z['status']=='RESTORED' for run in runs for cond in run['conditions'] for z in cond['nonlinear'].values()))
    gates={
        'complete_restoration_and_future_hits':all_restore and all(s['failed_or_censored_rows']==0 for s in summaries.values()),
        'fixed_representation_negative_control':control_max<=1e-8,
        'sign_prediction':enough('sign_accuracy',config['sign_accuracy_min']),
        'quantitative_skill':enough('skill',config['skill_min']),
        'bootstrap_beats_zero':primary.get('bootstrap_mse_improvement_CI95',[0])[0]>0,
        'median_relative_error':primary.get('median_relative_error',float('inf'))<=config['median_relative_error_max'],
        'selective_shape_skill':enough('shape_skill',config['shape_skill_min']),
        'selective_shape_energy':enough('shape_energy_fraction',config['shape_energy_fraction_min']),
        'response_mismatch_does_not_explain_contrast':primary.get('max_mismatch_fraction',float('inf'))<=config['mismatch_fraction_max'],
        'predeclared_numerical_refinement':len(refinements)==3*len(config['refinement_seeds']) and all(x['fraction']<=config['refinement_error_fraction_max'] for x in refinements),
        'event_derivative_independently_checked':fdmax<.002,
    }
    metrics={'scientific_decision':'GO_ON_INDEPENDENT_REVIEW_OF_THIS_PILOT' if all(gates.values()) else 'NO_GO_ON_FINITE_DOSE_PREDICTIVE_TRANSFER',
             'execution_is_not_independent_acceptance':True,'primary_duration':config['primary_duration'],'gates':gates,
             'summaries':summaries,'restore_max_L2':restore_max,'linear_control_max_relative_difference':control_max if np.isfinite(control_max) else None,
             'max_event_derivative_FD_relative_error':fdmax,'refinements':refinements,
             'all_seed_ids':config['seeds'],'excluded_seeds':[]}
    (ROOT/'METRICS.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    (ROOT/'RESULTS.json').write_text(json.dumps({'metrics':metrics,'runs':runs},ensure_ascii=False,indent=2,allow_nan=False)+'\n')
    with (ROOT/'COMPARISON.csv').open('w',newline='') as f:
        columns=sorted(set().union(*(x.keys() for x in flat)));w=csv.DictWriter(f,fieldnames=columns)
        w.writeheader();w.writerows(flat)
    make_plot(flat,config)
    report(metrics,forecasts)
    print(json.dumps(metrics,ensure_ascii=False,indent=2))


def make_plot(rows,config):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),layout='constrained')
    colors={'u':'#2369a0','v':'#c05429','z':'#438948'}
    for ax,h in zip(axes,config['durations']):
        valid=[r for r in rows if r['h']==h and r['status']=='MEASURED']
        vals=[]
        for task in colors:
            rr=[r for r in valid if r['task']==task]
            x=np.array([r['predicted_relative'] for r in rr])*100
            y=np.array([r['observed_relative'] for r in rr])*100;vals.extend(x);vals.extend(y)
            ax.scatter(x,y,label=task,s=34,color=colors[task],alpha=.8)
        lim=max(max(map(abs,vals),default=0)*1.12,1e-6)
        ax.plot([-lim,lim],[-lim,lim],color='#777',lw=1,ls='--',zorder=0)
        ax.axhline(0,color='#bbb',lw=.7);ax.axvline(0,color='#bbb',lw=.7)
        ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xlabel='Predicted time contrast (% of baseline)',
               ylabel='Measured time contrast (% of baseline)',title=f'h = {h:g}; {len(valid)}/36 measured')
        ax.legend(title='Next task',frameon=False);ax.set_aspect('equal',adjustable='box')
    fig.suptitle('Sealed forecast versus future first-hit time after response restoration',fontsize=11)
    fig.savefig(ROOT/'prediction_validation.svg');fig.savefig(ROOT/'prediction_validation.png',dpi=160);plt.close(fig)


def report(metrics,forecasts):
    lines=['# Предсказательный перенос в MLP: результат одного фиксированного опыта','',
           'Дата: 2026-10-02. Авторский вычислительный результат; внешний аудит не заявляется.',
           '', '**Решение: '+metrics['scientific_decision']+'**','',
           '12 заранее выбранных начальных состояний, одна архитектура 2–8–1 tanh, три будущие задачи, два размера истории. '
           'Предсказано именно конечное время первого попадания после возвращения нынешних ответов на шести входах. '
           'Метод использует следующую задачу и эталонную траекторию; дешёвого универсального выбора curriculum не получено.','',
           '## Предварительная фиксация','',
           'FORECASTS.json содержит прогнозы, состояния θ0 и SHA-256 источников; PREREG_RECEIPT.json — удалённый commit до историй. '
           'Код исполнения проверяет этот seal. Итоговые состояния не использованы для корректировки прогноза.','',
           '## Объявленные критерии','', '| Критерий | Итог |','|---|---|']
    lines += ['| '+name+' | '+('PASS' if okay else 'FAIL')+' |' for name,okay in metrics['gates'].items()]
    lines += ['', '## Числа','']
    for h,summary in metrics['summaries'].items():
        lines += ['### h='+h,'','```json',json.dumps(summary,ensure_ascii=False,indent=2),'```','']
    lines += ['## Что контролировалось','',
       '- Нынешняя функция: все шесть объявленных ответов; максимальная L2 невязка '+str(metrics['restore_max_L2'])+'.',
       '- Представление заморожено в линейном контроле; максимальная относительная разность времён '+str(metrics['linear_control_max_relative_difference'])+'.',
       '- Остаточное несовпадение: отдельная нормальная диагностическая оценка; она не добавлена к предварительному прогнозу.',
       '- Общая скорость: удалена общая компонента относительных контрастов внутри каждого состояния.',
       '- Статистическая единица: seed; все три задачи остаются одним bootstrap-блоком. Исключённых seeds нет.',
       '- Численный контроль: заранее назначены seeds 101,102; шаг истории уменьшен вдвое, будущая интеграция ужесточена.','',
       '## Граница вывода','',
       'Это обычный градиентный поток небольшой нелинейной сети на конечном синтетическом наборе. '
       'Равенство функций вне шести входов не проверялось. Нет шумового старения, биологического сна, '
       'выгодности обслуживания или восстановления всего класса задач. Общий механизм порядка, '
       'проекция возврата и чувствительность события известны. Неизвестность мирового приоритета '
       'не снимается успешным прогнозом. Непрохождение конечной дозы не спасается диагностикой h=.10.','',
       'Предшественник: [Sweeney, 2026](https://proceedings.mlr.press/v306/sweeney26a.html). '
       'Его observable — целевая потеря после порядка обновлений; здесь измерено последующее время после '
       'возвращения текущих ответов. Подстановка этого функционала в известный аппарат возможна: '
       'различие observable само по себе не доказывает новую теорему.','',
       '## Следующий шаг','',
       'При GO — независимая содержательная оценка этого фиксированного пилота и его остатка относительно литературы. '
       'При NO-GO — сохранить все числа, назвать проваленный критерий и завершить данную атаку без смены архитектуры/целей.']
    (ROOT/'REPORT_RU.md').write_text('\n'.join(lines)+'\n')


if __name__=='__main__':analyse()

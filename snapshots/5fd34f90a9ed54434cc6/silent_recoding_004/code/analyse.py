#!/usr/bin/env python3
"""004 only: deterministic design, independent R1 sums and finite binomial math."""
from pathlib import Path
import csv
import hashlib
import json
import math
import platform
import sys
import numpy
import scipy
from scipy.stats import beta, binom
from four_events import kl, confidence_endpoint, certify_counts, planned_bounds, sufficient_n

ROOT = Path(__file__).resolve().parents[1]
SOURCE, OUT = ROOT/'source', ROOT/'results'
CASE = 'D20_M23_GC102__seq7__seq6'
ROLES = ('ancestor', 'current', 'rewrite')


def read(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))


def write(name, rows):
    with (OUT/name).open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def js(name, data):
    (OUT/name).write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def cp(k, n, side, eta=.025/4):
    if side == 'lower':
        return 0. if k == 0 else float(beta.ppf(eta, k, n-k+1))
    return 1. if k == n else float(beta.ppf(1-eta, k+1, n-k))


def exact_power(probs, n, bias):
    out, tails, counts = [], [], []
    for p, side in zip(probs, ('lower', 'upper', 'lower', 'upper')):
        obs = max(0., p-bias) if side == 'lower' else min(1., p+bias)
        if side == 'lower':
            k = int(binom.ppf(.05/4, n, obs))
            tails.append(float(binom.cdf(k-1, n, obs)))
            out.append(max(0., cp(k, n, side)-bias))
        else:
            k = int(binom.ppf(1-.05/4, n, obs))
            tails.append(float(binom.sf(k, n, obs)))
            out.append(min(1., cp(k, n, side)+bias))
        counts.append(k)
    loss, restore = out[0]-out[1], out[2]-out[3]
    assert sum(tails) <= .05 + 1e-12
    return {'n': n, 'loss_lower': loss, 'restore_lower': restore,
            'conditional_power_lower': 1-sum(tails), 'tail_counts': counts,
            'pass': loss > 0 and restore >= 0}


def exact_grid(probs, bias):
    if min(probs[0]-probs[1], probs[2]-probs[3]) <= 4*bias:
        return None
    # First passing 1000-grid point, NOT a discrete global minimum theorem.
    for n in range(1000, 1000001, 1000):
        result = exact_power(probs, n, bias)
        if result['pass']:
            return result
    return None


def frequency_n(probs, bias=0.):
    def good(n):
        bounds = [confidence_endpoint(p, math.log(160)/n, side) + (-bias if side=='lower' else bias)
                  for p, side in zip(probs, ('lower', 'upper', 'lower', 'upper'))]
        return bounds[0] > bounds[1] and bounds[2] >= bounds[3]
    if min(probs[0]-probs[1], probs[2]-probs[3]) <= 2*bias:
        return None
    lo, hi = 1, 1000
    while not good(hi):
        lo, hi = hi, 2*hi
        if hi > 10**10:
            return None
    while hi-lo > 1:
        mid = (lo+hi)//2
        if good(mid): hi = mid
        else: lo = mid
    return hi


def quantities(weights, selected, unresolved_extra=frozenset()):
    masses = {}
    for role in ROLES:
        grouped = {'B': [], 'D': [], 'U': []}
        for row in weights:
            key = (int(row['position_chain']), row['target_aa'])
            if row['stop'] == 'True': group = 'D'
            elif key not in selected or row['expected_joint_status'] == '' or key in unresolved_extra: group = 'U'
            elif row['expected_joint_status'] == '1.0': group = 'B'
            else:
                assert row['expected_joint_status'] == '0.0'
                group = 'D'
            grouped[group].append(float(row['w2_'+role]))
        # NO_EXIT included in D, not removed or renormalized.
        noexit = 1-math.fsum(float(r['w2_'+role]) for r in weights)
        assert noexit >= -1e-12
        grouped['D'].append(max(0., noexit))
        masses[role] = {g: math.fsum(v) for g, v in grouped.items()}
        assert abs(math.fsum(masses[role].values())-1) < 1e-12
    probs = [masses['ancestor']['B'], masses['current']['B']+masses['current']['U'],
             masses['rewrite']['B'], masses['ancestor']['B']+masses['ancestor']['U']]
    return masses, probs


def controls():
    # Independently computed exact binomial masses; no simulation or RNG.
    rows = []
    for method in ('KL', 'CP'):
        worst = 0.
        for n in (1, 5, 20, 80):
            for p in (0., .0001, .001, .01, .03, .1, .5, .9, .99, .9999, 1.):
                pmf = [math.comb(n,k)*p**k*(1-p)**(n-k) for k in range(n+1)]
                for side in ('lower', 'upper'):
                    endpoints = [confidence_endpoint(k/n, math.log(160)/n, side) if method=='KL' else cp(k,n,side)
                                 for k in range(n+1)]
                    failure = math.fsum(v for v, b in zip(pmf, endpoints)
                                       if (b > p+2e-14 if side=='lower' else b < p-2e-14))
                    worst = max(worst, failure)
                    assert failure <= .025/4 + 1e-12
        rows.append({'control': method+'_exact_small_n_marginal_coverage',
                     'worst_noncoverage_on_declared_grid': worst, 'limit': .025/4, 'status': 'PASS'})
    good = {'a_B': (34229,10**6), 'x_BU': (27847,10**6), 'R_B': (41787,10**6), 'a_BU': (35528,10**6)}
    null = {e: (30000,10**6) for e in good}
    for method in ('KL','CP'):
        assert certify_counts(good,method=method)['arithmetic_certificate']
        assert not certify_counts(null,method=method)['arithmetic_certificate']
    bad = [dict(good, a_B=(40000,10**6)), dict(good, a_BU=(35528,999999)),
           dict(good, x_BU=(1.1,100)), dict(good, x_BU=(True,100)),
           dict(good, R_B=(101,100)), dict(good, R_B=(0,0)),
           {e:v for e,v in good.items() if e != 'a_B'}]
    for method in ('KL','CP'):
        for sample in bad:
            try: certify_counts(sample,method=method)
            except ValueError: pass
            else: raise AssertionError('Invalid input accepted')
    rows.append({'control': 'positive_planted_null_and_invalid_counts',
                 'worst_noncoverage_on_declared_grid': None, 'limit': None, 'status': 'PASS'})
    return rows


def main():
    OUT.mkdir(exist_ok=True)
    (OUT/'RUN_SUCCESS.json').unlink(missing_ok=True)
    lock = json.loads((ROOT/'SOURCE_LOCK.json').read_text())
    for name, digest in lock['sha256'].items():
        assert hashlib.sha256((SOURCE/name).read_bytes()).hexdigest() == digest, name
    weights = [r for r in read(SOURCE/'all_destination_weights.csv') if r['case_id']==CASE]
    old = read(SOURCE/'planned_phenotype_panel.csv')
    fixed = {(int(r['position_chain']), r['target_aa']) for r in old if r['role']=='FOCAL_RESERVE'}
    assert len(weights)==2140 and len(fixed)==85
    candidates = [r for r in weights if r['stop']!='True' and (int(r['position_chain']),r['target_aa']) not in fixed]
    candidates.sort(key=lambda r:(-max(float(r['w2_'+s]) for s in ROLES),int(r['position_chain']),r['target_aa']))
    # Sign-flipping expected statuses must not change the chosen order.
    rankkeys = [(int(r['position_chain']),r['target_aa']) for r in candidates]
    reference = json.loads((SOURCE/'advisor_reference/PROJECTED_DESIGN.json').read_text())
    advisor = {r['proposed_total_focal_statuses']:r for r in reference['frontier']}
    frontier, plans, exact, exact_rounded, sensitivity, panel_comparison, families = [], [], [], [], [], [], []
    for size in (85,600,650,800):
        selected = fixed | set(rankkeys[:size-85])
        chosen = [r for r in weights if (int(r['position_chain']),r['target_aa']) in selected]
        if size != 85:
            other = read(SOURCE/f'advisor_reference/PROPOSED_PANEL_{size}.csv')
            assert {(int(r['position_chain']),r['target_aa']) for r in other} == selected
            panel_comparison.append({'size':size,'same_keys_as_advisor':True,'retained_85':fixed<=selected})
        inventory = []
        for row in chosen:
            inventory.append({'position_chain':row['position_chain'],'current_aa':row['current_aa'],
                              'target_aa':row['target_aa'],'in_retained_85':(int(row['position_chain']),row['target_aa']) in fixed,
                              'planning_status_only':row['expected_joint_status'] or 'UNKNOWN',
                              'actual_status':'NOT_MEASURED','priority_probability_max':max(float(row['w2_'+s]) for s in ROLES)})
        write(f'PROPOSED_PANEL_{size}.csv',inventory)
        masses, probs = quantities(weights, selected)
        zero = sufficient_n(probs)
        biased = sufficient_n(probs,.0001)
        freq = frequency_n(probs)
        assert zero == advisor[size]['n_per_start_sufficient_for_95pct_conditional_power_zero_bias']
        assert biased == advisor[size]['n_per_start_sufficient_for_95pct_conditional_power_bias_1e_4']
        assert freq == advisor[size]['n_per_start_at_model_frequencies_zero_bias']
        g_loss, g_restore = probs[0]-probs[1], probs[2]-probs[3]
        for i,s in enumerate(ROLES):
            assert abs(masses[s]['U']-advisor[size]['unresolved_mass'][i])<1e-12
        objects = size+19+3
        frow = {'panel_focal':size,'additional_focal_vs_85':size-85,'protein_objects_including_QA_baselines':objects,
                'absolute_endpoint_family':2*objects,'alpha_F':.025,'alpha_F_per_endpoint':.025/(2*objects),
                'unknown_current_model_mass':masses['current']['U'],'coarse_loss_model_lower':g_loss,
                'coarse_restore_model_lower':g_restore,'KL_n_at_model_frequency_NOT_power':freq,
                'KL_n_sufficient_95pct_conditional_power':zero,'KL_n_same_power_bias_1e_4':biased,
                'independent_unit':'iid ordered histories per start, 3 starts; NOT reads/clonal descendants'}
        frontier.append(frow)
        for b in (0.,.0001,.00025,.0005,.001,.0015,.002):
            sensitivity.append({'panel':size,'certified_bias_per_start':b,'KL_n_for_95pct_conditional_power':sufficient_n(probs,b),
                                'limiting_worst_model_loss':g_loss-4*b,'limiting_worst_model_restore':g_restore-4*b})
        for n in ((310000,370000) if size==600 else (150000,170000) if size==650 else (150000,170000) if size==800 else (150000,)):
            for b in (0.,.0001):
                p=planned_bounds(probs,n,b)
                plans.append({'panel':size,'n_per_start':n,'bias':b,**{k:v for k,v in p.items() if k!='bounds'}})
        if size!=85:
            for b in (0.,.0001):
                r=exact_grid(probs,b)
                exact.append({'panel':size,'bias':b,'method':'exact CP confidence + exact binomial power quantiles',
                              'first_passing_n_on_1000_grid_NOT_global_minimum':r['n'],
                              'loss_lower':r['loss_lower'],'restore_lower':r['restore_lower'],
                              'conditional_power_lower_bound':r['conditional_power_lower'],
                              'KL_n_comparator':zero if b==0 else biased})
            for n in ((180000,220000) if size==600 else (90000,100000,110000)):
                for b in (0.,.0001):
                    r=exact_power(probs,n,b)
                    exact_rounded.append({'panel':size,'n':n,'bias':b,
                                          **{k:v for k,v in r.items() if k not in ('tail_counts','n')}})
        if size==650:
            for row in chosen:
                if row['expected_joint_status']=='1.0':
                    key=(int(row['position_chain']),row['target_aa'])
                    changed, pp=quantities(weights,selected,{key})
                    sensitivity.append({'panel':650,'certified_bias_per_start':0.,
                                        'KL_n_for_95pct_conditional_power':sufficient_n(pp),
                                        'limiting_worst_model_loss':pp[0]-pp[1],
                                        'limiting_worst_model_restore':pp[2]-pp[3],
                                        'one_planned_useful_status_unresolved':str(key)})
        controls_rows = [r for r in old if r['role']!='FOCAL_RESERVE']
        assert len(controls_rows)==19
        # The no-op background has no next-variant panel rows, but its common
        # baseline remains in the declared family. Read all starting roles.
        baselines = sorted({r['protein_id'] for r in read(SOURCE/'planned_starting_constructs.csv')})
        assert len(baselines)==3
        pid = weights[0]['protein_id']
        objects_list = [(f'{pid}:LC{r["position_chain"]}{r["target_aa"]}',pid,'FOCAL',f'{pid}:BASELINE') for r in chosen]
        objects_list += [(f'{r["protein_id"]}:LC{r["position_chain"]}{r["target_aa"]}',r['protein_id'],'QA',f'{r["protein_id"]}:BASELINE') for r in controls_rows]
        objects_list += [(p+':BASELINE',p,'BASELINE',p+':BASELINE') for p in baselines]
        assert len({o[0] for o in objects_list})==objects
        family = [{'object_id':o,'protein_background':p,'role':role,'baseline_id':bl,
                   'endpoint':endpoint,'lower':None,'upper':None,'state':'NOT_MEASURED',
                   'max_marginal_error':.025/(2*objects)}
                  for o,p,role,bl in objects_list for endpoint in ('binding_log10Ka','expression_ln_assay')]
        write(f'ABSOLUTE_ENDPOINT_FAMILY_{size}.csv',family)
        families.append({'panel':size,'objects':objects,'endpoints':len(family),'baseline_objects':3,
                         'state':'TEMPLATE_ONLY; valid assay intervals and scale calibration required'})
        print('panel',size,'KL conditional n',zero,'bias n',biased,flush=True)
    checks=controls()
    write('COST_FRONTIER.csv',frontier)
    write('PLANNED_POWER_AT_ROUNDED_N.csv',plans)
    write('EXACT_BINOMIAL_COMPARISON.csv',exact)
    write('EXACT_POWER_AT_ROUNDED_N.csv',exact_rounded)
    # Heterogeneous rows intentionally use a union of keys, not discarded fields.
    keys=list(dict.fromkeys(k for row in sensitivity for k in row))
    write('SENSITIVITY.csv',[{k:r.get(k) for k in keys} for r in sensitivity])
    js('COST_FRONTIER.json',{'scope':'before-new-data planning; not optimality or biological evidence','frontier':frontier,
                           'exact_binomial_step':exact,'model_masses':{s:quantities(weights,fixed|set(rankkeys[:s-85]))[0] for s in (85,600,650,800)}})
    js('AUDIT_CHECKS.json',{'marginal_coverage_checks':checks,'exact_advisor_panel_comparisons':panel_comparison,
                          'counts_fields_are_history_counts':True,'zero_exit_preserved':True,'source_unknowns_not_zeroed':True,
                          'old_R1_or_003_runner_executed':False,'mutation_simulations':0,'physical_data_generated':False})
    js('FUNCTIONAL_FAMILIES.json',families)
    js('ENVIRONMENT.json',{'python':sys.version,'platform':platform.platform(),'scipy':scipy.__version__,'numpy':numpy.__version__,
                          'model':'not exposed by this session metadata','effort':'not exposed by this session metadata'})
    js('RUN_SUCCESS.json',{'exit':0,'status':'NEW_004_CALCULATION_COMPLETE','biological_measurements':False,
                          'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                          'statistical_module_sha256':hashlib.sha256((ROOT/'code/four_events.py').read_bytes()).hexdigest()})


if __name__=='__main__':
    main()

"""Fixed-n one-sided Bernoulli KL bounds. No biological attestation implied."""
import math

EVENTS = ('a_B', 'x_BU', 'R_B', 'a_BU')
SIDES = ('lower', 'upper', 'lower', 'upper')


def kl(q, p):
    if not 0 <= q <= 1 or not 0 <= p <= 1:
        raise ValueError('Bernoulli probabilities outside [0,1]')
    if p == 0:
        return 0.0 if q == 0 else math.inf
    if p == 1:
        return 0.0 if q == 1 else math.inf
    first = q * math.log(q / p) if q else 0.0
    second = (1-q) * math.log((1-q)/(1-p)) if q < 1 else 0.0
    return max(0.0, first + second)


def confidence_endpoint(q, c, side):
    """Invert KL(q,p)<=c in p; conservative outward rounding."""
    if side not in ('lower', 'upper') or c <= 0 or not math.isfinite(c):
        raise ValueError('Invalid side or confidence budget')
    if not 0 <= q <= 1:
        raise ValueError('Invalid empirical probability')
    if side == 'lower' and q == 0:
        return 0.0
    if side == 'upper' and q == 1:
        return 1.0
    lo, hi = (0.0, q) if side == 'lower' else (q, 1.0)
    for _ in range(80):
        p = (lo + hi) / 2
        if (kl(q, p) > c) == (side == 'lower'):
            lo = p
        else:
            hi = p
    return lo if side == 'lower' else hi


def frequency_tail(p, c, side):
    """Invert KL(q,p)<=c in q, not in p. Used ONLY for planning power."""
    if not 0 <= p <= 1 or side not in ('lower', 'upper') or c <= 0:
        raise ValueError('Invalid tail input')
    if p in (0, 1):
        return float(p)
    lo, hi = (0.0, p) if side == 'lower' else (p, 1.0)
    if side == 'lower' and kl(0, p) <= c:
        return 0.0
    if side == 'upper' and kl(1, p) <= c:
        return 1.0
    for _ in range(80):
        q = (lo + hi) / 2
        if (kl(q, p) > c) == (side == 'lower'):
            lo = q
        else:
            hi = q
    return lo if side == 'lower' else hi


def certify_counts(counts, alpha_supply=0.025, bias=0.0, method='KL'):
    """Four declared event counts. Return algebraic bounds, not empirical rescue.

    Assumptions must be independently established: valid functional intervals,
    fixed bins independent of these supply observations, iid histories within
    start, fixed n, correctly observed first protein exit and certified bias.
    Ancestor B and BU use the SAME sample (not two separate n).
    """
    if set(counts) != set(EVENTS) or not 0 < alpha_supply < 1 or method not in ('KL','CP'):
        raise ValueError('Require exactly four declared events and valid alpha')
    if not math.isfinite(bias) or not 0 <= bias < 1:
        raise ValueError('Invalid certified bias budget')
    for k, n in counts.values():
        if type(k) is not int or type(n) is not int or n <= 0 or not 0 <= k <= n:
            raise ValueError('Counts must be integer histories with 0<=k<=n')
    if counts['a_B'][1] != counts['a_BU'][1] or counts['a_B'][0] > counts['a_BU'][0]:
        raise ValueError('Ancestor events must be nested on the same sample')
    bounds = []
    for event, side in zip(EVENTS, SIDES):
        k, n = counts[event]
        if method == 'KL':
            endpoint = confidence_endpoint(k/n, math.log(4/alpha_supply)/n, side)
        else:
            from scipy.stats import beta
            eta = alpha_supply/4
            endpoint = (0. if k==0 else float(beta.ppf(eta,k,n-k+1))) if side=='lower' else (1. if k==n else float(beta.ppf(1-eta,k+1,n-k)))
        bounds.append(max(0., endpoint-bias) if side == 'lower' else min(1., endpoint+bias))
    loss = bounds[0] - bounds[1]
    restore = bounds[2] - bounds[3]
    return {'method':method,'bounds': dict(zip(EVENTS, bounds)), 'loss_lower': loss,
            'restore_lower': restore, 'improvement_lower': loss+restore,
            'arithmetic_certificate': bool(loss > 0 and restore >= 0),
            'meaning': 'Conditional algebraic result; no assay/independence/biological attestation'}


def planned_bounds(probs, n, bias=0.0, alpha_supply=.025, failure_power=.05):
    if len(probs) != 4 or type(n) is not int or n <= 0 or not 0 <= bias < 1:
        raise ValueError('Invalid planning input')
    bounds = []
    for p, side in zip(probs, SIDES):
        # Worst allowed OBSERVATION law, then correction back to target law.
        observed = max(0., p-bias) if side == 'lower' else min(1., p+bias)
        q = frequency_tail(observed, math.log(4/failure_power)/n, side)
        bound = confidence_endpoint(q, math.log(4/alpha_supply)/n, side)
        bounds.append(max(0., bound-bias) if side == 'lower' else min(1., bound+bias))
    loss, restore = bounds[0]-bounds[1], bounds[2]-bounds[3]
    return {'loss_power_lower': loss, 'restore_power_lower': restore,
            'both_pass_sufficient_95pct_conditional_power': loss > 0 and restore >= 0,
            'bounds': bounds}


def sufficient_n(probs, bias=0.0, cap=10**9):
    # No finite n helps if the limiting conservative contrasts are negative.
    if min(probs[0]-probs[1], probs[2]-probs[3]) <= 4*bias:
        return None
    lo, hi = 1, 1000
    while hi < cap and not planned_bounds(probs, hi, bias)['both_pass_sufficient_95pct_conditional_power']:
        lo, hi = hi, min(cap, hi*2)
    if not planned_bounds(probs, hi, bias)['both_pass_sufficient_95pct_conditional_power']:
        return None
    while hi-lo > 1:
        mid = (lo+hi)//2
        if planned_bounds(probs, mid, bias)['both_pass_sufficient_95pct_conditional_power']:
            hi = mid
        else:
            lo = mid
    return hi

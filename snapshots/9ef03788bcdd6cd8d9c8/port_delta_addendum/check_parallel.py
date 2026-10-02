#!/usr/bin/env python3
"""Only the new parallel-splitting algebra, stdlib exact Fraction."""
from fractions import Fraction as Q
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
checks = []


def verify(name, condition):
    if not condition:
        raise AssertionError(name)
    checks.append(name)


def gram(cols):
    return [[sum(c[i]*c[j] for c in cols) for j in range(2)] for i in range(2)]


def solve(matrix, rhs):
    a = [list(row)+[val] for row,val in zip(matrix,rhs)]
    for i in range(len(a)):
        pivot = next(j for j in range(i,len(a)) if a[j][i])
        a[i],a[pivot] = a[pivot],a[i]
        scale = a[i][i]
        a[i] = [x/scale for x in a[i]]
        for j in range(len(a)):
            if j != i:
                scale = a[j][i]
                a[j] = [x-scale*y for x,y in zip(a[j],a[i])]
    return [row[-1] for row in a]


def physical(edges):
    # Independent Kirchhoff solve and differentiated node balances.
    free = ['h','u','v']; pins = {'s':Q(1),'g':Q(0)}
    L = [[Q(0) for _ in free] for _ in free]; rhs = [Q(0)]*3
    for x,y,k in edges:
        for node,other in ((x,y),(y,x)):
            if node in free:
                i = free.index(node); L[i][i] += k
                if other in free: L[i][free.index(other)] -= k
                else: rhs[i] += k*pins[other]
    p = dict(zip(free,solve(L,rhs))); p.update(pins)
    cols = []
    for x,y,_ in edges:
        force = [Q(0)]*3; drop = p[x]-p[y]
        if x in free: force[free.index(x)] -= drop
        if y in free: force[free.index(y)] += drop
        z = solve(L,force); cols.append((z[1],z[2]))
    return (p['u'],p['v']),cols


def main():
    (ROOT/'RESULTS.json').unlink(missing_ok=True)
    # Exact bad hub: sh=gh=uh=vh=3/2, uv=1/8. Outputs depend only on sh/(sh+gh).
    c = [Q(3,2)]*4+[Q(1,8)]
    base_columns = [(Q(1,6),Q(1,6)),(-Q(1,6),-Q(1,6)),(Q(0),Q(0)),
                    (Q(0),Q(0)),(Q(0),Q(0))]
    base_gram = gram(base_columns)
    N = 7
    parts = [[k/N]*N for k in c]
    endpoints = [('s','h'),('g','h'),('u','h'),('v','h'),('u','v')]
    q,derived = physical([(x,y,k) for (x,y),k in zip(endpoints,c)])
    verify('base_Kirchhoff_outputs',q == (Q(1,2),Q(1,2)))
    verify('base_columns_from_node_balances',derived == base_columns)
    split_q,split_columns = physical([(x,y,k) for (x,y),arr in zip(endpoints,parts) for k in arr])
    verify('split_Kirchhoff_outputs',split_q == q)
    verify('split_columns_from_node_balances',split_columns == [col for col in derived for _ in range(N)])
    for i,k in enumerate(c):
        verify('pair_sum_'+str(i),sum(parts[i]) == k)
    Kbad2=sum(k*k for k in c)
    Knew2=sum(k*k for arr in parts for k in arr)
    verify('norm_divides_by_N',Knew2 == Kbad2/N)
    # All star trace entries computed from totals, independent of a port solver.
    total=sum(c[:4]); total_new=sum(sum(arr) for arr in parts[:4])
    for i in range(4):
        for j in range(i+1,4):
            extra=c[4] if (i,j)==(2,3) else Q(0)
            extra_new=sum(parts[4]) if (i,j)==(2,3) else Q(0)
            old=c[i]*c[j]/total+extra
            new=sum(parts[i])*sum(parts[j])/total_new+extra_new
            verify('trace_'+str(i)+'_'+str(j),old == new)
    repeated=[col for col in base_columns for _ in range(N)]
    new_gram=gram(repeated)
    for i in range(2):
        for j in range(2):
            verify('Gram_'+str(i)+'_'+str(j),new_gram[i][j] == N*base_gram[i][j])
    verify('rank_one_achieved',new_gram[0][0]>0 and
           new_gram[0][0]*new_gram[1][1]-new_gram[0][1]*new_gram[1][0] == 0)
    K0=Q(2)
    verify('selected_N_sufficient',Knew2 <= K0*K0)
    # Store d²: no rounded square root is used as an equality certificate.
    dangling2=K0*K0-Knew2
    verify('dangling_nonnegative',dangling2>=0)
    verify('exact_target_norm',Knew2+dangling2 == K0*K0)
    verify('dangling_Gram_unchanged',gram(repeated+[(Q(0),Q(0))]) == new_gram)
    out={'status':'PASS','assertions':len(checks),'checks':checks,'N':N,
         'Kbad2':str(Kbad2),'Knew2':str(Knew2),'K0':str(K0),'dangling2':str(dangling2),
         'norm_equality_checked_as_exact_square':True,'old_359_rerun':False,
         'random_networks':0,'trajectories':0,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'scope':'new parallel splitting only; authored regressions, independent audit pending'}
    (ROOT/'RESULTS.json').write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ('status','assertions','old_359_rerun')}))


if __name__=='__main__':
    main()

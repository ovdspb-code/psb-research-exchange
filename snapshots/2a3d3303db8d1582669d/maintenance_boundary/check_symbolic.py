"""Exact algebra only: no trajectories, parameter search or new cycle experiment."""
import hashlib
import json
from pathlib import Path
import sympy as s

ROOT=Path(__file__).resolve().parent
k=s.symbols('k0:5',positive=True)
a,x,h,d,K=s.symbols('a x h d K',positive=True)
L=s.Matrix([[k[0]+k[1]+k[4],-k[4]],[-k[4],k[2]+k[3]+k[4]]])
q=L.inv()*s.Matrix([k[0],k[2]])
J=q.jacobian(k)
checks={}

def check(name,expr):
    elements=list(expr) if isinstance(expr,s.MatrixBase) else [expr]
    okay=all(s.cancel(e)==0 for e in elements)
    checks[name]=okay
    if not okay:raise AssertionError(name)

check('Kirchhoff',L*q-s.Matrix([k[0],k[2]]))
check('homogeneity_Jk',J*s.Matrix(k))
base=dict(zip(k,[a,a,a,a,h]))
check('fixed_current_response',q.subs(base)-s.ones(2,1)/2)
for mode in ['common','differential']:
    if mode=='common':
        state=[a+x,a-x,a+x,a-x,h]
        target=s.Matrix([s.Rational(1,2)+d]*2)
        A=a;c=1
    else:
        state=[a+x,a-x,a-x,a+x,h]
        target=s.Matrix([s.Rational(1,2)+d,s.Rational(1,2)-d])
        A=a+h;c=5
    F=-(J.T*(q-target))
    restricted=F.subs(dict(zip(k,state)),simultaneous=True)
    E=d-x/(2*A)
    da=-E*x/(4*A**2);dx=E/(4*A)
    dh=0 if mode=='common' else -E*x/A**2
    expected=s.Matrix([da+dx,da-dx,da+dx,da-dx,dh] if mode=='common'
                      else [da+dx,da-dx,da-dx,da+dx,dh])
    check(mode+'_full_physical_gradient',restricted-expected)
    check(mode+'_norm_conservation',s.Matrix(state).dot(expected))
    dA=da+dh if mode=='differential' else da
    check(mode+'_first_integral',2*A*dA+2*c*x*dx)
    check(mode+'_ratio_dynamics',(dx*A-x*dA)/A**2-E*(1+c*(x/A)**2)/(4*A**2))
phi=s.symbols('phi',real=True)
curve=s.Matrix([K*s.cos(phi)/2]*4+[K*s.sin(phi)])
check('actuator_resource_sphere',s.trigsimp(curve.dot(curve)-K**2))
check('actuator_constant_metric_speed',s.trigsimp(curve.diff(phi).dot(curve.diff(phi))-K**2))
check('healthy_A_squared',(5*K/s.sqrt(65))**2-s.Rational(5,13)*K**2)
check('damaged_A_squared',(s.sqrt(5)*K/2)**2-s.Rational(5,4)*K**2)
t=s.sqrt(s.Rational(3,2))-1
check('repair_threshold',(s.Rational(1,2)+t)**2/(1+t*t)-s.Rational(1,2))
out={'status':'EXACT_IDENTITIES_PASS_NOT_INDEPENDENT_AUDIT','sympy':s.__version__,
     'checks':checks,'trajectories_executed':0,'random_states':0,
     'healthy_worst_time_over_B':'10/13','damaged_worst_time_over_B':'5/2',
     'scope':'Predefined four tasks; fixed conductance norm; one declared coordinated actuator; exogenous damage',
     'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(ROOT/'RESULTS.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))

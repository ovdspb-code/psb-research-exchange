#!/usr/bin/env python3
"""Deterministic exact regression checks; not a universal-proof substitute."""
import hashlib
import json
import platform
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
if __name__ == "__main__":
    (ROOT/"RESULTS.json").unlink(missing_ok=True)

import sympy as s

COUNT = 0
CHECKS = []


def check(name, condition):
    global COUNT
    if not bool(condition):
        raise AssertionError(name)
    COUNT += 1
    CHECKS.append(name)


def zero(m):
    return all(s.simplify(x) == 0 for x in m)


def matrix_psd(m):
    from itertools import combinations
    for size in range(1, m.rows + 1):
        for ix in combinations(range(m.rows), size):
            if not bool(s.simplify(m.extract(ix, ix).det()) >= 0):
                return False
    return True


def network(edges, outputs=("u", "v")):
    vertices = ["s", "g", *outputs]
    vertices += sorted({a for e in edges for a in e[:2]} - set(vertices))
    ix = {a: i for i, a in enumerate(vertices)}
    L = s.zeros(len(vertices))
    incidence = []
    for a, b, k in edges:
        d = s.zeros(len(vertices), 1)
        d[ix[a]], d[ix[b]] = 1, -1
        L += k * d * d.T
        incidence.append(d)
    F = list(range(2, len(vertices)))
    P = list(range(2 + len(outputs)))
    hidden = list(range(len(P), len(vertices)))
    Linv = L.extract(F, F).inv()
    p = s.Matrix([1, 0, *(-Linv * L.extract(F, [0]))])
    S = s.zeros(len(outputs), len(F))
    for j in range(len(outputs)):
        S[j, j] = 1
    columns = [-S * Linv * d.extract(F, [0]) * (d.T*p)[0]
               for d in incidence]
    Jall = s.Matrix.hstack(*columns)
    active = [j for j, e in enumerate(edges) if e[2] > 0]
    J = Jall[:, active]
    Gamma = L.extract(P, P)
    if hidden:
        Gamma -= (L.extract(P, hidden) * L.extract(hidden, hidden).inv()
                  * L.extract(hidden, P))
    return {"vertices": vertices, "ix": ix, "L": L, "Linv": Linv,
            "p": p, "J": J, "Jall": Jall, "Gamma": Gamma,
            "H": Gamma[2:, 2:], "K2": sum(e[2]**2 for e in edges),
            "edges": edges, "incidence": incidence, "outputs": outputs}


def energy(n, v, w):
    return s.simplify(sum(k*(v[n["ix"][a]]-v[n["ix"][b]])
                            *(w[n["ix"][a]]-w[n["ix"][b]])
                          for a, b, k in n["edges"]))


def extension(n, boundary):
    nh = len(n["vertices"]) - len(boundary)
    if not nh:
        return s.Matrix(boundary)
    hp = list(range(len(boundary)))
    hh = list(range(len(boundary), len(n["vertices"])))
    x = -n["L"].extract(hh, hh).inv()*n["L"].extract(hh, hp)*s.Matrix(boundary)
    return s.Matrix([*boundary, *x])


def parts(v):
    return s.Matrix([max(0, x) for x in v]), s.Matrix([max(0, -x) for x in v])


def clipping(n, xi, label):
    t = n["p"][2]
    check(label+":equal_outputs", all(n["p"][2+j] == t
                                     for j in range(len(n["outputs"]))))
    zeta = extension(n, [0, 0, *xi])
    z = n["H"]*s.Matrix(xi)
    drops = s.Matrix([(d.T*n["p"])[0]*(d.T*zeta)[0]
                      for d, e in zip(n["incidence"], n["edges"]) if e[2] > 0])
    check(label+":adjoint", zero(n["J"].T*z+drops))
    U, V = parts(n["p"]-s.ones(len(n["p"]), 1)*t)
    Z, W = parts(zeta)
    nb = 2 + len(n["outputs"])
    hp = U-extension(n, list(U[:nb]))
    hp2 = V-extension(n, list(V[:nb]))
    hz = Z-extension(n, list(Z[:nb]))
    hz2 = W-extension(n, list(W[:nb]))
    check(label+":hp_same", zero(hp-hp2))
    check(label+":hz_same", zero(hz-hz2))
    ep, ez, ec = energy(n, hp, hp), energy(n, hz, hz), energy(n, hp, hz)
    g = -n["Gamma"][0, 1]
    a = t*(1-t)
    z_bound = sum(-n["Gamma"][2+i, 2+j]*max(xi[i], 0)*max(-xi[j], 0)
                  for i in range(len(xi)) for j in range(len(xi)))
    check(label+":hp_energy", 0 <= ep <= a*g)
    check(label+":hz_energy", 0 <= ez <= z_bound)
    check(label+":CS", ec**2 <= ep*ez)
    norm2 = sum(x*x for x in drops)
    for nm, vv, ww in (("UZ", U, Z), ("VW", V, W)):
        val = energy(n, vv, ww)
        check(label+":"+nm+"_physical_norm", val**2 <= n["K2"]*norm2)
        boundary_pair = (s.Matrix(vv[:nb]).T*n["Gamma"]*s.Matrix(ww[:nb]))[0]
        check(label+":"+nm+"_decomposition", s.simplify(val-boundary_pair-ec) == 0)


def two_bound(n, label):
    Gamma, H, t = n["Gamma"], n["H"], n["p"][2]
    A, B = -Gamma[0,2]-Gamma[1,2], -Gamma[0,3]-Gamma[1,3]
    a, g, f = t*(1-t), -Gamma[0,1], -Gamma[2,3]
    Delta = a*A*B-g*f
    r = s.sqrt(g*f/(a*A*B))
    upper = max(H.eigenvals())
    phi = s.simplify(Delta**2/(2*n["K2"]*max(A,B)**2*(1+r)**2*upper**2))
    check(label+":Delta_positive", Delta > 0)
    check(label+":bound_PSD", matrix_psd(n["J"]*n["J"].T-phi*s.eye(2)))
    return {"Delta": str(Delta), "K2": str(n["K2"]), "Phi": str(phi),
            "lambda_actual": str(min((n["J"]*n["J"].T).eigenvals()))}


def direct_derivatives(n, label):
    evar = s.Symbol("eps", real=True)
    for j, (a, b, k) in enumerate(n["edges"]):
        modified = list(n["edges"])
        modified[j] = (a, b, k+evar)
        # Solve all free vertices directly; independent of the adjoint expression.
        ix, L = n["ix"], s.zeros(len(n["vertices"]))
        for aa, bb, kk in modified:
            d = s.zeros(L.rows, 1)
            d[ix[aa]], d[ix[bb]] = 1, -1
            L += kk*d*d.T
        ff = list(range(2, L.rows))
        qq = -L.extract(ff, ff).inv()*L.extract(ff, [0])
        deriv = s.Matrix([s.diff(x, evar).subs(evar, 0) for x in qq[:len(n["outputs"])]])
        check(label+f":direct_derivative_{j}", zero(deriv-n["Jall"][:, j]))


def run():
    Q = s.Rational
    t, A, B, zz, f = s.symbols("t A B zz f", positive=True)
    R = zz+A+B
    weights = [t*R, (1-t)*R, A*R/zz, B*R/zz]
    total = sum(weights)
    expected = {(0,1):t*(1-t)*zz, (0,2):t*A, (1,2):(1-t)*A,
                (0,3):t*B, (1,3):(1-t)*B, (2,3):A*B/zz}
    for (i,j), target in expected.items():
        check(f"reverse_symbolic_trace_{i}_{j}", s.simplify(weights[i]*weights[j]/total-target) == 0)
    check("reverse_symbolic_norm_not_fixed", s.simplify(total-R**2/zz) == 0)

    # Delta<0, one hidden hub; exactly the proposed construction.
    reverse = network([("s","h",Q(3,2)), ("g","h",Q(3,2)),
                       ("u","h",Q(3,2)), ("v","h",Q(3,2)),
                       ("u","v",Q(1,8))])
    check("reverse_rank_one", reverse["J"].rank() == 1)
    check("reverse_outputs_half", reverse["p"][2:4] == [Q(1,2)]*2)
    for xi in ([1,-1], [1,0], [0,1], [2,-1], [1,1]):
        clipping(reverse, xi, "reverse_"+str(xi))

    branch = network([("s","u",Q(2)), ("u","g",Q(3)), ("u","v",Q(7))])
    expected_J = s.Matrix([[Q(3,25),-Q(2,25),0]]*2)
    check("branch_J_exact", zero(branch["J"]-expected_J))
    check("branch_A_or_B_zero", branch["Gamma"][0,3] == branch["Gamma"][1,3] == 0)
    clipping(branch, [0,1], "branch_zero_adjoint_output")

    port_edges = [("s","u",Q(1,2)), ("g","u",Q(1,2)),
                  ("s","v",Q(1,2)), ("g","v",Q(1,2)),
                  ("s","g",Q(1,5)), ("u","v",Q(1,5))]
    direct = network(port_edges)
    # Exact hidden realization of the same Gamma: star contributes 1/20 per pair.
    hidden_edges = [(a,b,k-Q(1,20)) for a,b,k in port_edges]
    hidden_edges += [(a,"h",Q(1,5)) for a in ("s","g","u","v")]
    hidden = network(hidden_edges)
    check("same_ports_hidden_realization", zero(direct["Gamma"]-hidden["Gamma"]))
    values = {"direct": two_bound(direct,"direct"), "hidden": two_bound(hidden,"hidden")}
    for name, n in (("direct",direct),("hidden",hidden)):
        direct_derivatives(n,name)
        for xi in ([1,-1], [1,0], [0,1], [2,-1], [-1,2], [1,1], [-2,-1]):
            clipping(n,xi,name+"_"+str(xi))
    active_face = network([(a,b,k) for a,b,k in port_edges if (a,b)!=("s","g")]
                          + [("s","g",Q(0))])
    check("zero_edge_excluded", active_face["J"].cols == 5)
    values["active_face"] = two_bound(active_face,"active_face")

    series_edges = []
    for j,(aa,bb,kk) in enumerate(port_edges):
        series_edges += [(aa,"x"+str(j),2*kk),("x"+str(j),bb,2*kk)]
    series = network(series_edges)
    check("series_same_Gamma",zero(series["Gamma"]-direct["Gamma"]))
    check("series_Gram_Nminus3",zero(series["J"]*series["J"].T-direct["J"]*direct["J"].T/8))
    check("series_norm_N3",series["K2"] == 8*direct["K2"])
    for j in range(len(port_edges)):
        check("series_derivative_"+str(j),zero(series["J"][:,2*j]-direct["J"][:,j]/4)
              and zero(series["J"][:,2*j+1]-direct["J"][:,j]/4))
    values["series_N2"] = two_bound(series,"series_N2")

    for scale in (Q(1,10**14), Q(2), Q(1000)):
        nn = network([(a,b,scale*k) for a,b,k in port_edges])
        check("scale_trace_"+str(scale), zero(nn["Gamma"]-scale*direct["Gamma"]))
        check("scale_J_"+str(scale), zero(nn["J"]-direct["J"]/scale))
        check("scale_K2_"+str(scale), nn["K2"] == direct["K2"]*scale**2)
        got = two_bound(nn,"scale_"+str(scale))
        check("scale_Phi_"+str(scale), s.sympify(got["Phi"])*scale**2 == s.sympify(values["direct"]["Phi"]))

    triple_edges = [(pin,o,Q(1,2)) for pin in ("s","g") for o in ("u","v","w")]
    triple_edges += [("s","g",Q(1,5)), ("u","v",Q(1,5)),
                     ("u","w",Q(1,5)), ("v","w",Q(1,5))]
    triple = network(triple_edges, ("u","v","w"))
    phi3 = s.simplify((Q(1,4)*(1-Q(2,5)))**2/(2*triple["K2"]*max(triple["H"].eigenvals())**2))
    check("three_outputs_full_rank", triple["J"].rank() == 3)
    check("three_outputs_bound", matrix_psd(triple["J"]*triple["J"].T-phi3*s.eye(3)))
    for xi in ([1,-1,0], [1,1,-2], [0,0,1], [-1,2,-3]):
        clipping(triple,xi,"triple_"+str(xi))
    values["three_outputs"] = {"Phi": str(phi3), "rank": triple["J"].rank()}

    # Noncoincident, same fixed example: check all boundary coefficients in section 5.
    perturbed = list(port_edges)
    perturbed[0] = ("s","u",Q(51,100))
    nn = network(perturbed)
    qu,qv = nn["p"][2:4]
    m, dq = (qu+qv)/2, (qu-qv)/2
    Gamma,H = nn["Gamma"],nn["H"]
    c = lambda i,j:-Gamma[i,j]
    Cp = m*(1-m)*c(0,1)+dq*((1-m)*c(0,3)+m*c(1,2))+dq**2*c(2,3)
    Ps = [(1-m)*c(0,2)-dq*H[0,0], (1-m)*c(0,3)+dq*c(2,3)]
    Qs = [m*c(1,3)-dq*H[1,1],m*c(1,2)+dq*c(2,3)]
    check("noncoincident_criterion", all(P>0 and Qv>0 and P*Qv>Cp*c(2,3) for P,Qv in zip(Ps,Qs)))
    up, vp = s.Matrix([1-m,0,dq,0]),s.Matrix([0,m,0,dq])
    check("noncoincident_Cp_boundary", (up.T*Gamma*vp)[0] == -Cp)
    for j in (0,1):
        basis = s.zeros(4,1); basis[2+j] = 1
        check("noncoincident_P_boundary_"+str(j),(up.T*Gamma*basis)[0] == -Ps[j])
        check("noncoincident_Q_boundary_"+str(j),(vp.T*Gamma*basis)[0] == -Qs[1-j])
    values["noncoincident"] = {"q": [str(qu),str(qv)],"Cp":str(Cp),"P":list(map(str,Ps)),"Q":list(map(str,Qs))}

    # CL: exact radial pairing, derivative of the port weight, and cubic identity.
    for label,n in (("direct",direct),("hidden",hidden),("active_face",active_face),
                    ("noncoincident",nn)):
        err = s.Matrix([Q(1,1000),-Q(1,2000)])
        y = extension(n,[0,0,*err])
        raw = s.Matrix([(d.T*n["p"])[0]*(d.T*y)[0] for d in n["incidence"]])
        velocity = s.Matrix([raw[j] if k>0 else max(raw[j],0)
                             for j,(_,_,k) in enumerate(n["edges"])])
        check(label+":CL_radial_norm", sum(k*velocity[j] for j,(_,_,k) in enumerate(n["edges"])) == 0)
        check(label+":CL_sign_adjoint", zero(raw+n["Jall"].T*n["H"]*err))
        check(label+":CL_projection_pairing", (err.T*n["H"]*n["Jall"]*velocity)[0] == -(velocity.T*velocity)[0])
        port_extensions = s.Matrix.hstack(*[extension(n,[1 if i==j else 0 for i in range(4)])
                                           for j in range(4)])
        Hdot = s.zeros(2)
        for vel,d in zip(velocity,n["incidence"]):
            row = d.T*port_extensions[:,2:]
            Hdot += vel*row.T*row
        remainder = (err.T*Hdot*err)[0]/2
        quadratic = (velocity.T*velocity)[0]
        M = len(n["edges"])
        # Square the absolute bound so every decision is rational and exact.
        check(label+":CL_cubic_bound", remainder**2 <= M*quadratic*(err.T*err)[0]**2)
        check(label+":CL_weighted_decrease_small_error", -quadratic+remainder < 0)
        values[label]["CL_cubic_remainder"] = str(remainder)
        if label == "noncoincident":
            check(label+":CL_cubic_nonzero", remainder != 0)

    # Closed port-radius formula, evaluated for the fixed direct fixture.
    n = direct; H,Gamma = n["H"],n["Gamma"]
    h,U,G = min(H.eigenvals()),max(H.eigenvals()),max(Gamma.eigenvals())
    t = Q(1,2); a=t*(1-t); A=B=Q(1); g=f=Q(1,5)
    pmin=pmax=a; gap=a*(a*A*B-g*f); cp0=a*g
    D=2*(1+s.sqrt(2))/h; Lp=1+3*D*G; Lc=1+4*D*G
    Lgap=3*pmax*Lp+Q(3,2)*G*Lc+cp0
    candidates=[h/2,min(t,1-t)/(2*D),1/(2*D),pmin/(2*Lp),gap/(2*Lgap)]
    eta=min(candidates,key=lambda x:float(x))
    check("bootstrap_eta_true_min",all(s.simplify(x-eta)>=0 for x in candidates))
    ell=pmin*gap/(18*s.sqrt(2)*pmax**2)
    C=ell/(Q(3,2)*s.sqrt(n["K2"])*U)
    hs=h/2; M=len(n["edges"])
    admission_ceiling=min(C*hs/(2*s.sqrt(M)),eta*C/(4*s.sqrt(M)),key=lambda x:float(x))
    check("bootstrap_positive_radius",eta>0 and C>0 and admission_ceiling>0)
    values["bootstrap_fixed_example"]={"eta_exact":str(s.simplify(eta)),"C_exact":str(s.simplify(C)),
        "eta_approx":float(eta),"C_approx":float(C),
        "E0_strict_ceiling_approx":float(admission_ceiling),
        "time_prefactor_approx":float(2/(C*C*hs)),
        "interpretation":"conservative analytic radius; not a simulated or practical speed estimate"}

    report = {"status":"PASS", "assertions":COUNT, "checks":CHECKS,
              "fixtures":values, "python":platform.python_version(),
              "python_executable":sys.executable,
              "sympy":s.__version__,"random_networks":0,"ODE_trajectories":0,
              "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "scope":"exact deterministic regressions; universal theorem is in PROOF_RU.md; no independent audit"}
    (ROOT/"RESULTS.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("status","assertions","python","sympy","random_networks","ODE_trajectories")},ensure_ascii=False))


if __name__ == "__main__":
    run()

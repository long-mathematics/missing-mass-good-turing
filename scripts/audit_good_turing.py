#!/usr/bin/env python3
"""Audits canonical Poissonized risk structure; does not certify phase diagrams."""
import json,math,random,time
from pathlib import Path
import numpy as np
import mpmath as mp
import sympy as sy
from scipy.optimize import minimize,minimize_scalar
ROOT=Path(__file__).resolve().parents[1];rng=random.Random(20261126);mp.mp.dps=60
z,T=sy.symbols('z T',positive=True)
phi=z**T*(z+T)*sy.exp(-z)
poly=z**3-(T+2)*z**2-T*(T-1)*z+T*T*(T-1)
assert sy.simplify(sy.diff(phi,z)-z**(T-1)*sy.exp(-z)*(T*T+z-z*z))==0
assert sy.simplify(sy.diff(phi,z,2)-z**(T-2)*sy.exp(-z)*poly)==0
for theta in range(2,51):assert sy.Poly(poly.subs(T,theta),z).count_roots(0,sy.oo)==2
cases=0;worst=mp.mpf(0)
for r in range(9):
 for t in map(mp.mpf,['.05',1,3,17,100]):
  for x in map(mp.mpf,['.0001','.01','.17','.5',1]):
   lam=t*x
   pr=mp.exp(-lam)*lam**r/mp.factorial(r)
   pr1=mp.exp(-lam)*lam**(r+1)/mp.factorial(r+1)
   bias=(r+1)/t*pr1-x*pr
   assert abs(bias)<mp.mpf('1e-55')
   direct=x*x*pr+((r+1)/t)**2*pr1
   reduced=lam**(r+1)*(lam+r+1)*mp.exp(-lam)/(mp.factorial(r)*t*t)
   rel=abs(direct-reduced)/max(1,abs(direct));worst=max(worst,rel);assert rel<mp.mpf('1e-55');cases+=1

def ph(x,theta):
 x=np.asarray(x)
 return x**theta*(x+theta)*np.exp(-x)
runs=0;max_adv=-1e99
for n in [2,3,5]:
 for r in [0,1,2,5]:
  theta=r+1
  for t in [.5,3,10,30]:
   best=float(ph(t,theta))
   for m in range(2,n+1):
    grid=np.linspace(0,1/(m-1),1601)
    val=(m-1)*ph(t*grid,theta)+ph(t*(1-(m-1)*grid),theta)
    best=max(best,float(max(val)))
    for i in range(1,len(grid)-1):
     if val[i]>=val[i-1] and val[i]>=val[i+1]:
      rr=minimize_scalar(lambda a:-float((m-1)*ph(t*a,theta)+ph(t*(1-(m-1)*a),theta)),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-13})
      best=max(best,-rr.fun)
   for _ in range(6):
    w=np.array([rng.expovariate(1) for __ in range(n)]);start=w/w.sum()
    rr=minimize(lambda p:-float(np.sum(ph(t*p,theta)))/max(1,best),start,method='SLSQP',bounds=[(0,1)]*n,constraints={'type':'eq','fun':lambda p:sum(p)-1},options={'ftol':1e-12,'maxiter':400})
    assert abs(sum(rr.x)-1)<1e-7
    adv=(-rr.fun*max(1,best)-best)/max(1,best);max_adv=max(max_adv,adv)
    assert adv<1e-8,(n,r,t,adv,rr.x)
    runs+=1
res={'seed':20261126,'symbolic_derivative_identities':2,'exact_positive_root_counts_theta':list(range(2,51)),'poisson_coordinate_risk_cases':cases,'max_relative_risk_error':str(worst),'full_simplex_optimizer_runs':runs,'max_relative_advantage_over_one_exceptional_family':max_adv,'numeric_search_is_not_proof':True,'scope':'Canonical unbiased Poissonized risk structure only; no fixed-count or complete phase certification.'}
(ROOT/'checks'/'good-turing-results.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))

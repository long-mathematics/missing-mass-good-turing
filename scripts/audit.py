#!/usr/bin/env python3
"""Independent exact and numerical checks. Numerical searches are not proofs."""
from __future__ import annotations
import argparse,itertools,json,math,random,time
from fractions import Fraction as F
from pathlib import Path
import mpmath as mp
import numpy as np
import scipy
from scipy.optimize import brentq,minimize,minimize_scalar
import sympy as sy

ROOT=Path(__file__).resolve().parents[1]
mp.mp.dps=70
rng=random.Random(20261126)
results={"seed":20261126,"precision_digits":70,"versions":{"sympy":sy.__version__,"numpy":np.__version__,"scipy":scipy.__version__,"mpmath":mp.__version__}}

def exact_checks():
    z,l,s=sy.symbols('z l s',positive=True)
    A=sy.exp(-(l-s)*z);B=sy.exp(-l*z);g=1+A-B
    identities={
      'mgf_curvature':sy.diff(g,z,2)*g-sy.diff(g,z)**2-((l-s)**2*A-l*l*B-s*s*A*B),
      'lambda_derivative':sy.diff(g,l)+z*(g-1),
      'affine_duality_exponentiated':g-sy.exp((s-l)*z)*(1+sy.exp(-(s-l)*z)-sy.exp(-s*z)),
    }
    for name,expr in identities.items(): assert sy.simplify(expr)==0,name
    a,p,q,P,Q=sy.symbols('a p q P Q')
    count=0
    for k in range(1,11):
        n=k+1;D=lambda j:1+P*p**j*(1-Q*q**j)
        gam=sum(p**(-j) for j in range(1,k+1));alpha=k*a/p;beta=(a+n)*p**k;S=alpha+beta-n*(a+1)
        E0=D(1)-q*D(0)
        H=gam*D(n)*E0+q*D(0)*(q**k*D(1)-D(n))
        N=alpha*D(1)*D(n)+beta*q**n*D(0)*D(1)-n*(a+1)*q*D(0)*D(n)
        assert sy.cancel((1-p)*(alpha-beta*gam)-S)==0
        assert sy.cancel(N-S*(q*D(0)*D(n)+D(n)*E0/(1-p))-beta*H)==0
        # Independent telescoping identity for the induction, treating D_j as free.
        j=k;Dj,Djp,D0,D1,E=sy.symbols('Dj Djp D0 D1 E')
        lhs=p**(-j)*E+q*D0*D1*(q**j/Djp-q**(j-1)/Dj)
        rhs=p**(-j)*E-q**j*D0*D1/(Dj*Djp)*(Djp-q*Dj)
        assert sy.cancel(lhs-rhs)==0
        count+=3
    # Derivative identity in spinodal proof, reduced modulo the root relation.
    d,t=sy.symbols('d t',positive=True)
    FF=sy.exp(d)*(d*d*sy.exp(t)-(d+t)**2)-t*t
    expr=sy.diff(FF,t)-sy.diff(FF,d)-(d*sy.exp(d+t)*(d-2)-t*t-2*t)
    assert sy.simplify(expr+FF)==0
    results['symbolic']={'identities_passed':len(identities)+count+1,'appendix_k':list(range(1,11))}

def comparison_stress():
    ntest=0;implies=0;maxrel=mp.mpf('0');minH=mp.inf
    for _ in range(900):
        k=rng.choice([1,2,3,5,10,30]);n=k+1
        a=mp.mpf(10)**mp.mpf(rng.uniform(-4,3))
        M=mp.mpf(10)**mp.mpf(rng.uniform(-5,1));R=mp.mpf(10)**mp.mpf(rng.uniform(-5,1))
        p=mp.exp(-M);q=mp.exp(-R)
        D=lambda j:1+mp.exp(-M*(a+j))*(-mp.expm1(-R*(a+j)))
        E=lambda j:D(j+1)-q*D(j)
        j=rng.choice([1,2,5,10,30])
        assert E(0)>0 and E(j)>0
        assert p**j*E(j)<E(0)
        assert D(j)*D(j+1)>q**j*D(0)*D(1)
        gam=sum(p**(-r) for r in range(1,k+1));S=k*a/p+(a+n)*p**k-n*(a+1)
        H=gam*D(n)*E(0)+q*D(0)*(q**k*D(1)-D(n))
        assert H>0
        N=k*a/p*D(1)*D(n)+(a+n)*p**k*q**n*D(0)*D(1)-n*(a+1)*q*D(0)*D(n)
        recon=S*(q*D(0)*D(n)+D(n)*E(0)/(1-p))+(a+n)*p**k*H
        maxrel=max(maxrel,abs(N-recon)/max(1,abs(N),abs(recon)))
        if S>=0: assert N>0;implies+=1
        minH=min(minH,H);ntest+=1
    results['comparison_numeric']={'cases':ntest,'hypothesis_satisfied':implies,'max_relative_identity_residual':str(maxrel),'min_remainder_H':str(minH),'not_a_proof':True}

def compositions(t,n):
    if n==1:yield (t,);return
    for i in range(t+1):
        for tail in compositions(t-i,n-1):yield (i,)+tail

def fixed_mgf(p,t,s):
    ans=mp.mpf('0')
    for counts in compositions(t,len(p)):
        prob=mp.factorial(t)
        for pi,c in zip(p,counts):prob*=pi**c/mp.factorial(c)
        unseen=sum(pi for pi,c in zip(p,counts) if c==0)
        ans+=prob*mp.exp(s*unseen)
    return ans

def egf(p,t,s):
    co=[mp.mpf(1)]+[mp.mpf(0)]*t
    for pi in p:
        f=[mp.exp(s*pi)]+[pi**j/mp.factorial(j) for j in range(1,t+1)]
        co=[sum(co[k]*f[j-k] for k in range(j+1)) for j in range(t+1)]
    return mp.factorial(t)*co[t]

def fixed_checks():
    worst=mp.mpf(0);cases=0
    for n in range(2,6):
      for t in range(1,6):
       for s0 in [-7,-1,mp.mpf('0.2'),3,11]:
        weights=[rng.randrange(1,12) for _ in range(n)];p=[mp.mpf(w)/sum(weights) for w in weights]
        left=fixed_mgf(p,t,s0);right=egf(p,t,s0)
        err=abs(left-right)/max(1,abs(left));worst=max(worst,err);assert err<mp.mpf('1e-60')
        # Independently compute frozen derivative via log-series coefficients.
        def ell(m,x):
          return x**m/mp.factorial(m)*sum((-1)**(r+1)*mp.factorial(r-1)*int(sy.functions.combinatorial.numbers.stirling(m,r))*mp.exp(-r*s0*x) for r in range(1,m+1))
        L=[mp.mpf(0)]+[sum(ell(m,x) for x in p) for m in range(1,t+1)]
        B=[mp.mpf(1)]
        for j in range(1,t+1): B.append(sum(m*L[m]*B[j-m] for m in range(1,j+1))/j)
        assert min(B)>0
        assert abs(mp.exp(s0)*mp.factorial(t)*B[t]-left)/max(1,abs(left))<mp.mpf('1e-58')
        phi=lambda x:sum(B[t-m]/B[t]*mp.diff(lambda xx:ell(m,xx),x) for m in range(1,t+1))
        direct=mp.diff(lambda e:mp.log(fixed_mgf([p[0]+e,p[1]-e]+p[2:],t,s0)),mp.mpf(0))
        assert abs(direct-(phi(p[0])-phi(p[1])))<mp.mpf('1e-53')
        cases+=1
    results['fixed_sample_mgf']={'cases':cases,'checks':['multinomial_vs_EGF','truncated_log_series','coefficient_positivity','frozen_directional_derivative'],'max_relative_mgf_error':str(worst)}

def tail_direct(p,t,a):
    ans=F(0)
    for counts in compositions(t,len(p)):
        prob=F(math.factorial(t))
        for pi,c in zip(p,counts):prob*=pi**c/F(math.factorial(c))
        if sum(pi for pi,c in zip(p,counts) if c==0)>=a:ans+=prob
    return ans

def tail_ie(p,t,a):
    total=F(0);n=len(p)
    for m in range(1,min(n,t)+1):
     for S in itertools.combinations(range(n),m):
      if sum(p[i] for i in S)>1-a:continue
      prob=F(0)
      for r in range(m+1):
       for R in itertools.combinations(S,r):prob+=(-1)**(m-r)*sum((p[i] for i in R),F(0))**t
      total+=prob
    return total

def tail_checks():
    cases=0
    for n in range(2,5):
     for ints in compositions(7,n):
      p=[F(i,7) for i in ints]
      for t in [1,2,3]:
       for a in [F(1,5),F(1,2),F(3,5),F(4,5)]:
        ans=tail_direct(p,t,a);assert ans==tail_ie(p,t,a)
        plateau=sum(sorted(p,reverse=True)[:min(t,n)])<=1-a
        assert (ans==1)==plateau
        if n<math.ceil(F(t)/(1-a)):assert ans<1
        cases+=1
    assert tail_direct([F(2,5),F(2,5),F(1,5)],2,F(3,5))==F(9,25)
    assert tail_direct([F(1,3)]*3,2,F(3,5))==F(1,3)
    results['fixed_tails_exact']={'rational_cases':cases,'denominator':7,'checks':['multinomial_vs_support_inclusion_exclusion','probability_one_iff_top_t_mass','nonsaturation','uniform_counterexample_9_25_vs_1_3']}

# Stable second derivative and factorized scalar excess for late lambda>=s.
quad_x,quad_w=np.polynomial.legendre.leggauss(40)
quad_r=(quad_x+1)/2;quad_w=quad_w/2*(1-quad_r)
def h(z,l,s):
    return np.log1p(np.exp(-(l-s)*np.asarray(z))*(-np.expm1(-s*np.asarray(z))))
def hpp(z,l,s):
    z=np.asarray(z);a=l-s;A=np.exp(-a*z);B=np.exp(-l*z);gg=1+A-B
    return (a*a*A-l*l*B-s*s*A*B)/(gg*gg)
def G(x,n,l,s):
    x=np.asarray(x);u=1/n;k=n-1
    zz=u+(x[...,None]-u)*quad_r
    return np.sum((k*hpp(zz,l,s)+k*k*hpp(1-k*zz,l,s))*quad_w,axis=-1)
def maxG(n,l,s):
    grid=np.linspace(0,1/n,101);vals=G(grid,n,l,s);best=float(max(vals))
    for i in range(1,len(grid)-1):
      if vals[i]>=vals[i-1] and vals[i]>=vals[i+1]:
        rr=minimize_scalar(lambda x:-float(G(x,n,l,s)),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-13})
        best=max(best,-float(rr.fun))
    return best

def sigma(n,s):
    t=s/n
    f=lambda d:d*d*math.exp(t)-(d+t)**2-t*t*math.exp(-d)
    return s+n*brentq(f,1e-15,2,xtol=1e-13)
def beta(n,s):
    hi=sigma(n,s)
    return brentq(lambda l:maxG(n,l,s),s,hi+max(1e-6,hi*1e-7),xtol=2e-10)

def phase_checks():
    rows=[];signs=0
    for n in [2,3,5,10]:
      ss=[.05,.3,1,3,8,16]
      bs=[beta(n,s) for s in ss]
      assert all(bs[i]<bs[i+1] for i in range(len(bs)-1))
      for s,b in zip(ss,bs):
        lo=sigma(1,s);hi=sigma(n,s)
        assert s<lo<=b+2e-8 and b<=hi+2e-8
        for l in np.linspace(s,b+max(1,.25*n),21):
          val=maxG(n,float(l),s)
          if l<b-1e-7: assert val<1e-8
          if l>b+1e-7: assert val>0
          signs+=1
        # Actual critical Hessian has the two claimed transverse signs.
        mp_s=mp.mpf(s);mp_l=mp.mpf(hi)
        KK=lambda l0,s0:(l0-s0)**2*mp.exp(l0/n)-l0*l0*mp.exp((l0-s0)/n)-s0*s0
        assert mp.diff(lambda z:KK(z,mp_s),mp_l)>0
        assert mp.diff(lambda z:KK(mp_l,z),mp_s)<0
        rows.append({'n':n,'s':s,'sigma_global_concavity':lo,'beta_numeric':b,'sigma_uniform_spinodal':hi})
    results['phase_exploration']={'rows':rows,'factorized_excess_sign_checks':signs,'not_a_proof':True,'method':'40-node Gauss-Legendre evaluation of quadratic factor; 101-point grid plus bounded refinement; Brent root'}

def simplex_stress():
    cases=0;max_adv=-1e99
    for n in [2,3,4,6]:
      for l,s in [(1,.3),(4,1),(8,3),(15,2),(20,12),(3,8)]:
        # Duality means optimization can be performed entirely in the late half-plane.
        ll,ss=max(l,s),min(l,s)
        xx=np.linspace(0,1/n,4001)
        vv=(n-1)*h(xx,ll,ss)+h(1-(n-1)*xx,ll,ss)
        best=float(max(vv))
        for i in range(1,len(xx)-1):
          if vv[i]>=vv[i-1] and vv[i]>=vv[i+1]:
            rr=minimize_scalar(lambda x:-float((n-1)*h(x,ll,ss)+h(1-(n-1)*x,ll,ss)),bounds=(xx[i-1],xx[i+1]),method='bounded',options={'xatol':1e-14})
            best=max(best,-rr.fun)
        starts=[np.ones(n)/n,np.eye(n)[0]]
        for _ in range(8):
            w=np.array([rng.expovariate(1) for __ in range(n)]);starts.append(w/w.sum())
        for start in starts:
          rr=minimize(lambda p:-float(np.sum(h(p,ll,ss))),start,method='SLSQP',bounds=[(0,1)]*n,constraints={'type':'eq','fun':lambda p:np.sum(p)-1},options={'ftol':1e-12,'maxiter':300})
          assert abs(sum(rr.x)-1)<1e-7
          adv=-rr.fun-best;max_adv=max(max_adv,adv)
          assert adv<2e-9,(n,l,s,adv,rr.x)
          cases+=1
    results['full_simplex_search']={'optimizer_runs':cases,'max_advantage_over_one_heavy_envelope':max_adv,'not_a_proof':True}

if __name__=='__main__':
    start=time.time()
    for func in [exact_checks,comparison_stress,fixed_checks,tail_checks,phase_checks,simplex_stress]:
        print(func.__name__,flush=True);func()
    results['elapsed_seconds']=time.time()-start
    out=ROOT/'checks'/'audit-results.json';out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(results,indent=2)+'\n')
    print('PASS',out,flush=True)

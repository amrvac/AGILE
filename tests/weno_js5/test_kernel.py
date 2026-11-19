#!/usr/bin/env python3
"""Compile production JS5 and compare with independent polynomial integrals.
Run with the repository's .venv/bin/python (numpy and fypp required).
"""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import numpy as np
from numpy.polynomial import Polynomial

ROOT = Path(__file__).resolve().parents[2]

def reference(v, epsilon):
    """Recover quadratics from cell averages, then integrate derivative squares."""
    q, beta = [], []
    for offset in range(3):
        centers = np.arange(offset-2, offset+1, dtype=float)
        matrix = np.array([[( (x+.5)**(p+1)-(x-.5)**(p+1))/(p+1)
                            for p in range(3)] for x in centers])
        polynomial = Polynomial(np.linalg.solve(matrix, v[offset:offset+3]))
        q.append(polynomial(.5))
        indicator = 0.
        for order in (1, 2):
            derivative = polynomial.deriv(order)
            integral = (derivative*derivative).integ()
            indicator += integral(.5)-integral(-.5)
        beta.append(indicator)
    beta = np.array(beta)
    alpha = np.array([.1,.6,.3])/(epsilon+beta)**2
    weights = alpha/alpha.sum()
    return np.r_[np.dot(weights,q), q, beta, weights]

DRIVER = '''program check_js
use mod_weno_js5
implicit none
double precision :: v(5), eps, val, q(3), b(3), w(3)
integer :: ios
do
  read(*,*,iostat=ios) eps,v
  if (ios/=0) exit
  call weno_js5_left(v,eps,val,q,b,w)
  write(*,'(10(es25.16,1x))') val,q,b,w
end do
end program
'''

def main():
    rng=np.random.default_rng(1996)
    cases=[]
    for eps in (1e-6,1e-18,1e-3):
        for scale in (1.,1e-10,1e4):
            cases.extend((eps,v*scale) for v in rng.normal(size=(40,5)))
        for v in ([1]*5,[-2,-1,0,1,2],[0,0,0,1,1],[1,1,1,0,0]):
            cases.append((eps,np.array(v,dtype=float)))
    x=np.arange(-2,3,dtype=float)
    # Exact cell averages of degree <= 2 polynomials.
    polys=[np.ones(5), 2+3*x, 2+3*x+4*(x*x+1/12)]
    cases.extend((1e-6,v) for v in polys)
    original=list(cases)
    cases.extend((eps,v[::-1]) for eps,v in original)
    # A smooth noncritical function and a critical-point function, exact averages.
    convergence={}
    convergence_start=len(cases)
    profiles=[('smooth',1e-6),('critical',1e-6),('critical-small-epsilon',1e-18)]
    for label,eps in profiles:
        for n in (20,40,80,160,320):
            h=1/n
            if label=='smooth':
                v=np.sin(.3+x*h)*np.sinc(h/(2*np.pi))
            else:
                # f=sin(x)+cos(2x)/2 has f'=0 at pi/6, f''' != 0.
                center=np.pi/6+x*h
                v=np.sin(center)*np.sinc(h/(2*np.pi))+.5*np.cos(2*center)*np.sinc(h/np.pi)
            cases.append((eps,v))
    with tempfile.TemporaryDirectory(prefix='agile-js5-') as tmp:
        tmp=Path(tmp)
        fypp=shutil.which('fypp') or str(ROOT/'.venv/bin/fypp')
        subprocess.run([fypp,str(ROOT/'src/limiter/mod_weno_js5.fpp'),str(tmp/'weno.f90')],check=True)
        (tmp/'driver.f90').write_text(DRIVER)
        subprocess.run([os.environ.get('FC','gfortran'),'-O0','-fcheck=all',
                        '-ffpe-trap=invalid,zero,overflow',str(tmp/'weno.f90'),
                        str(tmp/'driver.f90'),'-o',str(tmp/'check')],cwd=tmp,check=True)
        data=''.join(' '.join(map(str,[eps,*v]))+'\n' for eps,v in cases)
        result=subprocess.run([str(tmp/'check')],input=data,text=True,capture_output=True,check=True)
        actual=np.loadtxt(result.stdout.splitlines())
    expected=np.array([reference(v,eps) for eps,v in cases])
    np.testing.assert_allclose(actual[:,:7],expected[:,:7],rtol=2e-12,atol=2e-12)
    # At a critical point with tiny epsilon, roundoff in derivative differences
    # is amplified in the weights. Bound it using the indicator's amplitude,
    # rather than relaxing the candidate/beta/reconstruction checks.
    scales=np.array([max(abs(v)) for eps,v in cases])
    epsilons=np.array([eps for eps,v in cases])
    weight_tol=np.maximum(2e-12,16*np.finfo(float).eps*scales/
                          np.sqrt(epsilons+expected[:,4:7].min(axis=1)))
    assert np.all(abs(actual[:,7:]-expected[:,7:])<=weight_tol[:,None])
    assert np.isfinite(actual).all() and (actual[:,7:]>=0).all()
    np.testing.assert_allclose(actual[:,7:].sum(axis=1),1.,atol=4e-16)
    np.testing.assert_allclose(actual[len(original)-3:len(original),0],[1,3.5,4.5],atol=3e-15)
    for j,(label,eps) in enumerate(profiles):
        error=[]
        for k,n in enumerate((20,40,80,160,320)):
            xf=np.pi/6+.5/n
            exact=np.sin(.3+.5/n) if label=='smooth' else np.sin(xf)+.5*np.cos(2*xf)
            error.append(abs(actual[convergence_start+j*5+k,0]-exact))
        convergence[label]={'epsilon':eps,'errors':error,'orders':np.log2(np.array(error[:-1])/error[1:]).tolist()}
    print(json.dumps({'cases':len(cases),'status':'PASS','convergence':convergence},indent=2))

if __name__=='__main__': main()

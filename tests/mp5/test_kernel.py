#!/usr/bin/env python3
"""Independent cell-integral polynomial / interval reference for production MP5."""
import json,subprocess,tempfile,shutil
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[2]
def reference(v,tol):
 x=np.arange(-2,3,dtype=float)
 A=np.array([[((c+.5)**(k+1)-(c-.5)**(k+1))/(k+1) for k in range(5)] for c in x])
 q=np.polynomial.polynomial.polyval(.5,np.linalg.solve(A,v))
 def mm(*a):
  return min(a) if min(a)>0 else max(a) if max(a)<0 else 0.
 mp=v[2]+mm(v[3]-v[2],4*(v[2]-v[1]))
 if (q-v[2])*(q-mp)<=tol:return q
 d=np.diff(v,n=2)
 def curvature(a,b):return mm(4*a-b,4*b-a,a,b)
 md=(v[2]+v[3]-curvature(d[1],d[2]))/2
 lc=v[2]+(v[2]-v[1])/2+4*curvature(d[0],d[1])/3
 ul=v[2]+4*(v[2]-v[1])
 lo=max(min(v[2],v[3],md),min(v[2],ul,lc))
 hi=min(max(v[2],v[3],md),max(v[2],ul,lc))
 return np.clip(q,lo,hi)
def main():
 rng=np.random.default_rng(1997);cases=[]
 x=np.arange(-2,3.)
 for tol in (0.,1e-10):
  for scale in (1e-12,1.,1e4):
   for v in rng.normal(size=(300,5))*scale:cases.extend([(tol,v),(tol,v[::-1])])
  for v in (np.ones(5),x,x*x+1/12,np.zeros(5),np.array([0,0,1,1,1.])):cases.append((tol,v))
 # Explicitly straddle the early-acceptance product without relying on
 # exact floating-point equality at the branch boundary.
 for v in (np.array([0.,0.,1.,1.,0.]),np.array([0.,1.,1.,0.,0.])):
  q=(2*v[0]-13*v[1]+47*v[2]+27*v[3]-3*v[4])/60
  a,b=v[3]-v[2],4*(v[2]-v[1])
  mp=v[2]+(min(a,b) if min(a,b)>0 else max(a,b) if max(a,b)<0 else 0.)
  product=(q-v[2])*(q-mp)
  assert product>0
  cases.extend([(product*(1-1e-6),v),(product*(1+1e-6),v)])
 driver='''program test
use mod_mp5
implicit none
real(8) :: v(5),t,q
integer :: ios
do
read(*,*,iostat=ios)t,v
if(ios/=0)exit
call mp5_left(v,t,q)
write(*,'(es25.16)')q
end do
end program
'''
 with tempfile.TemporaryDirectory() as td:
  td=Path(td);(td/'test.f90').write_text(driver)
  subprocess.run([shutil.which('fypp') or str(ROOT/'.venv/bin/fypp'),str(ROOT/'src/limiter/mod_mp5.fpp'),str(td/'mp5.f90')],check=True)
  subprocess.run(['gfortran','-O0','-fcheck=all','-ffpe-trap=invalid,zero,overflow',str(td/'mp5.f90'),str(td/'test.f90'),'-o',str(td/'test')],cwd=td,check=True)
  r=subprocess.run([str(td/'test')],input=''.join(' '.join(map(str,[t,*v]))+'\n' for t,v in cases),text=True,capture_output=True,check=True)
  actual=np.loadtxt(r.stdout.splitlines())
 expected=np.array([reference(v,t) for t,v in cases]);scale=np.array([max(1.,max(abs(v))) for t,v in cases])
 assert np.all(abs(actual-expected)<2e-12*scale),np.max(abs(actual-expected)/scale)
 result={'cases':len(cases),'max_scaled_error':float(np.max(abs(actual-expected)/scale)),'status':'PASS'}
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()

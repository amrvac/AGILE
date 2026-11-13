#!/usr/bin/env python3
"""Test actual FV face reconstruction using a test-only public build of its module.
No production API is exported for tests. Requires a built HD validation case.
"""
import argparse
import json
from pathlib import Path
import subprocess
import tempfile
import numpy as np
from test_kernel import reference

ROOT=Path(__file__).resolve().parents[2]
DRIVER='''program faces
use mod_finite_volume, only: weno_reconstruct_face, highorder_face_primitive
use mod_variables
use mod_physics_vars
use mod_global_parameters, only: weno_js_epsilon
implicit none
double precision :: u(5,7), l(5,2), r(5,2), donor(5), bad(5)
integer :: ios,mode
iw_rho=1
iw_e=5
allocate(iw_mom(3))
iw_mom=[2,3,4]
gamma_1=.4d0
inv_gamma_1=2.5d0
do
 read(*,*,iostat=ios) mode,weno_js_epsilon,u
 if(ios/=0) exit
 if(mode==0) then
   call weno_reconstruct_face(u(:,1:6),l(:,1),r(:,1))
   call weno_reconstruct_face(u(:,2:7),l(:,2),r(:,2))
   write(*,'(20(es25.16,1x))') l(:,1),r(:,1),l(:,2),r(:,2)
 else
   bad=u(:,1)
   donor=u(:,2)
   call highorder_face_primitive(bad,donor)
   write(*,'(5(es25.16,1x))') bad
 end if
end do
end program
'''

def primitive(x):
    rho=1+.2*np.sin(x)
    v=.4+.1*np.cos(1.3*x)
    p=1+.1*np.sin(.7*x)
    return np.array([rho,v,np.zeros_like(x),np.zeros_like(x),p])

def conserved(x):
    q=primitive(x)
    q[4]=q[4]/.4+.5*q[0]*q[1]**2
    q[1]*=q[0]
    return q

def to_primitive(q):
    q=q.copy()
    q[1:4]/=q[0]
    q[4]=.4*(q[4]-.5*q[0]*np.sum(q[1:4]**2))
    return q

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--build',type=Path,default=ROOT/'tests/weno_js5/hd/build/latest')
    args=ap.parse_args()
    build=args.build.resolve()
    samples=[]
    exact=[]
    nodes,weights=np.polynomial.legendre.leggauss(12)
    for h in (.2,.1,.05,.025,.0125):
        centers=.37+np.arange(-3,4)*h
        u=np.stack([conserved(c+.5*h*nodes)@weights/2 for c in centers],axis=1)
        samples.append(u)
        exact.append(np.stack([primitive(.37-h/2),primitive(.37-h/2),
                               primitive(.37+h/2),primitive(.37+h/2)],axis=1))
    rng=np.random.default_rng(96)
    for i in range(30):
        u=np.zeros((5,7))
        u[0]=1+rng.random(7)
        u[1:4]=rng.normal(0,.1,(3,7))
        u[4]=3+rng.random(7)
        samples.append(u)
    with tempfile.TemporaryDirectory(prefix='agile-faces-') as tmp:
        tmp=Path(tmp)
        source=(build/'f90/mod_finite_volume.f90').read_text()
        # Only visibility is changed; the production routines are compiled intact.
        (tmp/'fv.f90').write_text(source.replace('  private\n','  public\n',1))
        (tmp/'driver.f90').write_text(DRIVER)
        flags=['mpif90','-O0','-g','-fcheck=all','-ffree-line-length-none','-cpp','-I'+str(build/'f90')]
        subprocess.run([*flags,'-c',str(tmp/'fv.f90'),'-o',str(tmp/'fv.o')],cwd=tmp,check=True)
        objects=[str(p) for p in (build/'obj').glob('*.o') if p.name not in ('agile.o','mod_finite_volume.o')]
        subprocess.run([*flags,str(tmp/'driver.f90'),str(tmp/'fv.o'),*objects,'-o',str(tmp/'check')],cwd=tmp,check=True)
        data=''.join(' '.join(map(str,[0,1e-6,*u.ravel(order='F')]))+'\n' for u in samples)
        fallback=[]
        for bad in ([0,0,0,0,2.5],[-1,0,0,0,2.5],[1,0,0,0,-1]):
            u=np.zeros((5,7));u[:,0]=bad;u[:,1]=[1,.2,0,0,2.52]
            fallback.append(u)
            data+=' '.join(map(str,[1,1e-6,*u.ravel(order='F')]))+'\n'
        result=subprocess.run([str(tmp/'check')],input=data,text=True,capture_output=True,check=True)
    lines=result.stdout.splitlines()
    errors=[]
    for i,u in enumerate(samples):
        row=np.fromstring(lines[i],sep=' ')
        actual=row.reshape((5,4),order='F')
        expected=[]
        for indices in (np.arange(5),np.arange(5,0,-1),np.arange(1,6),np.arange(6,1,-1)):
            expected.append(to_primitive(np.array([reference(u[j,indices],1e-6)[0] for j in range(5)])))
        np.testing.assert_allclose(actual,np.array(expected).T,rtol=3e-13,atol=3e-13)
        if i<len(exact): errors.append(np.max(np.abs(actual-exact[i])))
    for line in lines[len(samples):]:
        row=np.fromstring(line,sep=' ')
        np.testing.assert_allclose(row,[1,.2,0,0,1],atol=1e-15)
    orders=np.log2(np.array(errors[:-1])/errors[1:])
    assert np.all(orders[-2:]>4.5),orders
    print(json.dumps({'status':'PASS','stencils':len(samples),'fallback_cases':3,
                      'nonlinear_face_errors':errors,'orders':orders.tolist()},indent=2))
if __name__=='__main__': main()

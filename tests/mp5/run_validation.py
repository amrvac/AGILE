#!/usr/bin/env python3
"""Reproducible integration checks; outputs go outside tracked source files.
Build tests/weno_js5/hd first. --quick runs the infrastructure subset;
--convergence adds four resolutions and a halved time-step check.
"""
import sys
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import time
import numpy as np
import f90nml

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'tests/weno_js5'))
BASE=ROOT/'tests/weno_js5/hd/agile.par'

def run_case(binary, folder, *, n=32, block=16, axis=1, limiter='mp5',
             ghost=2, ranks=1, levels=1, stop=.001, dt=.0001, flux='hllc',
             epsilon=1e-6, tolerance=1e-10, reflect=False, expected_error=None, stepper='threestep', energy=True, layout='per_cell'):
    folder.mkdir(parents=True,exist_ok=True)
    par=f90nml.read(BASE)
    par['methodlist']['limiter']=[limiter]*20 if isinstance(limiter,str) else limiter
    par['methodlist']['flux_scheme']=[flux]*20
    par['methodlist']['weno_js_epsilon']=epsilon
    par['methodlist']['mp5_tolerance']=tolerance
    par['methodlist']['time_stepper']=stepper
    if layout is not None: par['methodlist']['fv_face_flux']=layout
    par['hd_list']['hd_energy']=energy
    if not energy:
        for key in par['boundlist']:
            if key.startswith('typeboundary'): par['boundlist'][key]=par['boundlist'][key][:4]
    par['boundlist']['nghostcells']=ghost
    par['stoplist']['time_max']=stop
    par['paramlist']['dtpar']=dt
    mesh=par['meshlist']
    mesh.update(iprob=axis,refine_max_level=levels,ditregrid=1)
    for d in range(1,4):
        mesh[f'domain_nx{d}']=n if d==abs(axis) or axis==4 or (axis>4 and d==1) else 8
        mesh[f'block_nx{d}']=block if d==abs(axis) or axis==4 or (axis>4 and d==1) else 8
    if reflect:
        for d in range(1,4):
            bc=['symm']*5
            bc[d]='asymm'
            par['boundlist'][f'typeboundary_min{d}']=bc
            par['boundlist'][f'typeboundary_max{d}']=bc
    if axis in (10,11):
        par['boundlist']['typeboundary_min1']=['cont']*5
        par['boundlist']['typeboundary_max1']=['cont']*5
    if axis==11:
        mesh['xprobmin1']=-5.
        mesh['xprobmax1']=5.
    par.write(folder/'agile.par',force=True)
    for p in folder.glob('cells-*.txt'): p.unlink()
    (folder/'diagnostics.txt').unlink(missing_ok=True)
    command=[str(binary),'-i','agile.par']
    if ranks>1: command=['mpirun','--oversubscribe','-np',str(ranks),*command]
    start=time.monotonic()
    with (folder/'run.log').open('w') as log:
        result=subprocess.run(command,cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=600,
                              env={**os.environ,"GFORTRAN_UNBUFFERED_ALL":"y"})
    log=(folder/'run.log').read_text()
    if expected_error:
        assert expected_error in log,(folder,log[-2000:])
        return {'rejected':expected_error}
    if result.returncode:
        raise AssertionError(f'{folder}: exit {result.returncode}\n{log[-3000:]}')
    stats=np.atleast_2d(np.loadtxt(folder/'diagnostics.txt'))
    cells=np.concatenate([np.atleast_2d(np.loadtxt(p)) for p in folder.glob('cells-*.txt')])
    cells=cells[np.lexsort((cells[:,2],cells[:,1],cells[:,0]))]
    if levels>1: assert len(np.unique(cells[:,3]))>1, (folder,'no coarse/fine interface')
    assert np.isfinite(cells).all() and np.all(cells[:,4]>0) and np.all(cells[:,-1]>0), folder
    match=re.search(r'(?:WENO-JS5|MP5) donor-state fallbacks:\s*(\d+)',log)
    fallbacks=int(match[1]) if match else None
    if axis<10 and fallbacks is not None: assert fallbacks==0,(folder,fallbacks)
    change=stats[-1,3:8]-stats[0,3:8]
    if axis>=10:
        if axis==10:
            boundary_flux=np.array([0.,.9,0.,0.,0.])
        else:
            rho,v,p=3.857143,2.629369,10.33333
            energy=p/.4+.5*rho*v*v
            boundary_flux=np.array([rho*v,rho*v*v+p-1,0,0,(energy+p)*v])
        change=change-stop*boundary_flux
    drift=np.max(np.abs(change))
    assert drift<(2e-11 if axis<10 else 1e-8),(folder,drift)
    if axis==5 and levels>1:
        import importlib.util
        spec=importlib.util.spec_from_file_location('datutil',ROOT/'tools/python/amrvac_pytools/datfiles/reading/datfile_utilities.py')
        datutil=importlib.util.module_from_spec(spec);spec.loader.exec_module(datutil)
        snapshots=sorted(folder.glob('solution*.dat'))
        trees=[]
        for snapshot in (snapshots[0],snapshots[-1]):
            with snapshot.open('rb') as stream:
                level,index,_=datutil.get_tree_info(stream)
            trees.append(set(map(tuple,index[level==1])))
        assert trees[0]-trees[1] and trees[1]-trees[0], 'Expected both refinement and coarsening'
    return {**({'l1':stats[-1,1],'linf':stats[-1,2]} if axis<8 else {}),'conservation_drift':drift,
            'ghost':int(stats[-1,-1]),'fallbacks':fallbacks,
            'seconds':time.monotonic()-start,'cells':cells,'rows':len(cells)}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--binary',type=Path,default=ROOT/'tests/weno_js5/hd/agile')
    ap.add_argument('--output',type=Path,default=Path('/tmp/agile-weno-validation'))
    ap.add_argument('--convergence',action='store_true')
    ap.add_argument('--quick',action='store_true')
    ap.add_argument('--layout',choices=['per_cell','shared'],default='per_cell')
    args=ap.parse_args()
    results={}
    def run(name,**kw):
        data=run_case(args.binary.resolve(),args.output/name,layout=args.layout,**kw)
        results[name]={k:v for k,v in data.items() if k!='cells'}
        args.output.mkdir(parents=True,exist_ok=True)
        (args.output/'results.json').write_text(json.dumps(results,indent=2))
        print(name,results[name],flush=True)
        return data
    for lim in ('vanleer','weno5','mp5'):
        reference=run(f'{lim}-g3',limiter=lim,ghost=3)
        for name,kw in [('g4',{'ghost':4}),('mpi',{'ranks':2,'ghost':3}),('blocks',{'block':8,'ghost':3})]:
            other=run(f'{lim}-{name}',limiter=lim,**kw)
            np.testing.assert_allclose(other['cells'],reference['cells'],rtol=3e-13,atol=3e-13)
    for axis in (2,3):
        directional=run(f'axis-{axis}',axis=axis)
        np.testing.assert_allclose([directional['l1'],directional['linf']],
                                   [reference['l1'],reference['linf']],atol=1e-13,rtol=1e-10)
    run('negative-advection',axis=-1)
    run('auto-halo')
    run('reflect',axis=7,reflect=True)
    ref=run('reflect-wave-g3',axis=8,reflect=True,ghost=3)
    other=run('reflect-wave-g4',axis=8,reflect=True,ghost=4,ranks=2)
    np.testing.assert_allclose(ref['cells'],other['cells'],rtol=2e-12,atol=2e-12)
    for ghost in (3,4):
        ref=run(f'amr-g{ghost}',levels=2,ghost=ghost,axis=6,n=32,block=8,stop=.002)
        other=run(f'amr-mpi-g{ghost}',levels=2,ghost=ghost,axis=6,n=32,block=8,stop=.002,ranks=2)
        np.testing.assert_allclose(ref['cells'],other['cells'],rtol=2e-12,atol=2e-12)
    run('mixed',levels=2,limiter=['vanleer','mp5']+['weno5']*18,axis=6,block=8)
    run('bad-epsilon',epsilon=0.,expected_error='weno_js_epsilon must be finite and positive')
    run('no-energy',energy=False,expected_error='mp5 requires the energy equation')
    for label,eps in [('negative-epsilon',-1.),('nan-epsilon',float('nan')),('infinite-epsilon',float('inf'))]:
        run(label,epsilon=eps,expected_error='weno_js_epsilon must be finite and positive')
    for tol in (-1.,float('nan'),float('inf')):
        run('bad-tolerance-'+str(tol),tolerance=tol,expected_error='mp5_tolerance must be finite and nonnegative')
    run('zero-tolerance',tolerance=0.)
    run('bad-limiter',limiter='wenoz5',expected_error='Unsupported FV limiter')
    run('bad-block',levels=2,block=4,expected_error='AMR requires block_nx')
    if not args.quick:
        run('multidimensional',axis=4,n=16,block=8)
        run('dynamic-amr',levels=2,axis=5,n=32,block=8,stop=.08,dt=.001)
        for flux in ('tvdlf','hll','hllc'):
            run('sod-'+flux,axis=10,n=128,stop=.1,dt=.0002,flux=flux)
        run('shu-osher',axis=11,n=256,stop=.2,dt=.0005)
    if args.convergence:
        errors=[]
        for n in (32,64,128,256):
            case=run(f'convergence-{n}',n=n,stop=.005,dt=.1/n**2)
            errors.append([case['l1'],case['linf']])
        orders=np.log2(np.array(errors[:-1])/errors[1:])
        fine=run('convergence-half-dt',n=256,stop=.005,dt=.05/256**2)
        np.testing.assert_allclose([fine['l1'],fine['linf']],errors[-1],rtol=.1,atol=2e-14)
        assert np.all(orders[-1]>=4.5),orders
        results['orders']=orders.tolist()
        (args.output/'results.json').write_text(json.dumps(results,indent=2))
        print('orders',orders,flush=True)

if __name__=='__main__': main()

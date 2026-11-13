#!/usr/bin/env python3
"""Spatial RHS accuracy after integrating the single-face layouts."""
import argparse,json
import numpy as np
from pathlib import Path
from run_validation import run_case,ROOT
p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);a=p.parse_args();summary={}
for method in ('weno5','mp5'):
 errors=[]
 for n in (32,64,128,256):
  reference=None
  for layout in ('per_cell','shared'):
   r=run_case((ROOT/'tests/weno_js5/hd/agile').resolve(),a.output/f'{method}-{layout}-{n}',n=n,limiter=method,layout=layout,stepper='onestep',stop=1e-4,dt=1e-4)
   cells=r['cells'];x=cells[:,0];initial=1+.2*np.sinc(1/n)*np.sin(2*np.pi*x)
   rhs=(cells[:,4]-initial)/1e-4;exact=-.4*np.pi*np.sinc(1/n)*np.cos(2*np.pi*x)
   if reference is None:reference=cells
   else:np.testing.assert_array_equal(reference,cells)
  e=np.abs(rhs-exact);errors.append([float(e.mean()),float(e.max())])
 orders=np.log2(np.array(errors[:-1])/errors[1:]);assert np.all(orders[-1]>4.5),orders
 summary[method]={'errors':errors,'orders':orders.tolist(),'layouts_bitwise_identical':True};print(method,summary[method],flush=True)
(a.output/'summary.json').write_text(json.dumps(summary,indent=2))

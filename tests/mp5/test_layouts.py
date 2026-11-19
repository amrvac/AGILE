#!/usr/bin/env python3
"""Compare face layouts and optionally the pre-integration executable."""
import argparse,json
from pathlib import Path
import numpy as np
from run_validation import run_case,ROOT
p=argparse.ArgumentParser();p.add_argument('--previous',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
binary=(ROOT/'tests/weno_js5/hd/agile').resolve();rows=[]
for method in ('vanleer','koren','weno5','mp5'):
 for axis in (1,2,3):
  for flux in ('tvdlf','hll','hllc'):
   key=f'{method}-d{axis}-{flux}';kw=dict(limiter=method,axis=axis,flux=flux)
   per=run_case(binary,a.output/(key+'-per_cell'),layout='per_cell',**kw)
   shared=run_case(binary,a.output/(key+'-shared'),layout='shared',**kw)
   np.testing.assert_array_equal(per['cells'],shared['cells'])
   row=dict(case=key,layout_max_difference=0.,ghost=per['ghost'])
   if a.previous:
    old=run_case(a.previous.resolve(),a.output/(key+'-previous'),layout=None,**kw)
    np.testing.assert_allclose(per['cells'],old['cells'],atol=3e-13,rtol=3e-13)
    row['previous_max_difference']=float(np.max(abs(per['cells']-old['cells'])))
   rows.append(row);print(row,flush=True)
run_case(binary,a.output/'invalid-layout',layout='invalid',expected_error='unknown fv_face_flux')
(a.output/'summary.json').write_text(json.dumps(rows,indent=2))

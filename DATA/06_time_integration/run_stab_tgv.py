#!/usr/bin/env python3
# stability limit of the RK3 + sixth order KEP scheme: TGV 32^3 (acoustic modes present), fixed time
import os, sys, numpy as np
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import run_rk as R
U0 = 1.1832159566199232; c = np.sqrt(1.4 * 100.0); dx = 2 * np.pi / 32
S = (2 * (U0 + c) + c) / dx
jobs = []
for cfl in np.arange(0.6, 2.01, 0.1):
    dt = cfl / S; n = int(np.ceil(2.0 / dt)); dt = 2.0 / n
    jobs.append((f'stabtgv_cfl{cfl:.2f}', R.TG, dt, n, 'KEEP', 0, '', 4))
with ThreadPoolExecutor(max_workers=2) as ex:
    for name, rc in ex.map(R.run, jobs): print(name, rc, flush=True)

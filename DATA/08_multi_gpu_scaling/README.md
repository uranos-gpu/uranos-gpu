# Multi-GPU scaling (Sec. 4.2, Figs. 26-27, Table 9)

Triperiodic flow at rest, `(2 pi)^3`, 1024 x 1024 x 512 nodes, viscous terms on,
CFL = 0.5, no output in the time loop. Four GPUs per node.

- strong scaling `scal_S_<C>_n<nodes>.dat`: the grid stays at 1024 x 1024 x 512;
- weak scaling `scal_W_<C>_n<nodes>.dat`: 1024 x 1024 x (512 nodes);
- `*_w.dat`: 10-iteration warm-up of the same configuration;
- `nodes` = 1, 2, 4, ..., 256.

| config | convective | gas | grid | viscous |
|---|---|---|---|---|
| C1 | JP | CP | Cartesian | staggered |
| C2 | JP + WENO-Z 5 on every cell | CP | Cartesian | staggered |
| C3 | JP + WENO-Z 5 on every cell | TP | Cartesian | staggered |
| C4 | JP + WENO-Z 5 on every cell | CP | curved | staggered |
| C5 | AEC + TENO 7 on every cell | TP | curved | conservative |
| C6 | JP | CP | curved | staggered |

Curved grid: `grids/wavy_1024.txt` (`gen_curved_mesh.py`, amplitude 0.2).

```bash
mpirun -np 16 ./Uranos.exe DATA/08_multi_gpu_scaling/decks/scal_W_C5_n4_w.dat   # warm-up
mpirun -np 16 ./Uranos.exe DATA/08_multi_gpu_scaling/decks/scal_W_C5_n4.dat     # timed
```

The time per iteration is printed at the end of the run. The paper reports the
mean of three timed runs of 100 iterations after the warm-up.
`make_case.py` regenerates decks and mesh.

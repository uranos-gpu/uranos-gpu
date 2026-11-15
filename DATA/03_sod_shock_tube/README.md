# Sod shock tube (Sec. 3.1.3, Figs. 6-7)

`rho_R/rho_L = 0.125`, `p_R/p_L = 0.1`, gas at rest, `-5 <= x <= 5`, 200 x 5 x 5
nodes on the fake-Cartesian mesh `mesh_fake_cartesian.txt`, `Dt = 1e-3` up to
`t = 1.4`. Sixth-order JP flux blended with each shock-capturing reconstruction
through the density-jump sensor.

| decks | gas |
|---|---|
| `weno3 weno5 weno7 teno5 teno7 tenoa5 .dat` | CP, gamma = 1.4 (Fig. 6) |
| `*_tp2000.dat` | TP, `T_L = 2000 K` (Fig. 7) |

```bash
python3 DATA/03_sod_shock_tube/run_sod.py                 # CP, six runs
python3 DATA/03_sod_shock_tube/run_sod.py 1 2000 _tp2000  # TP at T_L = 2000 K
python3 DATA/03_sod_shock_tube/plot_sod.py                # sod_profiles.pdf
python3 DATA/03_sod_shock_tube/plot_sod_cp_tp.py          # sod_cp_tp.pdf
```

`run_sod.py` builds each deck from `shock_tube_base.dat` plus the overrides
listed at its end, and the mesh with `gen_fake_cartesian_mesh.py`.

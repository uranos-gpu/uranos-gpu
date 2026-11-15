# Time integration (Sec. 3.3, Fig. 11)

Temporal self-convergence of the three-stage Runge-Kutta scheme: grid and spatial
operators fixed, `Dt` refined, error against a reference run with `Dt0/128`.

| decks | case |
|---|---|
| `dwave_{cart,wavy}_{KEEP,PEP}_r<k>.dat` | density wave, `u = 0.5`, 32 x 32, KEEP-CP and PEP-TP, Cartesian and wavy (`mesh_wavy_32.txt`) |
| `tgv_cart_KEEP_r<k>.dat` | viscous Taylor-Green vortex 32^3, Ma = 0.2, KEEP, Laplacian |
| `stab_{KEEP,PEP}_cfl<c>.dat`, `stabtgv_cfl<c>.dat` | stability limit at increasing CFL |

`r<k>` is the refinement of the time step, `Dt = Dt0/k`.

```bash
python3 DATA/06_time_integration/run_rk.py        # every dwave, tgv and stab run
python3 DATA/06_time_integration/run_stab_tgv.py  # stability scan on the TGV
python3 DATA/06_time_integration/plot_rk2.py      # rk_order.pdf (uses plot_rk.py)
```

`rk_lib.py` and `kolmo_lib.py` are helper modules (deck template, wavy mesh,
binary reader); `setup.json` stores the reference step sizes.

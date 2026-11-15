# Isentropic vortex (Sec. 3.1.1, Figs. 2-3)

Two-dimensional isentropic vortex (strength 5) advected by a uniform stream in a
periodic box `[-5, 5]^2`, 100 x 100 x 4 nodes, CFL = 0.5, sixth-order KEP fluxes.
Wavy grid: `mesh_wavy.txt`, Eq. (60) of the paper with `a_w = 0.6`.

| decks | figure |
|---|---|
| `{jp,keep,aec,pep}_{cart,wavy}_x2.dat` | Fig. 3: two crossings along x; JP, KEEP, AEC with the CP gas, PEP with the TP gas |
| `jp_wavy.dat`, `jp_wavy_y.dat` | Fig. 2: one crossing along x and along y, JP on the wavy grid |

```bash
python3 DATA/01_isentropic_vortex/run_vortex.py 2     # the eight *_x2 runs (Fig. 3)
python3 DATA/01_isentropic_vortex/run_vortex.py       # the eight one-crossing runs
python3 DATA/01_isentropic_vortex/plot_vortex.py 2    # vortex_profiles.pdf
python3 DATA/01_isentropic_vortex/plot_storyboard.py  # vortex_storyboard.pdf (runs jp_wavy_y if missing)
```

The run scripts rewrite the decks and the mesh before running. A single run:
`mpirun -np 1 ./Uranos.exe DATA/01_isentropic_vortex/pep_wavy_x2.dat`.
`vortex_lib.py` and `make_gif.py` are helper modules (deck template, mesh builder,
binary reader).

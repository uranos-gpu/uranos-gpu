# Kolmogorov flow (Sec. 3.2.1, Fig. 8)

Steady Kolmogorov flow (`n = 4`, forcing `1e-3`, Ma = 0.05, Re = 10), convective
terms off, sixth-order viscous schemes at `N_g = 48`, on the Cartesian grid and on
the wavy grid `mesh_wavy.txt`.

Decks: `{laplacian,staggered,conservative}_{real_cart,wavy}.dat`.

```bash
python3 DATA/04_kolmogorov_flow/run_kolmo.py    # six runs
python3 DATA/04_kolmogorov_flow/plot_kolmo.py   # kolmogorov_profiles.pdf
```

The order-of-accuracy study of the same flow (Table 4, Fig. 9) is block C of
`../02_order_of_accuracy`. `kolmo_lib.py` is the helper module (deck template,
mesh builder, exact solution).

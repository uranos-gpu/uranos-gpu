# Order of accuracy (Secs. 3.1.2 and 3.2.1, Tables 3-4, Figs. 4, 5, 9)

Error-scaling campaign over every pure scheme: 324 convergence curves, 2349 runs.

| block | case | schemes |
|---|---|---|
| `A_ep` | harmonic density wave, `rho = 1 + 0.2 sin(pi x)`, `u = 0.05` | JP, KEEP, AEC (CP and TP), PEP (TP), fd order 2, 4, 6 |
| `B_weno` | same wave, shock-capturing flux on every cell | WENO-Z 3/5/7, TENO 5/7, TENO-A 5, CP and TP |
| `C_viscous` | steady Kolmogorov flow, convection off | laplacian, staggered, conservative, fd order 2, 4, 6 |

Each curve runs on three mesh families: `real_cart` (analytic grid), `fake_cart`
(the same nodes read from a mesh file) and `wavy` (Eq. (60), `a_w = 0.1`), along
x, y and z for blocks A and B.

- `decks/<block>/<curve>/N<nnnn>.dat`: one deck per run.
- `grids/<mesh>_<nx>x<ny>x<nz>.txt`: one mesh file per geometry.

```bash
python3 DATA/02_order_of_accuracy/campaign.py plan                    # counts
python3 DATA/02_order_of_accuracy/campaign.py list --block A          # curve names
python3 DATA/02_order_of_accuracy/campaign.py run --curve A_KEEP_CP_O6_x_wavy \
        --launcher "mpirun -np 1"                                     # one curve
python3 DATA/02_order_of_accuracy/campaign.py run --block B --index 7 # same, by index
python3 DATA/02_order_of_accuracy/aggregate.py                        # orders tables
python3 DATA/02_order_of_accuracy/campaign.py prepare                 # rewrite decks and meshes
python3 DATA/02_order_of_accuracy/campaign.py gate                    # fake_cart == real_cart
```

`run` executes the whole ladder of a curve, reduces each run to its L1, L2 and
Linf density (or velocity) error and writes `runs/<block>/<curve>/scaling.txt`
and `results/<curve>.json`. `aggregate.py` collects the fitted orders into
`results/{A,B,C}_orders.txt` and `results/orders_table.tex`. Every run uses
one rank. The paper campaign ran one curve per A100 GPU.

# Taylor-Green vortex (Sec. 3.2.2, Fig. 10, Table 5)

Re = 1600, Ma = 0.1, Pr = 0.71, `(2 pi)^3` periodic box, constant viscosity,
CFL = 0.5, up to `t u_inf / delta_0 = 20`. 168 runs:

- grid: `cart` (analytic) or `wavy` (`grids/wavy_tgv_<N>.txt`, amplitude 0.2);
- convective flux and gas: JP, KEEP, AEC with CP and TP, PEP with TP;
- viscous scheme: `lap`, `stag`, `cons`;
- resolution: 64^3, 128^3, 256^3, 512^3.

Deck names: `decks/<grid>_<gas>_<flux>_<viscous>_N<nnnn>.dat`; `campaign.json`
lists every case. Run from the repository root:

```bash
mpirun -np 4  ./Uranos.exe DATA/05_taylor_green_vortex/decks/wavy_TP_PEP_cons_N0128.dat
mpirun -np 16 ./Uranos.exe DATA/05_taylor_green_vortex/decks/cart_CP_KEEP_lap_N0512.dat
python3 DATA/05_taylor_green_vortex/plot_tgv.py     # tgv_validation.pdf
```

Each run writes `out/<tag>/TGV_MONITOR.txt` (enstrophy, dissipation, mean total
energy). A run stopped by the wall-clock limit restarts from its last field:
append `DATA/05_taylor_green_vortex/out/<tag>/BINARY/<last>.bin` to the command.
The paper used 4 A100 GPUs up to 256^3 and 16 at 512^3.
`reference/dns_vanrees2011_Re1600_512.txt` is the reference DNS.
`make_case.py` regenerates decks, meshes and `campaign.json`.

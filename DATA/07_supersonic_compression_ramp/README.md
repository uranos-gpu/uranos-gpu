# Supersonic compression ramp (Sec. 3.5, Table 7, DNS)

24 deg compression ramp at Ma = 2.91, `T_inf = 107.1 K`, isothermal wall
`T_w/T_inf = 2.87`. Single body-fitted block `-29 <= x/delta_0 <= 12`,
wall-normal extent `6 delta_0`, span `2.2 delta_0`, rounded corner;
2560 x 384 x 192 nodes (`mesh/dns/MeshXY_ramp.txt`). Digital-filter turbulent
inflow, NSCBC outflow, supersonic top boundary, Ducros sensor at 0.2, CFL = 0.5.

| case | numerics |
|---|---|
| `dns_keep_cp` | KEEP-CP, WENO-Z 5, Laplacian |
| `dns_aec_tp` | AEC-TP, WENO-Z 7, conservative |
| `dns_pep_tp` | PEP-TP, WENO-Z 7, conservative |

Each case has a transient deck (`_init`, statistics off) and a statistics deck
(`_stats`) that restarts from it. Both use 32 ranks (`cart_dims = 8, 1, 4`):

```bash
mpirun -np 32 ./Uranos.exe DATA/07_supersonic_compression_ramp/dns_aec_tp_init.dat
mpirun -np 32 ./Uranos.exe DATA/07_supersonic_compression_ramp/dns_aec_tp_stats.dat \
       DATA/07_supersonic_compression_ramp/dns_aec_tp/BINARY/<last>.bin
```

`make_decks.py` regenerates the six decks. The mesh comes from the Fortran
generator in `mesh/`:

```bash
cd DATA/07_supersonic_compression_ramp/mesh
python3 make_ramp_mesh.py --retau-dir dns --inflow-retau 380 --view
```

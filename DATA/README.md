# URANOS-3.0 paper cases

Every folder reproduces one section of

> F. De Vanna, S. Regazzo, A. Aiello, S. Baaziz, G. Baldan, G. Coppola, E. Benini,
> *URANOS-3.0: GPU-accelerated structure-preserving solver for hypersonic flows on
> curvilinear/body-fitted grids*.

Each folder is self-contained: decks (`*.dat`), meshes, the script that generated
them and the scripts that post-process the runs into the paper figures.

| folder | paper | content |
|---|---|---|
| `01_isentropic_vortex` | Sec. 3.1.1, Figs. 2-3 | isentropic vortex, four KEP fluxes, Cartesian and wavy grid |
| `02_order_of_accuracy` | Secs. 3.1.2, 3.2.1, Tables 3-4, Figs. 4, 5, 9 | density-wave and Kolmogorov-flow order of accuracy (2349 runs) |
| `03_sod_shock_tube` | Sec. 3.1.3, Figs. 6-7 | Sod tube, six shock-capturing reconstructions, CP and TP gas |
| `04_kolmogorov_flow` | Sec. 3.2.1, Fig. 8 | steady Kolmogorov flow, three viscous schemes |
| `05_taylor_green_vortex` | Sec. 3.2.2, Fig. 10, Table 5 | Taylor-Green vortex at Re = 1600 (168 runs) |
| `06_time_integration` | Sec. 3.3, Fig. 11 | temporal order and stability of the RK3 scheme |
| `07_supersonic_compression_ramp` | Sec. 3.5, Table 7 (DNS) | 24 deg ramp at Ma = 2.91, body-fitted DNS |
| `08_multi_gpu_scaling` | Sec. 4.2, Figs. 26-27, Table 9 | strong and weak scaling, triperiodic flow at rest |

## Running

Build the solver (see the top-level `Makefile`), then run **from the repository
root**: every path in a deck (`mesh_file`, `data_dir`) is relative to it.

```bash
make comp=nvhpc mode=gpu          # NVIDIA GPUs (OpenACC), as in the paper
make comp=gnu   mode=cpu          # CPU, enough for folders 01, 03, 04, 06

mpirun -np <N> ./Uranos.exe DATA/<case>/<deck>.dat [restart.bin]
```

The solver writes its output in `DATA/<data_dir>`, i.e. in a subfolder of the case
(`DATA/<case>/<run>/`). The Python scripts need `numpy`, `matplotlib` and, for
some figures, `scipy` and a LaTeX installation.

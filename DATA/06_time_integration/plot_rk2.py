#!/usr/bin/env python3
# Paper figure fig:rk_order. (a) temporal error scaling of all cases against CFL;
# (b) stability region of the three stage Runge--Kutta scheme with the spectrum of the sixth
#     order central operators at the default CFL and Fourier numbers.
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location('prk', os.path.join(HERE, 'plot_rk.py'))
plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'xtick.top': True, 'ytick.right': True,
    'lines.linewidth': 0.8, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})
# reuse the error sweeps of plot_rk.py without plotting
src = open(os.path.join(HERE, 'plot_rk.py')).read()
src = src[:src.index("COL = {'KEEP'")]
ns = {'__file__': os.path.join(HERE, 'plot_rk.py')}; exec(compile(src, 'plot_rk.py', 'exec'), ns)
sweep, spectral_radius = ns['sweep'], ns['spectral_radius']

COL = {'KEEP': '#374a95', 'PEP': '#ab3013'}; LS = {'KEEP': '-', 'PEP': '--'}
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw=dict(wspace=0.25, width_ratios=[1, 1]))

# ---------------- (a) error scaling
for grid, mk in (('cart', 'o'), ('wavy', 's')):
    Sr = spectral_radius(grid)
    for ep in ('KEEP', 'PEP'):
        dt, err = sweep(f'dwave_{grid}_{ep}', grid)
        ax.loglog(dt * Sr, err, color=COL[ep], ls=LS[ep], marker=mk, ms=3.0, mfc='white', mew=0.6, zorder=3 if ep == 'PEP' else 2)
U0 = 1.1832159566199232; c = np.sqrt(1.4 * 100.0); dx = 2 * np.pi / 32
dt, err = sweep('tgv_cart_KEEP', 'cart')
ax.loglog(dt * (2 * (U0 + c) + c) / dx, err, color=COL['KEEP'], ls='-', marker='D', ms=2.8, mfc='white', mew=0.6)
# slope triangle, order 3, right below the wavy grid curve
cw = lambda x: 4.76e-7 * (x / 1.783)**3
x0, x1 = 0.28, 0.56; ytop = 0.42 * cw(x1); y0 = ytop / (x1 / x0)**3
ax.fill([x0, x1, x1], [y0, y0, ytop], color='0.92', lw=0, zorder=1)
ax.plot([x0, x1, x1, x0], [y0, y0, ytop, y0], color='k', lw=0.6, zorder=2)
ax.text(x1 * 1.07, np.sqrt(y0 * ytop), r'$3$', fontsize=7, va='center')
ax.text(np.sqrt(x0 * x1), y0 * 0.62, r'$1$', fontsize=7, ha='center', va='top')
ax.set_xlim(0.05, 2.5); ax.set_ylim(4e-11, 1.5e-6)
ax.set_xlabel(r'$\mathrm{CFL}$'); ax.set_ylabel(r'$\mathcal{L}_2$')
ax.text(0.03, 0.95, '(a) temporal error', transform=ax.transAxes, va='top')

# ---------------- (b) stability region, upper half plane (symmetric), a few CFL numbers
R = lambda z: 1 + z + z**2 / 2 + z**3 / 6
X, Y = np.meshgrid(np.linspace(-2.8, 0.4, 900), np.linspace(-0.1, 3.0, 800))
A = np.abs(R(X + 1j * Y))
bx.contourf(X, Y, A, levels=[0, 1], colors=['#ececec'])
bx.contour(X, Y, A, levels=[1], colors='k', linewidths=1.0, zorder=8)
th = np.linspace(0, np.pi, 600)
kp = 2 * (3/4 * np.sin(th) - 3/20 * np.sin(2*th) + 1/60 * np.sin(3*th))
kpp = 49/18 - 2 * (3/2 * np.cos(th) - 3/20 * np.cos(2*th) + 1/90 * np.cos(3*th))
def reach(cfl, fo):
    # eigenvalues reached with the harmonic sum of the two limits: z = -b k''^2 + i a k', a/cfl + b/fo = 1
    xs = np.linspace(-fo * kpp.max(), 0, 500); ym = []
    for x in xs:
        ok = fo * kpp >= -x - 1e-12
        with np.errstate(divide='ignore', invalid='ignore'):
            v = cfl * kp[ok] * (1 + x / (fo * kpp[ok]))
        v = v[np.isfinite(v)]
        ym.append(v.max() if v.size else cfl * kp.max())
    return xs, np.maximum(np.array(ym), 0)
# CFL and Fo increased together from the default (0.5, 0.1); the limit is CFL = 1.09, Fo = 0.218
# CFL and Fo increased together from the default (0.5, 0.1); the limit is CFL = 1.09, Fo = 0.218
regs = [(1.3, 0.26, '#ab3013', '#f1d2cb', '-', 'unstable'), (1.09, 0.218, 'k', '#dcdcdc', '--', 'limit'),
        (0.8, 0.16, '#3db81d', '#d3efc8', '-', ''), (0.5, 0.1, '#374a95', '#c9d2ea', '-', 'default')]
bounds = []
for k, (cfl, fo, col, tint, ls, tag) in enumerate(regs):
    xs, ym = reach(cfl, fo); bounds.append((xs, ym))
    bx.fill_between(xs, 0, ym, color=tint, lw=0, zorder=2 + k)
    bx.plot(np.r_[xs, 0], np.r_[ym, 0], color=col, ls=ls, lw=1.0, zorder=7)
# part of the unstable area outside the stability region, in solid red
xs_u, ym_u = bounds[0]
inside_u = (X >= xs_u[0]) & (X <= 0) & (Y <= np.interp(X, xs_u, ym_u))
bx.contourf(X, Y, (inside_u & (A > 1)).astype(float), levels=[0.5, 1.5], colors=['#ab3013'], zorder=6)
bx.text(-1.25, 1.45, r'$|R(z)|\le 1$', fontsize=7.5, ha='center', zorder=10)
bx.xaxis.set_minor_locator(plt.MultipleLocator(0.1)); bx.yaxis.set_minor_locator(plt.MultipleLocator(0.1))
bx.tick_params(which='minor', length=1.5, width=0.5)
# label every area inside its own band, rotated along it, as (CFL, Fo)
bx.set_xlim(-1.8, 0.3); bx.set_ylim(0, 2.75)
fig.canvas.draw()
ylab = [0.62, 0.70, 0.62, 0.30]
for k, (cfl, fo, col, tint, ls, tag) in enumerate(regs):
    xs, ym = bounds[k]; y0 = ylab[k]
    xo = float(np.interp(y0, ym, xs))
    xi = float(np.interp(y0, bounds[k + 1][1], bounds[k + 1][0])) if k < 3 else 0.0
    x0 = 0.5 * (xo + xi)
    i0 = int(np.argmin(abs(ym - y0)))
    p1 = bx.transData.transform((xs[i0 - 20], ym[i0 - 20])); p2 = bx.transData.transform((xs[i0 + 20], ym[i0 + 20]))
    ang = np.degrees(np.arctan2(p2[1] - p1[1], p2[0] - p1[0]))
    txt = rf'$({cfl:g},\,{fo:g})$' + (f' {tag}' if tag and tag != 'default' else '')
    if tag == 'default': txt = r'\textbf{default}' + '\n' + txt
    bx.text(x0, y0, txt, rotation=ang, rotation_mode='anchor', ha='center', va='center', linespacing=1.1, multialignment='center',
            fontsize=6.6, color=col if col != 'k' else '0.1', zorder=10)
bx.text(0.97, 0.97, r'labels: $(\mathrm{CFL},\,\mathrm{Fo})$', transform=bx.transAxes, va='top', ha='right', fontsize=7)
xr = -2.5127453266
bx.plot([0], [np.sqrt(3)], 'ko', ms=2.4, zorder=4)
bx.text(0.04, np.sqrt(3), r'$\mathrm{i}\sqrt{3}$', fontsize=7, va='center')
bx.annotate(r'$-2.51$', xy=(-1.8, 0.16), xytext=(-1.6, 0.16), fontsize=7, va='center',
            arrowprops=dict(arrowstyle='->', lw=0.5, color='k'))


bx.xaxis.set_major_locator(plt.MultipleLocator(0.5))
bx.set_xlabel(r'$\mathrm{Re}(\lambda\,\Delta t)$'); bx.set_ylabel(r'$\mathrm{Im}(\lambda\,\Delta t)$')
bx.text(0.03, 0.97, '(b) stability region', transform=bx.transAxes, va='top', ha='left')

h = [plt.Line2D([], [], color=COL['KEEP'], ls='-', label=r'\texttt{KEEP}, calorically perfect'),
     plt.Line2D([], [], color=COL['PEP'], ls='--', label=r'\texttt{PEP}, thermally perfect'),
     plt.Line2D([], [], color='0.3', ls='none', marker='o', ms=3.0, mfc='white', mew=0.6, label='wave, Cartesian'),
     plt.Line2D([], [], color='0.3', ls='none', marker='s', ms=3.0, mfc='white', mew=0.6, label='wave, wavy'),
     plt.Line2D([], [], color='0.3', ls='none', marker='D', ms=2.8, mfc='white', mew=0.6, label=r'Taylor--Green, $32^3$'),
]
leg = fig.legend(handles=h, loc='upper center', ncol=3, frameon=True, fancybox=True, edgecolor='k',
                 bbox_to_anchor=(0.5, -0.01), handlelength=2.8, columnspacing=1.6)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'rk_order.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
#print('max k\'', kp.max(), 'max k\'\'^2', kpp.max(), 'CFL_max', yi / kp.max(), 'Fo_max', -xr / kpp.max())

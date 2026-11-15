#!/usr/bin/env python3
# Paper figure: isentropic vortex, central (KEP) schemes on the Cartesian and wavy grids.
import os, sys, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.ticker
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, HERE)
import vortex_lib as rv
import make_gif as mg

plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in', 'lines.linewidth': 0.8,
    'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})

G = 1.4; BETA = 5.0; ALPHA = BETA**2*(G-1)/(8*G*np.pi**2)
NCROSS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
SUF = '' if NCROSS == 1 else f'_x{NCROSS}'
N = rv.N_INPLANE; T1 = rv.TMAX; T = NCROSS*T1; L = rv.LX
TSTOP = float(sys.argv[2]) if len(sys.argv) > 2 else NCROSS   # plot the runs only up to this crossing
TLAB = f'{TSTOP:g}'
# colors of the sixth order central schemes in the order of accuracy figure; dashed = thermally perfect
CASES = [('jp', r'\texttt{JP}', '#4b0082', '-'), ('keep', r'\texttt{KEEP}', '#374a95', '-'),
         ('aec', r'\texttt{AEC}', '#3db81d', '-'), ('pep', r'\texttt{PEP}', '#ab3013', '--')]

def grid(kind):
    XX, YY = rv.build_wavy(N, N) if kind == 'wavy' else (None, None)
    if kind == 'cart':
        amp = rv.AMP; rv.AMP = 0.0; XX, YY = rv.build_wavy(N, N); rv.AMP = amp
    ng = rv.NG
    X, Y = XX[ng:-ng, ng:-ng], YY[ng:-ng, ng:-ng]
    xe, xn = np.gradient(XX, axis=0), np.gradient(XX, axis=1)
    ye, yn = np.gradient(YY, axis=0), np.gradient(YY, axis=1)
    dV = (xe*yn - xn*ye)[ng:-ng, ng:-ng]
    return X, Y, dV

tab = np.loadtxt(os.path.join(HERE, 'pep_cart', 'thermo_profiles.txt'))
cpr = np.interp(1.0, tab[:, 1], tab[:, 2]); G_TP = cpr/(cpr-1.0)   # gamma_ref of the thermally perfect runs

def exact_rho(X, Y, t, g=G):
    alpha = BETA**2*(g-1)/(8*g*np.pi**2); uinf = np.sqrt(g)*rv.MACH
    xc = uinf*t
    dx = X - xc; dx -= L*np.round(dx/L)
    return (1 - alpha*np.exp(1 - dx**2 - Y**2))**(1/(g-1))

def series(name, X, Y, dV, g=G):
    files = sorted(glob.glob(os.path.join(HERE, name, 'BINARY', '*.bin')))
    tt, e2, snaps = [], [], []
    for f in files:
        t, nx, ny, nz, phi = mg.read_binary(f)
        if t/T1 > TSTOP + 0.05: break
        r = phi[0].mean(axis=2)
        err = r - exact_rho(X, Y, t, g)
        tt.append(t/T1); e2.append(np.sqrt(np.sum(err**2*dV)/np.sum(dV))); snaps.append((abs(t/T1 - TSTOP), t, r))
    _, t, r = min(snaps, key=lambda q: q[0])                   # snapshot closest to TSTOP
    return np.array(tt), np.array(e2), r, t

data = {}
for kind in ('cart', 'wavy'):
    X, Y, dV = grid(kind)
    j0 = np.argmin(np.abs(Y[N//2, :]))                    # grid line through the center
    for key, *_ in CASES:
        tt, e2, rho, tf = series(f'{key}_{kind}{SUF}', X, Y, dV, G_TP if key == 'pep' else G)
        data[(key, kind)] = (tt, e2, rho[:, j0], X[:, j0], Y[:, j0], tf)

def style_inset(a):
    a.tick_params(labelsize=6, length=2, direction='in', top=True, right=True)
    for s_ in a.spines.values(): s_.set_linewidth(0.5)

LW = {'jp': 0.8, 'keep': 0.8, 'aec': 0.8, 'pep': 0.8}
Z = {'jp': 3, 'keep': 2, 'aec': 1, 'pep': 4}
OFF = {'jp': 0, 'keep': 1, 'aec': 2, 'pep': 3}
MRK = {'cart': 'o', 'wavy': 's'}
def mk(key, kind, start, step):
    # open markers as in the error scaling figure: circle Cartesian, square wavy; staggered per scheme
    return dict(marker=MRK[kind], ms=2.8, mfc='white', mew=0.6, markevery=(start + OFF[key]*step/4, step))
fig = plt.figure(figsize=(7.0, 6.6))
gs = fig.add_gridspec(3, 2, hspace=0.36, wspace=0.24)
labels = iter('abcdef')
for col, kind, title in ((0, 'cart', 'Cartesian grid'), (1, 'wavy', 'wavy grid')):
    ax = fig.add_subplot(gs[0, col])
    xs, ys = data[('jp', kind)][3], data[('jp', kind)][4]
    xf = np.linspace(xs.min(), xs.max(), 2001); yf = np.interp(xf, xs, ys)
    ax.plot(xf, exact_rho(xf, yf, data[('jp', kind)][5]), color='k', lw=0.8)
    for key, lab, c, ls in CASES:
        ax.plot(data[(key, kind)][3], data[(key, kind)][2], color=c, ls=ls, lw=LW[key], zorder=Z[key], **mk(key, kind, 0.03, 0.26))
    ax.set_xlim(-5, 5); ax.set_ylim(0.45, 1.07)
    ax.set_xlabel(r'$x/\delta_0$'); ax.set_ylabel(r'$\rho/\rho_\infty$')
    a = ax.inset_axes([0.66, 0.10, 0.31, 0.46])
    a.plot(xf, exact_rho(xf, yf, data[('jp', kind)][5]), color='k', lw=0.8)
    for key, lab, c, ls in CASES:
        a.plot(data[(key, kind)][3], data[(key, kind)][2], color=c, ls=ls, lw=0.8, zorder=Z[key], **mk(key, kind, 0.05, 0.3))
    tj = data[('jp', kind)][5]; xc = (np.sqrt(G)*rv.MACH*tj - rv.XMIN) % L + rv.XMIN   # exact core position
    a.set_xlim(xc - 0.35, xc + 0.35); a.set_ylim(0.492, 0.512)
    style_inset(a); mark_inset(ax, a, loc1=2, loc2=3, fc='none', ec='0.45', lw=0.5)
    ax.text(0.03, 0.10, f'({next(labels)}) {title}' + '\n' + rf'$t\,u_\infty/L_x={TLAB}$', transform=ax.transAxes, linespacing=1.4)
for col, kind, title in ((0, 'cart', 'Cartesian grid'), (1, 'wavy', 'wavy grid')):
    ax = fig.add_subplot(gs[1, col])
    for key, lab, c, ls in CASES:
        xs, ys = data[(key, kind)][3], data[(key, kind)][4]
        g = G_TP if key == 'pep' else G
        ax.semilogy(xs, np.abs(data[(key, kind)][2] - exact_rho(xs, ys, data[(key, kind)][5], g)), color=c, ls=ls, lw=LW[key], zorder=Z[key], **mk(key, kind, 0.03, 0.26))
    ax.set_xlim(-5, 5); ax.set_ylim(1e-9, 10 if TSTOP > 2 else 1e-3)
    ax.set_xlabel(r'$x/\delta_0$'); ax.set_ylabel(r'$|\rho-\rho_{\mathrm{exact}}|/\rho_\infty$')
    ax.text(0.03, 0.88, f'({next(labels)}) {title}, ' + rf'$t\,u_\infty/L_x={TLAB}$', transform=ax.transAxes)
for col, kind, title in ((0, 'cart', 'Cartesian grid'), (1, 'wavy', 'wavy grid')):
    ax = fig.add_subplot(gs[2, col])
    for key, lab, c, ls in CASES:
        tt, e2 = data[(key, kind)][:2]
        ax.semilogy(tt, e2, color=c, ls=ls, lw=LW[key], zorder=Z[key], **mk(key, kind, 0.04, 0.26))
    ax.set_xlim(0, TSTOP); ax.set_ylim(*((1e-9, 1e-1) if TSTOP > 2 else (1e-10, 1e-4)))
    if TSTOP <= 2: ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
    ax.set_xlabel(r'$t\,u_\infty/L_x$'); ax.set_ylabel(r'$\mathcal{L}_2$')
    ax.text(0.03, 0.06 if TSTOP > 2 else 0.88, f'({next(labels)}) {title}', transform=ax.transAxes)
    # inset on the end of the run, linear axis with the order of magnitude on top
    t0 = 0.8*TSTOP
    a = ax.inset_axes([0.50, 0.09, 0.46, 0.36])
    lo, hi = np.inf, 0.0
    for key, lab, c, ls in CASES:
        tt, e2 = data[(key, kind)][:2]; m = tt >= t0
        a.plot(tt[m], e2[m], color=c, ls=ls, lw=0.8, zorder=Z[key], **mk(key, kind, 0.05, 0.3))
        lo, hi = min(lo, e2[m].min()), max(hi, e2[m].max())
    pad = 0.08*(hi - lo); a.set_xlim(t0, TSTOP); a.set_ylim(lo - pad, hi + pad)
    fm = mticker.ScalarFormatter(useMathText=True); fm.set_powerlimits((0, 0)); a.yaxis.set_major_formatter(fm)
    a.yaxis.get_offset_text().set_fontsize(6)
    a.set_xticks(np.linspace(t0, TSTOP, 3)); a.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.1f'))
    a.yaxis.set_major_locator(plt.MaxNLocator(3))
    style_inset(a); mark_inset(ax, a, loc1=1, loc2=2, fc='none', ec='0.45', lw=0.5)
handles = [plt.Line2D([], [], color='k', lw=0.8, label='Exact')] + \
          [plt.Line2D([], [], color=c, ls=ls, lw=0.8, label=l) for key, l, c, ls in CASES] + \
          [plt.Line2D([], [], color='0.3', ls='none', marker=m, ms=2.8, mfc='white', mew=0.6, label=l)
           for m, l in (('o', 'Cartesian'), ('s', 'wavy'))]
leg = fig.legend(handles=handles, loc='upper center', ncol=7, frameon=True, fancybox=True, edgecolor='k',
                 bbox_to_anchor=(0.5, 0.035), handlelength=2.8, columnspacing=1.6)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'vortex_profiles.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
for (key, kind), v in data.items():
    print(key, kind, 'final t/T', round(v[0][-1], 4), 'L2 final', f'{v[1][-1]:.3e}', 'rho min', round(v[2].min(), 4))

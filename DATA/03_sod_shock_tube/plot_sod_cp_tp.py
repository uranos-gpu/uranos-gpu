#!/usr/bin/env python3
# Paper figure: Sod shock tube, calorically vs thermally perfect gas (T_L = 2000 K).
import os, glob, runpy
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
base = runpy.run_path(os.path.join(HERE, 'plot_sod.py'), run_name='lib')   # reuse style, reader, CP exact
read_bin, x, CASES, exact_cp, G = base['read_bin'], base['x'], base['CASES'], base['exact_sod'], base['G']
TAG = '_tp2000'

# ---- thermally perfect thermodynamics from the solver's own table (dimensionless, R = 1) ----
tab = np.loadtxt(os.path.join(HERE, 'weno5' + TAG, 'thermo_profiles.txt'))
TH, CP, H, CHI = tab[:, 1], tab[:, 2], tab[:, 3], tab[:, 4]
h   = lambda th: np.interp(th, TH, H)
cp  = lambda th: np.interp(th, TH, CP)
chi = lambda th: np.interp(th, TH, CHI)
e_of = lambda th: h(th) - th
th_of_e   = lambda e: np.interp(e, H - TH, TH)
th_of_chi = lambda c: np.interp(c, CHI, TH)
gam = lambda th: cp(th)/(cp(th) - 1.0)

def exact_tp(t):
    rL, pL, rR, pR = 1.0, 1.0, 0.125, 0.1
    thL, thR = pL/rL, pR/rR
    cL = np.sqrt(gam(thL)*thL)
    sL = chi(thL) - np.log(pL)
    # left rarefaction: isentrope, u(p) = int_p^pL dp'/(rho c)
    pp = np.linspace(pL, 0.05, 20001)
    thp = th_of_chi(sL + np.log(pp)); rp = pp/thp; cpp = np.sqrt(gam(thp)*thp)
    up = np.concatenate(([0.0], np.cumsum(0.5*(1/(rp[1:]*cpp[1:])+1/(rp[:-1]*cpp[:-1]))*(pp[:-1]-pp[1:]))))
    uL_of_p = lambda p: np.interp(-p, -pp, up)
    # right shock: Rankine-Hugoniot with h(theta)
    def shock(p2):
        f = lambda th2: h(th2)-h(thR) - 0.5*(p2-pR)*(1/rR + th2/p2)
        th2 = brentq(f, thR, 5*thR); r2 = p2/th2
        m = np.sqrt((p2-pR)/(1/rR-1/r2))
        return r2, (p2-pR)/m, m/rR           # density, velocity jump, shock speed
    pm = brentq(lambda p: uL_of_p(p) - shock(p)[1], 0.11, 0.99)
    um = uL_of_p(pm); rmR, _, SS = shock(pm)
    sel = pp >= pm
    xf = (up[sel] - cpp[sel])*t
    xs = np.concatenate(([-5.0], xf, [um*t, um*t, SS*t, SS*t, 5.0]))
    rho = np.concatenate(([rL], rp[sel], [rp[sel][-1], rmR, rmR, rR, rR]))
    vel = np.concatenate(([0.0], up[sel], [um, um, um, 0.0, 0.0]))
    tem = np.concatenate(([thL], thp[sel], [thp[sel][-1], pm/rmR, pm/rmR, thR, thR]))
    return xs, np.column_stack((rho, vel, tem)), dict(pm=pm, um=um, xc=um*t, xs=SS*t, xt=(um-np.sqrt(gam(thp[sel][-1])*thp[sel][-1]))*t)

def profiles(key, tp):
    t, phi = read_bin(sorted(glob.glob(os.path.join(HERE, key, 'BINARY', '*.bin')))[-1])
    rho = phi[0].mean(axis=(0, 1)); u = (phi[1]/phi[0]).mean(axis=(0, 1))
    e = (phi[4].mean(axis=(0, 1)) - 0.5*rho*u**2)/rho
    th = th_of_e(e) if tp else (G-1)*e
    return t, rho, u, th

t = profiles('weno5', False)[0]
xcp, excp, Wcp = exact_cp(t)
excp = excp.copy(); excp[:, 2] = excp[:, 2]/excp[:, 0]           # p -> theta = p/rho
xtp, extp, Wtp = exact_tp(t)
CPd = {k: profiles(k, False)[1:] for k, *_ in CASES}
TPd = {k: profiles(k + TAG, True)[1:] for k, *_ in CASES}

def draw(ax, j, lw=0.8):
    ax.plot(xcp, excp[:, j], color='k', lw=lw, ls='-', zorder=1)
    ax.plot(xtp, extp[:, j], color='k', lw=lw, ls='--', zorder=1)
    for key, lab, col, _ in CASES:
        ax.plot(x, TPd[key][j], color=col, lw=lw, ls='--', zorder=2)

def inset(ax, j, box, xlim, ylim, loc1, loc2):
    a = ax.inset_axes(box); draw(a, j, lw=0.8)
    a.set_xlim(*xlim); a.set_ylim(*ylim); a.tick_params(labelsize=6, length=2)
    for s_ in a.spines.values(): s_.set_linewidth(0.5)
    mark_inset(ax, a, loc1=loc1, loc2=loc2, fc='none', ec='0.45', lw=0.5)

fig = plt.figure(figsize=(7.0, 5.8))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.30, wspace=0.24)
ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, 0]); ax3 = fig.add_subplot(gs[1, 1])
xcC, xcT, xsC, xsT = Wcp['xc'], Wtp['xc'], Wcp['xs'], Wtp['xs']
tailC, tailT = Wcp['xt'], Wtp['xt']

# (a) density: rarefaction tail, contact, shock
draw(ax1, 0); ax1.set_xlim(-5, 5); ax1.set_ylim(0.05, 1.17)
ax1.set_xlabel(r'$x/\delta_0$'); ax1.set_ylabel(r'$\rho/\rho_L$')
inset(ax1, 0, [0.09, 0.08, 0.24, 0.42], (min(tailC, tailT)-0.25, max(tailC, tailT)+0.25), (0.395, 0.445), 2, 1)
inset(ax1, 0, [0.42, 0.56, 0.24, 0.40], (min(xcC, xcT)-0.35, max(xcC, xcT)+0.35), (0.25, 0.45), 3, 4)
inset(ax1, 0, [0.75, 0.56, 0.23, 0.40], (min(xsC, xsT)-0.25, max(xsC, xsT)+0.25), (0.10, 0.30), 3, 4)
ax1.text(0.012, 0.965, '(a)', transform=ax1.transAxes, va='top', ha='left')

# (b) velocity: post shock plateau and shock
draw(ax2, 1); ax2.set_xlim(-5, 5); ax2.set_ylim(-0.05, 1.12)
ax2.set_xlabel(r'$x/\delta_0$'); ax2.set_ylabel(r'$u/u_0$')
inset(ax2, 1, [0.15, 0.28, 0.20, 0.47], (min(xsC, xsT)-0.6, max(xsC, xsT)+0.2), (0.88, 1.00), 2, 1)
ax2.text(0.012, 0.965, '(b)', transform=ax2.transAxes, va='top', ha='left')

# (c) temperature: end of rarefaction, contact and shock
draw(ax3, 2); ax3.set_xlim(-5, 5); ax3.set_ylim(0.62, 1.55)
ax3.set_xlabel(r'$x/\delta_0$'); ax3.set_ylabel(r'$T/T_L$')
inset(ax3, 2, [0.14, 0.06, 0.20, 0.28], (-0.3, max(xcC, xcT)+0.1), (0.69, 0.79), 1, 4)
inset(ax3, 2, [0.50, 0.62, 0.47, 0.34], (min(xcC, xcT)-0.2, max(xsC, xsT)+0.2), (1.02, 1.20), 3, 4)
ax3.text(0.012, 0.965, '(c)', transform=ax3.transAxes, va='top', ha='left')

blank = plt.Line2D([], [], color='none', label=' ')
sch = [plt.Line2D([], [], color=c, lw=0.8, ls='--', label=l) for _, l, c, _ in CASES]
handles = [plt.Line2D([], [], color='k', lw=0.8, ls='-', label='Exact, calorically perfect'),
           plt.Line2D([], [], color='k', lw=0.8, ls='--', label='Exact, thermally perfect'), blank] + sch[:3] + sch[3:]
leg = fig.legend(handles=handles, loc='upper center', ncol=3, frameon=True, fancybox=True,
                 edgecolor='k', bbox_to_anchor=(0.5, 0.035), handlelength=2.8, columnspacing=1.6)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'sod_cp_tp.pdf'); fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
print('CP star', Wcp['pm'], Wcp['um'], Wcp['xc'], Wcp['xs'])
print('TP star', Wtp['pm'], Wtp['um'], Wtp['xc'], Wtp['xs'])
for k, lab, *_ in CASES:
    print(lab, 'TP rho range', TPd[k][0].min(), TPd[k][0].max(), 'L1 vs TP exact',
          np.mean(np.abs(TPd[k][0]-np.interp(x, xtp, extp[:, 0])))*10)

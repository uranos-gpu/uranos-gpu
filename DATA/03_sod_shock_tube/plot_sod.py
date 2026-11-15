#!/usr/bin/env python3
# Paper figure: Sod shock tube, WENO/TENO reconstructions vs exact solution.
import os, struct, glob
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

HERE = os.path.dirname(os.path.abspath(__file__))
NG, G = 7, 1.4
NX, NY = 200, 5

plt.rcParams.update({
    'text.usetex': True, 'font.family': 'serif', 'font.serif': ['Computer Modern Roman'],
    'font.size': 8, 'axes.labelsize': 8, 'legend.fontsize': 7,
    'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'lines.linewidth': 0.8, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02})

CASES = [('weno3', r'WENO-Z\,3',  '#73ade0', '-'),     # colors as in the order of accuracy figure
         ('weno5', r'WENO-Z\,5',  '#557bbb', '-'),
         ('weno7', r'WENO-Z\,7',  '#374a95', '-'),
         ('teno5', r'TENO\,5',    '#ca6026', '-'),
         ('teno7', r'TENO\,7',    '#ab3013', '-'),
         ('tenoa5', r'TENO-A\,5', '#4dca23', '-')]

def read_bin(path):
    with open(path, 'rb') as f:
        nx, ny, nz, it = struct.unpack('<4i', f.read(16)); t, dt = struct.unpack('<2d', f.read(16))
        phi = np.frombuffer(f.read(5*nx*ny*nz*8), dtype='<f8').reshape((5, nz, ny, nx))
    return t, phi

def read_x(path):
    v = [float(x) for x in open(path).read().split()[3:]]
    m = (NX+2*NG)*(NY+2*NG)
    xx = np.array(v[:m]).reshape((NY+2*NG, NX+2*NG))
    return xx[NG:-NG, NG:-NG].mean(axis=0)

def exact_sod(t):
    """Exact solution of the Sod Riemann problem (left rarefaction, contact, right shock),
    returned as a piecewise curve with the discontinuities as exact vertical jumps."""
    rL, pL, rR, pR = 1.0, 1.0, 0.125, 0.1
    cL, cR = np.sqrt(G*pL/rL), np.sqrt(G*pR/rR); Gm = (G-1)/(2*G)
    fL = lambda p: (2*cL/(G-1))*((p/pL)**Gm-1)                          # rarefaction
    fR = lambda p: (p-pR)*np.sqrt((2/((G+1)*rR))/(p+pR*(G-1)/(G+1)))    # shock
    dfL = lambda p: (1/(rL*cL))*(p/pL)**(-(G+1)/(2*G))
    dfR = lambda p: np.sqrt((2/((G+1)*rR))/(p+pR*(G-1)/(G+1)))*(1-(p-pR)/(2*(p+pR*(G-1)/(G+1))))
    pm = 0.3
    for _ in range(60):
        pm -= (fL(pm)+fR(pm))/(dfL(pm)+dfR(pm))
    um = 0.5*(fR(pm)-fL(pm))
    rmL = rL*(pm/pL)**(1/G)
    rmR = rR*((pm/pR)+(G-1)/(G+1))/(((G-1)/(G+1))*(pm/pR)+1)
    cmL = cL*(pm/pL)**Gm
    SH, ST = -cL, um-cmL                                   # rarefaction head and tail speeds
    SS = cR*np.sqrt((G+1)/(2*G)*pm/pR+Gm)                   # shock speed
    xf = np.linspace(SH*t, ST*t, 4001)                     # exact fan
    uf = (2/(G+1))*(cL+xf/t); cf = (2/(G+1))*(cL-(G-1)/2*xf/t)
    rf = rL*(cf/cL)**(2/(G-1)); pf = pL*(cf/cL)**(2*G/(G-1))
    xs = np.concatenate(([-5.0], xf, [um*t, um*t, SS*t, SS*t, 5.0]))
    rho = np.concatenate(([rL], rf, [rmL, rmR, rmR, rR, rR]))
    vel = np.concatenate(([0.0], uf, [um, um, um, 0.0, 0.0]))
    prs = np.concatenate(([pL], pf, [pm, pm, pm, pR, pR]))
    return xs, np.column_stack((rho, vel, prs)), dict(pm=pm, um=um, xh=SH*t, xt=ST*t, xc=um*t, xs=SS*t)

x = read_x(os.path.join(HERE, 'mesh_fake_cartesian.txt'))
data = {}
for key, *_ in CASES:
    t, phi = read_bin(sorted(glob.glob(os.path.join(HERE, key, 'BINARY', '*.bin')))[-1])
    rho = phi[0].mean(axis=(0, 1)); u = (phi[1]/phi[0]).mean(axis=(0, 1))
    p = ((G-1)*(phi[4]-0.5*(phi[1]**2+phi[2]**2+phi[3]**2)/phi[0])).mean(axis=(0, 1))
    data[key] = (rho, u, p)
xe, ex, W = exact_sod(t)
um, pm = W['um'], W['pm']
xc, xs, xt, xh = W['xc'], W['xs'], W['xt'], W['xh']

# ---- density jump sensor (as coded), evaluated on the WENO-Z 5 solution ----
MACH, EPS1, EPS2, WNUM = 1.0, 0.01, 0.05, 3
r5, u5, p5 = data['weno5']
Mloc = np.abs(u5)/np.sqrt(G*p5/r5)
toll = MACH*EPS1 + MACH*(EPS2-EPS1)/(1.0+(Mloc+0.2)**20)
jump = np.abs(np.append(np.diff(r5), 0.0))
shock_node = jump > toll
iface = np.array([shock_node[max(i-WNUM+1,0):i+WNUM+1].any() for i in range(len(x)-1)])
spans = []
i = 0
while i < len(iface):
    if iface[i]:
        j = i
        while j+1 < len(iface) and iface[j+1]: j += 1
        spans.append((x[i], x[j+1])); i = j+1
    else: i += 1

def shade(ax):
    for a_, b_ in spans:
        ax.axvspan(a_, b_, color='0.88', lw=0, zorder=0)

def draw(ax, j, lw=0.8):
    shade(ax)
    ax.plot(xe, ex[:, j], color='k', lw=lw, zorder=1)
    for key, lab, col, ls in CASES:
        ax.plot(x, data[key][j], color=col, ls=ls, lw=lw, zorder=2)

def inset(ax, j, box, xlim, ylim, loc1=2, loc2=4):
    a = ax.inset_axes(box)
    draw(a, j, lw=0.8)
    a.set_xlim(*xlim); a.set_ylim(*ylim)
    a.tick_params(labelsize=6, length=2)
    for s in a.spines.values(): s.set_linewidth(0.5)
    mark_inset(ax, a, loc1=loc1, loc2=loc2, fc='none', ec='0.45', lw=0.5)
    return a

fig = plt.figure(figsize=(7.0, 5.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.30, wspace=0.24)
ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, 0]); ax3 = fig.add_subplot(gs[1, 1])

# (a) density
draw(ax1, 0)
ax1.set_xlim(-5, 5); ax1.set_ylim(0.05, 1.17)
ax1.set_xlabel(r'$x/\delta_0$'); ax1.set_ylabel(r'$\rho/\rho_L$')
inset(ax1, 0, [0.07, 0.10, 0.22, 0.38], (xt-0.30, xt+0.30), (0.418, 0.442), 2, 1)
inset(ax1, 0, [0.42, 0.56, 0.22, 0.40], (xc-0.45, xc+0.45), (0.24, 0.45), 3, 4)
inset(ax1, 0, [0.76, 0.56, 0.22, 0.40], (xs-0.30, xs+0.30), (0.11, 0.29), 3, 4)
ax1.text(0.012, 0.965, '(a)', transform=ax1.transAxes, va='top', ha='left')

# (b) velocity
draw(ax2, 1)
ax2.set_xlim(-5, 5); ax2.set_ylim(-0.05, 1.12)
ax2.set_xlabel(r'$x/\delta_0$'); ax2.set_ylabel(r'$u/u_0$')
inset(ax2, 1, [0.15, 0.36, 0.23, 0.40], (xh-0.35, xh+0.35), (-0.02, 0.10), 3, 4)
ax2.text(0.012, 0.965, '(b)', transform=ax2.transAxes, va='top', ha='left')

# (c) pressure
draw(ax3, 2)
ax3.set_xlim(-5, 5); ax3.set_ylim(0.05, 1.17)
ax3.set_xlabel(r'$x/\delta_0$'); ax3.set_ylabel(r'$p/p_L$')
inset(ax3, 2, [0.58, 0.45, 0.38, 0.40], (xs-0.35, xs+0.25), (0.27, 0.32), 3, 4)
ax3.text(0.012, 0.965, '(c)', transform=ax3.transAxes, va='top', ha='left')

blank = plt.Line2D([], [], color='none', label=' ')
exact = plt.Line2D([], [], color='k', lw=0.8, label='Exact')
band = matplotlib.patches.Patch(color='0.88', label='Shock capturing active')
sch = [plt.Line2D([], [], color=c, lw=0.8, ls=ls, label=l) for _, l, c, ls in CASES]
handles = [exact, band, blank] + sch[:3] + sch[3:]
leg = fig.legend(handles=handles, loc='upper center', ncol=3, frameon=True, fancybox=True,
                 edgecolor='k', bbox_to_anchor=(0.5, 0.035), handlelength=2.8, columnspacing=1.6)
leg.get_frame().set_linewidth(0.6)
out = os.path.join(HERE, 'sod_profiles.pdf')
fig.savefig(out); fig.savefig(out.replace('.pdf', '.png'), dpi=200)
print('spans', [(round(a,3),round(b,3)) for a,b in spans]); print('t =', t, ' x_c =', xc, ' x_s =', xs, ' x_tail =', xt, ' x_head =', xh)
for key, lab, *_ in CASES:
    r = data[key][0]; print(f'{lab:10s} rho in [{r.min():.4f},{r.max():.4f}]')

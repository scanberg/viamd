import sys, math, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mdtraj as md, pickle
import story
from movie_lib import *

OUT = os.path.join(ROOT, 'docs/images/aspirin_movie')
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 110})

A = pickle.load(open(os.path.join(WORK, 'asp_analysis.pkl'), 'rb')); E = A['edr']
traj = md.load(os.path.join(DATA, 'aspirin-phospholipase.xtc'), top=os.path.join(DATA, 'aspirin-phospholipase.gro'))
top = traj.topology; xyz = traj.xyz * 10 - 40.2078
frames = np.arange(601)
def edr(k, w=5):
    v = np.array([E[k][5 * f] for f in frames]); return np.convolve(np.pad(v, (w // 2, w // 2), mode='edge'), np.ones(w) / w, mode='valid')

# geometry
pocket_res = [5, 9, 18] + list(range(21, 24)) + list(range(27, 32)) + [44, 47, 48, 63, 100]
pocket = np.array([a.index for a in top.atoms if a.residue.is_protein and a.residue.resSeq in pocket_res])
ain = np.array([a.index for a in top.atoms if a.residue.name == 'AIN'])
ca = [a.index for a in top.atoms if a.residue.name == 'CA'][0]
wO = np.array([a.index for a in top.atoms if a.residue.is_water and a.element.symbol == 'O'])
cog = xyz[:, ain, :].mean(1)
d_pocket = np.linalg.norm(cog - xyz[:, pocket, :].mean(1), axis=1)
d_ca = np.linalg.norm(cog - xyz[:, ca, :], axis=1)
nw = np.array([(np.linalg.norm(xyz[f][wO][:, None, :] - xyz[f][ain][None, :, :], axis=2).min(1) < 4.0).sum() for f in frames])

EV = [(330, 'touchdown'), (288, 'Coulomb snap'), (120, 'in the pocket'), (4, 'Ca$^{2+}$ handshake')]
PH = [(600, 350, '#dfe9f5', 'approach'), (350, 288, '#fbe5d0', 'dolly zoom'), (288, 120, '#e3f2e1', 'descent'), (120, 0, '#f3dff0', 'handshake')]

def decorate(ax, labels=False):
    ax.set_xlim(600, 0)
    for a, b, c, n in PH:
        ax.axvspan(a, b, color=c, alpha=0.6, lw=0)
    for f, n in EV:
        ax.axvline(f, color='k', lw=0.6, ls=':')
        if labels:
            ax.text(f, 1.02, n, transform=ax.get_xaxis_transform(), rotation=0, ha={330: 'right', 288: 'left'}.get(f, 'center'), va='bottom', fontsize=7)

# ---- Figure 1: energies and geometry along the (reversed) trajectory ----
fig, ax = plt.subplots(3, 1, figsize=(7.2, 6.6), sharex=True)
for k, lab, c in [('Coul-SR:Protein-AIN', 'protein', '#e41a1c'), ('Coul-SR:CA-AIN', 'Ca$^{2+}$', '#e6a100'), ('Coul-SR:AIN-SOL', 'water', '#1f77b4')]:
    ax[0].plot(frames, edr(k), c=c, lw=1.1, label=lab)
ax[0].set_ylabel('Coul-SR with aspirin\n(kJ/mol)'); ax[0].legend(ncol=3, frameon=False, loc='center left', bbox_to_anchor=(0.0, 0.55), fontsize=8)
for k, lab, c in [('LJ-SR:Protein-AIN', 'protein', '#e41a1c'), ('LJ-SR:CA-AIN', 'Ca$^{2+}$', '#e6a100'), ('LJ-SR:AIN-SOL', 'water', '#1f77b4')]:
    ax[1].plot(frames, edr(k), c=c, lw=1.1, label=lab)
ax[1].set_ylabel('LJ-SR with aspirin\n(kJ/mol)')
ax[2].plot(frames, d_pocket, c='#17becf', lw=1.1, label='pocket (centres)')
ax[2].plot(frames, d_ca, c='#e6a100', lw=1.1, label='Ca$^{2+}$ (centres)')
ax[2].set_ylabel('distance from aspirin (Å)'); ax[2].legend(frameon=False, loc='lower left', fontsize=8)
a2 = ax[2].twinx(); a2.plot(frames, nw, c='#1f77b4', lw=0.9, alpha=0.8); a2.set_ylabel('waters within 4 Å', color='#1f77b4'); a2.spines['right'].set_visible(True)
for i, a in enumerate(ax):
    decorate(a, labels=(i == 0))
ax[2].set_xlabel('trajectory frame (5 ps each); the movie runs from 600 to 0')
plt.tight_layout(); plt.savefig(f'{OUT}/fig1_energies.png', dpi=130); plt.close()

# ---- Figure 2: time schedule and the dolly zoom ----
rows = run_tool(story.keys, story.ANCH, (0, 60, 0.05))
ts = rows[:, 0]; fr = rows[:, 1]; fov = np.degrees(rows[:, 8]); dist = rows[:, 9]
fig, ax = plt.subplots(3, 1, figsize=(7.2, 6.8), sharex=True)
ax[0].plot(ts, fr, c='k', lw=1.3)
for t, f in story.PINS:
    ax[0].plot(t, f, 'o', c='#d62728', ms=4)
ax[0].set_ylabel('trajectory frame'); ax[0].invert_yaxis()
ax[1].semilogy(ts, dist, c='#1f77b4', label='camera distance $d$ (Å)'); ax[1].set_ylabel('distance (Å)', color='#1f77b4')
b = ax[1].twinx(); b.plot(ts, fov, c='#ff7f0e'); b.set_ylabel('field of view (deg)', color='#ff7f0e'); b.spines['right'].set_visible(True)
ax[2].plot(ts, dist * np.tan(np.radians(fov) / 2), c='k'); ax[2].set_ylabel('$d\\,\\tan(\\theta/2)$ (Å)\nhalf-height at the subject'); ax[2].set_xlabel('movie time (s)')
for a in ax:
    for t, n in [(4, ''), (16, ''), (25, ''), (30.5, ''), (44, ''), (57, '')]:
        a.axvline(t, color='gray', lw=0.5, ls=':')
    a.axvspan(16, 23, color='#fbe5d0', alpha=0.6, lw=0); a.axvspan(25, 30.5, color='#cfe3f7', alpha=0.7, lw=0)
ax[0].text(19.5, 590, 'dolly zoom', ha='center', fontsize=8); ax[0].text(27.75, 590, 'bullet time', ha='center', fontsize=8)
plt.tight_layout(); plt.savefig(f'{OUT}/fig2_schedule_dolly.png', dpi=130); plt.close()

# ---- Figure 3: bullet time ----
m = (ts >= 24.0) & (ts <= 31.0)
look = rows[m, 2:5]; eye = rows[m, 5:8]; rel = eye - look
az = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0])); el = np.degrees(np.arcsin(rel[:, 2] / np.linalg.norm(rel, axis=1)))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3))
sc = ax[0].scatter(rel[:, 0], rel[:, 1], c=ts[m], s=6, cmap='viridis'); ax[0].plot(0, 0, 'o', c='crimson', ms=7); ax[0].set_aspect('equal')
ax[0].set_xlabel('x (Å)'); ax[0].set_ylabel('y (Å)'); ax[0].set_title('eye around aspirin, seen from +z', fontsize=9)
cb = plt.colorbar(sc, ax=ax[0], fraction=0.046); cb.set_label('movie time (s)')
ax[1].plot(ts[m], np.degrees(az - az[0]), c='k'); ax[1].set_xlabel('movie time (s)'); ax[1].set_ylabel('turned (deg)')
a = ax[1].twinx(); a.plot(ts[m], fr[m], c='#d62728'); a.set_ylabel('trajectory frame', color='#d62728'); a.set_ylim(284, 300); a.invert_yaxis(); a.spines['right'].set_visible(True)
plt.tight_layout(); plt.savefig(f'{OUT}/fig3_bullet.png', dpi=130); plt.close()

# ---- Figure 4: representations and overlays on the timeline ----
T = 2.5

reps = [('protein-cartoon', [(0, 1), (15.5, 0), (31.5, 1)], '#4c72b0'), ('protein-vdw', [(0, 0), (15.5, 1), (24.5, 0)], '#dd8452'),
        ('protein-licorice', [(0, 0), (24.5, 1), (31.5, 0)], '#55a868'), ('pocket', [(0, 0), (27.5, 1), (45.5, 0)], '#c44e52'),
        ('ligand (ball-stick)', [(0, 1), (52, 0)], '#8172b3'), ('ligand-vdw', [(0, 0), (52, 1)], '#937860'),
        ('ion', [(0, 1)], '#e6a100'), ('water shell', [(0, 0), (8.5, 1), (45.5, 0)], '#1f77b4')]
tt = np.linspace(0, 60, 1201)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for i, (n, pr, c) in enumerate(reps):
    y = len(reps) - i
    vals = []
    prev = pr[0][1]
    cur = prev; last_t = None; target = prev
    for tv in tt:
        # latest key at or before tv
        k = None
        for kt, kv in pr:
            if kt <= tv: k = (kt, kv)
        if k is None or k == pr[0]:
            vals.append(pr[0][1]); continue
        idx = pr.index(k); frm = pr[idx - 1][1]
        # from wherever the previous change had got to (all gaps here exceed T)
        u = min((tv - k[0]) / T, 1.0); s = u * u * (3 - 2 * u)
        vals.append(frm + (k[1] - frm) * s)
    ax.fill_between(tt, y - 0.4, y - 0.4 + 0.8 * np.array(vals), color=c, alpha=0.85, lw=0)
ax.set_yticks([len(reps) - i for i in range(len(reps))]); ax.set_yticklabels([r[0] for r in reps], fontsize=8)
for t, n in [(4, 'fly in'), (16, 'dolly'), (25, 'bullet'), (30.5, ''), (44, 'pocket'), (57, 'finale')]:
    ax.axvline(t, color='gray', lw=0.5, ls=':')
ax.set_xlabel('movie time (s)'); ax.set_xlim(0, 60)
plt.tight_layout(); plt.savefig(f'{OUT}/fig4_representations.png', dpi=130); plt.close()

# ---- Figure 5: placement: the pocket opens to +z ----
f288 = 288
fig, ax = plt.subplots(1, 2, figsize=(7.4, 3.6), gridspec_kw={'width_ratios': [1.1, 1]})
P = PROT[f288]
ax[0].scatter(P[:, 0], P[:, 2], s=3, c='#b0b7c3', lw=0)
sc = ax[0].scatter(COG[::4, 0], COG[::4, 2], c=np.arange(0, 601, 4), s=10, cmap='plasma_r', zorder=3)
ax[0].plot(*[ION[0, 0]], ION[0, 2], '*', c='#e6a100', ms=12, mec='k', zorder=4)
ax[0].plot(xyz[0, pocket, 0].mean(), xyz[0, pocket, 2].mean(), 's', c='#17becf', ms=5, zorder=4)
for t_, c_ in [(16, '#ff7f0e'), (23, '#2ca02c')]:
    kk = [k for k in story.keys if abs(k.t - t_) < 1e-6][0]; e = kk.eye()
    ax[0].plot([e[0], kk.look[0]], [e[2], kk.look[2]], '-', c=c_, lw=1)
    ax[0].plot(e[0], e[2], 'v', c=c_, ms=6)
ax[0].set_xlim(-45, 45); ax[0].set_ylim(-25, 90); ax[0].set_aspect('equal')
ax[0].set_xlabel('x (Å)'); ax[0].set_ylabel('z (Å)'); cb = plt.colorbar(sc, ax=ax[0], fraction=0.04); cb.set_label('frame'); cb.ax.invert_yaxis()
# clear-view map for the bound state
tgt = 0.5 * (COG[0] + ION[0]); azs = np.arange(-180, 181, 10); els = np.arange(10, 81, 5)
Pm = PROT[0]; M = np.zeros((len(els), len(azs)))
for i, el_ in enumerate(els):
    for j, az_ in enumerate(azs):
        eye = tgt + 19 * eye_dir(az_, el_); v = tgt - eye; L_ = np.linalg.norm(v); v /= L_
        w = Pm - eye; s_ = w @ v; perp = np.linalg.norm(w - np.outer(s_, v), axis=1)
        M[i, j] = ((s_ > 0) & (s_ < L_ - 4) & (perp < 2.0)).sum()
im = ax[1].imshow(M, origin='lower', extent=[-185, 185, 7.5, 82.5], aspect='auto', cmap='Greys', vmax=15)
ax[1].set_xlabel('camera azimuth (deg)'); ax[1].set_ylabel('camera elevation (deg)'); ax[1].set_title('atoms blocking the view (bound state)', fontsize=8)
plt.colorbar(im, ax=ax[1], fraction=0.046)
plt.tight_layout(); plt.savefig(f'{OUT}/fig5_placement.png', dpi=130); plt.close()

# numbers quoted in the text
def mean_window(k, a, b): return np.mean([E[k][5 * f] for f in range(a, b)])
tot = lambda a, b: sum(mean_window(k, a, b) for k in ('Coul-SR:Protein-AIN', 'LJ-SR:Protein-AIN', 'Coul-SR:CA-AIN', 'LJ-SR:CA-AIN', 'Coul-SR:AIN-SOL', 'LJ-SR:AIN-SOL'))
print('window totals (frames): ', {w: round(tot(*w)) for w in [(0, 30), (30, 120), (120, 270), (270, 306), (306, 400), (400, 600)]})
print('water Coul+LJ', {w: round(mean_window('Coul-SR:AIN-SOL', *w) + mean_window('LJ-SR:AIN-SOL', *w)) for w in [(0, 30), (30, 120), (120, 270), (270, 306), (306, 400), (400, 600)]})
print('protein+Ca Coul+LJ', {w: round(mean_window('Coul-SR:Protein-AIN', *w) + mean_window('LJ-SR:Protein-AIN', *w) + mean_window('Coul-SR:CA-AIN', *w) + mean_window('LJ-SR:CA-AIN', *w)) for w in [(0, 30), (30, 120), (120, 270), (270, 306), (306, 400), (400, 600)]})
print('d_pocket frames 0,120,288,330,600:', d_pocket[[0, 120, 288, 330, 600]].round(1), 'd_ca', d_ca[[0, 120, 288, 330, 600]].round(1), 'nw', nw[[0, 120, 288, 330, 600]], 'mean nw 400-600', nw[400:].mean().round(1), 'frames 0-30', nw[:30].mean().round(1))

import os
import sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from story import *
times = [float(x) for x in sys.argv[1:]] or [2, 8, 12, 17, 20, 23, 25, 27.5, 33, 40, 47, 52, 57, 59]
rows = run_tool(keys, ANCH, (0, 60, 0.05))
def row(t): return rows[int(round(t / 0.05))]
n = len(times); cols = 4; rws = (n + cols - 1) // cols
fig, axs = plt.subplots(rws, cols, figsize=(cols * 4.8, rws * 2.7), facecolor='white')
axs = np.atleast_2d(axs)
for ax in axs.ravel(): ax.axis('off')
for ax, t in zip(axs.ravel(), times):
    r = row(t); fr = r[1]; look, eye, fov = r[2:5], r[5:8], r[8]
    fwd = look - eye; fwd /= np.linalg.norm(fwd); right = np.cross(fwd, [0, 0, 1]); right /= np.linalg.norm(right); up = np.cross(right, fwd)
    def proj(P):
        a = P - eye; z = a @ fwd; h = math.tan(fov / 2)
        return (a @ right) / (z * h) * 9 / 16, (a @ up) / (z * h), z
    px, py, pz = proj(at_frame(PROT, fr))
    m = pz > 1
    ax.scatter(px[m], py[m], s=np.clip(900 / (pz[m] * math.tan(fov/2)) ** 1.0 * 0.25, 0.3, 60), c=pz[m], cmap='Blues_r', alpha=0.6, vmin=pz[m].min(), vmax=pz[m].max() + 20)
    ax_, ay_, az_ = proj(at_frame(ASP, fr) if False else np.array([at_frame(ASP[:, i, :], fr) for i in range(20)]))
    sc = 1.0 / (np.mean(az_) * math.tan(fov / 2))
    ax.scatter(ax_, ay_, s=np.clip(60 * 25 * sc ** 1.0 * 1.0, 8, 400), c='crimson', zorder=5)
    ix, iy, iz = proj(at_frame(ION, fr)[None, :])
    ax.scatter(ix, iy, s=np.clip(2200 * sc, 10, 600), c='gold', edgecolors='k', zorder=6)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect('equal'); ax.set_xlim(-1.0, 1.0); ax.set_ylim(-1.0, 1.0)
    ax.add_patch(plt.Rectangle((-1, -1), 2, 2, fill=False, ec='k'))
    ax.set_title('t=%.1f f=%.0f fov=%.0f d=%.0f' % (t, fr, math.degrees(fov), r[9]), fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(WORK, 'preview.png'), dpi=80)

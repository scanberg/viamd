import os
import sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from movie_lib import *

DUR = 60.0
ANCH = (0.0, 600.0, 60.0, 0.0)
PINS = [(0, 600), (4, 600), (16, 350), (25, 288), (30.5, 288), (44, 120), (50, 30), (57, 0), (60, 0)]

_frames = frame_schedule(PINS, ANCH, (0, 60, 0.05))
def FR(t):
    return float(np.interp(t, _frames[:, 0], _frames[:, 1]))

def smooth(s):
    return s * s * (3 - 2 * s)

def L(t, sigma=3.0, off=(0, 0, 0)):
    return asp_s(FR(t), sigma) + np.array(off, float)

def fov_for(d, width):
    return math.degrees(2 * math.atan(width / 2 / d))

def mid(t):
    f = FR(t)
    return 0.5 * (at_frame(COG, f) + at_frame(ION, f))

keys = []
add = keys.append

add(Key(0.0, (0, 0, 6), -60, 18, 125, 34, name='Intro'))
add(Key(4.0, (-5, -6, 14), -52, 20, 88, 34, name='Fly in'))
add(Key(8.0, 0.65 * L(8) + 0.35 * np.array([0, 0, 8]), -40, 24, 70, 34, name='Approach'))
add(Key(12.0, 0.85 * L(12) + 0.15 * np.array([0, 0, 8]), -25, 32, 48, 32))
W, D0, D1 = 22.0, 14.0, 80.0
for i, t in enumerate(np.arange(16.0, 23.01, 1.0)):
    s = smooth((t - 16.0) / 7.0)
    d = D0 * (D1 / D0) ** s
    add(Key(float(t), L(t, 2.0), -90, 55, d, fov_for(d, W), name='Touchdown' if i == 0 else ''))
look_b = COG[288].astype(float)
add(Key(25.0, look_b, -90, 28, 22, 40, roll=7, ease=1, name='Bullet time'))
add(Key(30.5, look_b, -90, 36, 22, 40, roll=7, ease=1, spin=1, spin_axis=0, const=0))
add(Key(33.0, L(33), -60, 55, 18, 36, name='Descent'))
add(Key(38.0, L(38), -30, 60, 16, 34))
add(Key(44.0, L(44), 20, 60, 16, 32, name='Pocket'))
add(Key(47.0, mid(47), 150, 55, 20, 32, name='Handshake'))
add(Key(50.0, mid(50), 120, 62, 19, 32))
add(Key(54.0, mid(54), 90, 68, 19, 32))
add(Key(57.0, mid(57), 60, 62, 19, 32))
add(Key(60.0, (0, 0, 3), 70, 50, 105, 36, name='Finale'))
pin = dict(PINS)
for k in keys:
    if k.t in pin:
        k.frame = pin[k.t]

def diagnostics():
    rows = run_tool(keys, ANCH, (0, 60, 0.25))
    print(' t     frame  eyeClear(A)  LOSblock  asp_screen(x,y in half-heights)  size(d*tan)')
    for r in rows[::2]:
        t, fr = r[0], r[1]
        look, eye, fov, d = r[2:5], r[5:8], r[8], r[9]
        P = at_frame(PROT, fr)
        clear = np.linalg.norm(P - eye, axis=1).min()
        fwd = (look - eye); fwd /= np.linalg.norm(fwd)
        right = np.cross(fwd, [0, 0, 1]); right /= np.linalg.norm(right); up = np.cross(right, fwd)
        a = at_frame(COG, fr) - eye
        z = a @ fwd
        sx = (a @ right) / (z * math.tan(fov / 2)); sy = (a @ up) / (z * math.tan(fov / 2))
        v = at_frame(COG, fr) - eye; L_ = np.linalg.norm(v); v /= L_
        w = P - eye; s = w @ v
        perp = np.linalg.norm(w - np.outer(s, v), axis=1)
        block = ((s > 0) & (s < L_ - 4) & (perp < 2.0)).sum() if d < 40 else 0
        print('%5.2f %6.1f  %8.1f  %6d    (%5.2f %5.2f)   %6.2f' % (t, fr, clear, block, sx, sy, d * math.tan(fov / 2)))

if __name__ == '__main__':
    diagnostics()

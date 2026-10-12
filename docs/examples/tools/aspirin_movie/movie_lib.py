import numpy as np, subprocess, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..', '..'))
WORK = os.environ.get('ASPIRIN_WORK', os.path.join(HERE, 'work'))
DATA = os.environ.get('ASPIRIN_DATA', os.path.expanduser('~/Desktop/bioexcel/aspirin'))

D = np.load(os.path.join(WORK, 'asp_tracks.npz'))
PROT, ASP, ION = D['prot'], D['asp'], D['ion']
COG = ASP.mean(1)
NF = PROT.shape[0]
TOOL = os.path.join(WORK, 'camtool')


def asp_s(frame, sigma=4.0):
    """Smoothed aspirin center of geometry at a (fractional) frame."""
    fs = np.arange(max(0, int(frame - 4 * sigma)), min(NF, int(frame + 4 * sigma) + 2))
    w = np.exp(-0.5 * ((fs - frame) / sigma) ** 2)
    return (COG[fs] * w[:, None]).sum(0) / w.sum()


def at_frame(arr, frame):
    f = min(max(frame, 0.0), NF - 1.0)
    f0 = int(math.floor(f)); f1 = min(f0 + 1, NF - 1); w = f - f0
    return arr[f0] * (1 - w) + arr[f1] * w


def quat_from_R(R):
    m = R
    tr = m[0, 0] + m[1, 1] + m[2, 2]
    if tr > 0:
        s = math.sqrt(tr + 1) * 2
        w = 0.25 * s; x = (m[2, 1] - m[1, 2]) / s; y = (m[0, 2] - m[2, 0]) / s; z = (m[1, 0] - m[0, 1]) / s
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        s = math.sqrt(1 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
        w = (m[2, 1] - m[1, 2]) / s; x = 0.25 * s; y = (m[0, 1] + m[1, 0]) / s; z = (m[0, 2] + m[2, 0]) / s
    elif m[1, 1] > m[2, 2]:
        s = math.sqrt(1 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
        w = (m[0, 2] - m[2, 0]) / s; x = (m[0, 1] + m[1, 0]) / s; y = 0.25 * s; z = (m[1, 2] + m[2, 1]) / s
    else:
        s = math.sqrt(1 + m[2, 2] - m[0, 0] - m[1, 1]) * 2
        w = (m[1, 0] - m[0, 1]) / s; x = (m[0, 2] + m[2, 0]) / s; y = (m[1, 2] + m[2, 1]) / s; z = 0.25 * s
    q = np.array([x, y, z, w]); return q / np.linalg.norm(q)


def eye_dir(az, el):
    a, e = math.radians(az), math.radians(el)
    return np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])


def orientation(az, el, up=(0, 0, 1)):
    back = eye_dir(az, el)          # from the look-at point to the eye
    f = -back
    r = np.cross(f, np.array(up, float)); r /= np.linalg.norm(r)
    u = np.cross(r, f)
    return quat_from_R(np.column_stack([r, u, back]))


class Key:
    def __init__(self, t, look, az, el, d, fov_deg, roll=0.0, frame=None, ease=0, spin=0, spin_axis=0, const=0, name=''):
        self.t = t; self.look = np.array(look, float); self.az = az; self.el = el; self.d = d
        self.fov = math.radians(fov_deg); self.roll = roll; self.frame = frame; self.ease = ease
        self.spin = spin; self.spin_axis = spin_axis; self.const = const; self.name = name

    def eye(self):
        return self.look + self.d * eye_dir(self.az, self.el)

    def quat(self):
        return orientation(self.az, self.el)

    def tool_line(self):
        e = self.eye(); q = self.quat()
        return '%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %d %.6f %d %d %d %d 0 0 0 0 -1 %.6f' % (
            self.t, self.fov, self.d, e[0], e[1], e[2], q[0], q[1], q[2], q[3],
            1 if self.frame is not None else 0, self.frame if self.frame is not None else 0.0,
            self.spin, self.spin_axis, self.const, self.ease, math.radians(self.roll))


def run_tool(keys, anchors, times, upright=True):
    txt = [str(len(keys))] + [k.tool_line() for k in keys]
    txt.append('%f %f %f %f %d 0 0 1' % (anchors[0], anchors[1], anchors[2], anchors[3], 1 if upright else 0))
    txt.append('1 %d' % NF)
    txt.append('0 0 0' + ' 0 0 0' * (NF - 1))
    t0, t1, dt = times
    txt.append('%f %f %f' % (t0, t1, dt))
    out = subprocess.run([TOOL], input='\n'.join(txt), capture_output=True, text=True).stdout
    rows = np.array([[float(x) for x in l.split()] for l in out.strip().split('\n')])
    return rows  # t frame look(3) eye(3) fov dist


def frame_schedule(pins, anchors, times):
    keys = [Key(t, (0, 0, 0), 0, 0, 10, 30, frame=f) for t, f in pins]
    return run_tool(keys, anchors, times, upright=False)


def frame_at(pins, anchors, t):
    return float(frame_schedule(pins, anchors, (t, t, 1.0))[0, 1])

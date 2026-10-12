import os
import sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import story
from story import keys, PINS, DUR, ANCH

from movie_lib import ROOT
SRC = os.path.join(ROOT, 'docs/examples/aspirin_phospholipase_movie.via')
OUT = os.path.join(ROOT, 'docs/examples/aspirin_binding_movie.via')

ex = open(SRC, encoding='utf-8').read()
head = ex[:ex.index('[Animation]')]                       # banner + [Files]
render = ex[ex.index('[RenderSettings]'):ex.index('[Camera]')]
tail = ex[ex.index('[ActiveSelection]'):]
dist_view = ex[ex.index('[TimelineView]'):ex.index('[RenderSettings]')]

def f(x):
    return ('%.6f' % x).rstrip('0').rstrip('.') if abs(x) < 1e7 else '%g' % x

def vec(*v):
    return ','.join(f(x) for x in v)

out = [head.rstrip('\n') + '\n\n']

out.append('''[Animation]
Frame=600
Fps=2.94000006
Tension=0
Interpolation=2

''')

# ---- plots of the Timelines and Distributions windows ----
SUBPLOTS = {1: 'Coulomb (SR), aspirin with ...', 2: 'Lennard-Jones (SR), aspirin with ...', 3: 'Distances (A)', 4: 'Waters within 4 A of aspirin'}
out.append('[Timeline]\nNumSubplots=4\n\n')
for i in range(10):
    out.append('[TimelineSubplot]\nSubplot=%d\nId=%d\nName=%s\n\n' % (i, i + 1, SUBPLOTS.get(i + 1, '')))

def series(kind, sub, ident, color):
    return ('[%sSeries]\nSubplot=%d\nSource=script\nVariant=values\nPath=script/%s\nColor=%s,1\nPlotType=line\nMarker=1\nMarkerSize=1\nBarWidth=1\n'
            'NumBins=128\nUseColormap=0\nColormap=5\nColormapAlpha=1\n\n') % (kind, sub, ident, color)

RED, GOLD, BLUE = '0.890196145,0.101960793,0.10980393', '0.96,0.65,0.05', '0.121568635,0.470588267,0.70588237'
CYAN, PURPLE = '0.258823544,0.807843208,0.890196145', '0.55,0.2,0.85'
for i, (ident, col) in enumerate([('Coul_protein', RED), ('Coul_Ca', GOLD), ('Coul_water', BLUE)]):
    out.append(series('Timeline', 0, ident, col))
for ident, col in [('LJ_protein', RED), ('LJ_Ca', GOLD), ('LJ_water', BLUE)]:
    out.append(series('Timeline', 1, ident, col))
out.append(series('Timeline', 2, 'd_pocket', CYAN))
out.append(series('Timeline', 2, 'd_Ca', GOLD))
out.append(series('Timeline', 3, 'n_water', BLUE))

out.append('[Distributions]\nNumSubplots=1\n\n')
for i in range(10):
    out.append('[DistributionSubplot]\nSubplot=%d\nId=%d\nName=%s\n\n' % (i, 11 + i, 'Carboxylate twist (deg)' if i == 0 else ''))
out.append(series('Distribution', 0, 'dih', '0.698039234,0.874509871,0.541176498'))

out.append(dist_view)
out.append(render)

# ---- camera at the start ----
k0 = keys[0]
e0, q0 = k0.eye(), k0.quat()
out.append('[Camera]\nPosition=%s\nOrientation=%s\nDistance=%s\nMode=0\nFovY=%s\n\n' % (vec(*e0), vec(*q0), f(k0.d), f(k0.fov)))

# ---- movie ----
EASE_HOLD = 3
lanes = 1 | 2 | 4 | 32 | 64 | 128 | 256 | 512 | 1024 | 2048
m = ['[Movie]', 'Resolution=1', 'ResX=1920', 'ResY=1080', 'Fps=24', 'StartFrame=600', 'EndFrame=0', 'Timeline=2', 'Duration=%d' % DUR,
     'TrajectoryBegin=0', 'TrajectoryEnd=%d' % DUR, 'Playhead=0', 'FilenamePrefix=aspirin_binding', 'AnimateCamera=1', 'ShowPath=0',
     'PathOptions=15', 'Output=1', 'Crf=18', 'ResScale=100', 'AaSamples=0', 'SaveCopy=1', 'RepTransition=2.5', 'SnapFrames=1',
     'Lanes=%d' % lanes, 'LaneHeight=200', 'FitLanes=1', 'RepEqualRows=1', 'Overlays=1', 'RenderRange=0,0,%d' % DUR, 'Loop=0',
     'KeepUpright=1', 'UpAxis=2', 'AnimateParams=1']
for k in keys:
    e, q = k.eye(), k.quat()
    uf = 1 if k.frame is not None else 0
    vals = [k.t, k.fov, k.d, e[0], e[1], e[2], q[0], q[1], q[2], q[3], uf, k.frame if uf else 0, k.spin, k.spin_axis, k.const, k.ease, 0, 0, 0, 0, -1]
    m.append('KeyframeV3=' + vec(*vals))
    if k.name:
        m.append('KeyframeName=' + k.name)
    if k.roll:
        m.append('KeyframeRoll=' + f(k.roll))

# Look parameters: 0 background color, 1 background intensity, 9 DOF aperture (%).
P = []
def param(p, t, v, ease=0):
    v = list(v) + [0] * (3 - len(v))
    P.append('ParamKey=' + vec(p, t, *v, ease))
for t, v in [(0, 0), (14.5, 0), (19, 1.5), (23, 1.5), (25, 2), (30.5, 2), (33, 0)]:
    param(9, t, [v])
for t, c in [(0, (1, 1, 1)), (24.5, (1, 1, 1)), (26.5, (0.55, 1.0, 0.7)), (30, (0.55, 1.0, 0.7)), (33, (1, 1, 1))]:
    param(0, t, c)
for t, v in [(0, 24), (24.5, 24), (26.5, 14), (30, 14), (33, 24)]:
    param(1, t, [v])
m += P

# Representations: id -> (name, ...). Visible keys use hold; a change starts at its key and takes RepTransition seconds.
R = []
def rep_key(rep, prop, t, v, ease=0, v2=0, v3=0):
    R.append('RepKey=' + vec(rep, prop, t, v, ease, v2, v3))
def visible(rep, pairs):
    for t, v in pairs:
        rep_key(rep, 0, t, v, EASE_HOLD)
visible(1, [(0, 1), (15.5, 0), (31.5, 1)])        # protein-cartoon: leaves for the surface and comes back after the bullet time
visible(2, [(0, 0), (15.5, 1), (24.5, 0)])        # protein-vdw: the backdrop of the dolly zoom
visible(3, [(0, 0), (24.5, 1), (31.5, 0)])        # protein-licorice: the frozen world
visible(4, [(0, 0), (27.5, 1), (45.5, 0)])        # pocket
visible(5, [(0, 1), (52, 0)])                     # ligand as ball and stick, then as spheres
visible(6, [(0, 0), (52, 1)])
visible(7, [(0, 1)])
visible(8, [(0, 0), (8.5, 1), (45.5, 0)])         # the waters around aspirin, which leave when it docks
for t, v in [(0, 2.4), (16, 2.4), (25, 1.6), (33, 1.25), (60, 1.25)]:
    rep_key(5, 1, t, v)                           # aspirin ball scale
for t, v in [(0, 1.1), (52, 1.1), (54.8, 1.6), (58, 1.1)]:
    rep_key(7, 1, t, v)                           # the calcium ion swells at the strongest attraction
for t, v in [(0, 1.0), (17, 1.0), (19.5, 0.35), (23, 0.35), (24.5, 1.0)]:
    rep_key(2, 5, t, v)                           # the surface loses its colour behind the dolly zoom
m += R
out.append('\n'.join(m) + '\n\n')

# ---- overlays ----
def overlay(typ, begin, end, fin, fout, anchor, size, color=(0, 0, 0, 1), bg=(0, 0, 0, 0), text='', extra=''):
    return ('[MovieOverlay]\nType=%d\nEnabled=1\nRange=%s\nAnchor=%d\nSize=%s\nSizeUnit=0\nColor=%s\nBackground=%s\nLength=0\nText=%s\n%s\n'
            % (typ, vec(begin, end, fin, fout), anchor, f(size), vec(*color), vec(*bg), text, extra))
TL, TC, TR, ML, C, MR, BL, BC, BR = range(9)
O = []
O.append(overlay(3, 0, DUR, 0, 0, TL, 0.08, (1, 1, 1, 1)))                                  # logo
O.append(overlay(0, 0.4, 5.0, 1.0, 1.0, TC, 0.07, text='Aspirin finds its pocket'))
O.append(overlay(0, 1.2, 5.0, 1.0, 1.0, BC, 0.032, text='Phospholipase A2 with Ca2+. 3 ns of molecular dynamics, played in reverse'))
for b, e, txt in [(5.5, 13.5, '1 | The approach'), (16.0, 23.5, '2 | First contact: the dolly zoom'),
                  (25.0, 30.5, '3 | Bullet time: the Coulomb snap'), (33.0, 43.5, '4 | Into the pocket'),
                  (47.0, 56.5, '5 | The calcium handshake')]:
    O.append(overlay(0, b, e, 0.7, 0.7, TC, 0.05, text=txt))
O.append(overlay(0, 57.6, DUR, 0.8, 0.0, TC, 0.07, text='Aspirin finds its pocket'))
O.append(overlay(1, 5.0, DUR, 0.7, 0.0, BL, 0.045))                                         # time stamp
O.append(overlay(2, 8.0, 57.0, 1.0, 1.0, BR, 0.05))                                          # scale bar
O.append(overlay(5, 5.0, DUR, 0.7, 0.0, BC, 0.035, (0.649350643, 0, 1, 1), (0, 0, 0, 0.258823544), extra='Width=0.45\nTimeBarLabels=3\n'))
plot_common = 'PlotAxis=0\nPlotFlags=15\nFontPoints=21\nLinePoints=0\nPalette=0\n'
O.append(overlay(6, 6.0, 58.0, 1.0, 1.0, ML, 0.62, (1, 1, 1, 1), (0, 0, 0, 0.5),
                 extra='Width=0.30\n' + plot_common + 'Panel=0,1\nPanelTime=6,0\nPanel=0,2\nPanelTime=9.5,0\n'))
O.append(overlay(6, 11.0, 58.0, 1.0, 1.0, MR, 0.46, (1, 1, 1, 1), (0, 0, 0, 0.5),
                 extra='Width=0.28\n' + plot_common + 'Panel=0,3\nPanelTime=11,0\nPanel=0,4\nPanelTime=15,0\n'))
O.append(overlay(7, 40.0, 58.0, 1.0, 1.0, TR, 0.25, (1, 1, 1, 1), (0, 0, 0, 0.5),
                 extra='Width=0.20\n' + plot_common + 'NumBins=32\nPanel=1,11\nPanelTitle=Carboxylate twist\n'))
for b, e, ident in [(9.0, 14.5, 'n_water'), (17.0, 23.5, 'd_pocket'), (36.0, 42.0, 'dih'), (47.5, 56.0, 'd_Ca')]:
    O.append(overlay(8, b, e, 0.8, 0.8, BL, 0.05, (1, 1, 1, 1), text=ident))
out.append(''.join(O))

# ---- markers ----
for t, lab in [(18.0, 'Touchdown'), (25.0, 'Coulomb snap'), (44.0, 'In the pocket'), (54.8, 'Ca2+ handshake')]:
    out.append('[MovieMarker]\nTime=%s\nLabel=%s\nSubplot=0\n\n' % (f(t), lab))

# ---- operations, representations, script ----
ops = ex[ex.index('[Operations]'):ex.index('[Representation]')]
out.append(ops)

SS = ('SecondaryStructureColorUnknown=0.5,0.5,0.5,1\nSecondaryStructureColorCoil=0.860000014,0.860000014,0.860000014,1\n'
      'SecondaryStructureColorHelix=0.119999997,0.860000014,0.119999997,1\nSecondaryStructureColorSheet=0.119999997,0.119999997,0.860000014,1\n')
def rep(i, name, filt, enabled, typ, cm, param, base=(1, 1, 1, 1), tint=(1, 0, 0, 1), tscale=0, bond=0, dyn=0):
    return ('[Representation]\nId=%d\nName=%s\nFilter=%s\nEnabled=%d\nType=%d\nColorMapping=%d\nBaseColor=%s\nSaturation=1\nTintColor=%s\nTintScale=%s\n%s'
            'BondColor=%d\nBondSharpness=0.5\nBondBaseColor=1,1,1,1\nParam=%s\nDynamicEval=%d\n\n'
            % (i, name, filt, enabled, typ, cm, vec(*base), vec(*tint), f(tscale), SS, bond, vec(*param), dyn))
POCKET = 'residue({5,9,18,21:23,27:31,44,47:48,63,100});'
LIG = 'not (protein or nucleic or water or ion)'
# types: 0 spacefill, 1 licorice, 2 ball and stick, 4 cartoon; color mappings: 0 uniform, 1 type, 5 residue index, 8 secondary structure
out.append(rep(1, 'protein-cartoon', 'protein', 1, 4, 8, (1, 1, 1, 1)))
out.append(rep(2, 'protein-vdw', 'protein', 0, 0, 5, (0.85, 1, 1, 1)))
out.append(rep(3, 'protein-licorice', 'protein', 0, 1, 1, (1.0, 1, 1, 1)))
out.append(rep(4, 'pocket', POCKET, 0, 0, 5, (0.7, 1, 1, 1)))
out.append(rep(5, 'ligand', LIG, 1, 2, 1, (2.4, 1.195, 1, 1), bond=1))
out.append(rep(6, 'ligand-vdw', LIG, 0, 0, 1, (0.85, 1.195, 1, 1), bond=1))
out.append(rep(7, 'ion', 'ion', 1, 0, 0, (1.1, 1, 1, 1), base=(1.0, 0.78, 0.1, 1)))
out.append(rep(8, 'water', 'residue(resname("SOL") and within(4.0,resname("AIN")))', 0, 2, 1, (0.8, 0.538, 1, 1), dyn=1))

out.append('''[Script]
Text=\"\"\"
sel1 = residue({5,9,18,21:23,27:31,44,47:48,63,100});
d_pocket = distance(sel1,resname("AIN"));
d_Ca = distance(resname("CA"),resname("AIN"));
n_water = count(resname("SOL") and within(4.0,resname("AIN")),'residue');
dih = dihedral(1198, 1207, 1208, 1209);

Coul_protein = attr("edr/coul_sr_protein_ain");
Coul_Ca = attr("edr/coul_sr_ca_ain");
Coul_water = attr("edr/coul_sr_ain_sol");
LJ_protein = attr("edr/lj_sr_protein_ain");
LJ_Ca = attr("edr/lj_sr_ca_ain");
LJ_water = attr("edr/lj_sr_ain_sol");
\"\"\"

''')
out.append(tail)
open(OUT, 'w', encoding='utf-8').write(''.join(out))
print('wrote', OUT, len(''.join(out).splitlines()), 'lines')

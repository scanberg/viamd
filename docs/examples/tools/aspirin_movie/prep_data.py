"""Builds the two data files the movie scripts read, from the aspirin-phospholipase .gro, .xtc and .edr.

  python3 prep_data.py        (needs numpy, mdtraj, pyedr; ASPIRIN_DATA is the folder of the three files)

Order of the whole pipeline: prep_data.py, build.sh, story.py (diagnostics), gen_movie.py (writes the .via),
roundtrip.py (reads it back), make_figs.py (figures of the article). The files go to ASPIRIN_WORK (default: ./work).
"""
import os, pickle
import numpy as np
import mdtraj as md
import pyedr

HERE = os.path.dirname(os.path.abspath(__file__))
WORK = os.environ.get('ASPIRIN_WORK', os.path.join(HERE, 'work'))
DATA = os.environ.get('ASPIRIN_DATA', os.path.expanduser('~/Desktop/bioexcel/aspirin'))
os.makedirs(WORK, exist_ok=True)

traj = md.load(os.path.join(DATA, 'aspirin-phospholipase.xtc'), top=os.path.join(DATA, 'aspirin-phospholipase.gro'))
top = traj.topology
BOX_CENTER = 40.2078   # Å; VIAMD puts the centre of the cell at the origin

# View space: positions in Å with the cell centre at the origin (atoms 1..1189 protein, 1190 Ca2+, 1191..1210 aspirin)
xyz = traj.xyz * 10 - BOX_CENTER
np.savez(os.path.join(WORK, 'asp_tracks.npz'), prot=xyz[:, :1189].astype(np.float32), asp=xyz[:, 1190:1210].astype(np.float32),
         ion=xyz[:, 1189].astype(np.float32), time=traj.time)

# Distances in Å (minimum atom distances) and the energy file
ain = top.select('resname AIN'); ca = top.select('resname CA'); prot = top.select('protein')
pocket_res = [5, 9, 18] + list(range(21, 24)) + list(range(27, 32)) + [44, 47, 48, 63, 100]
pocket = np.array([a.index for a in top.atoms if a.residue.is_protein and a.residue.resSeq in pocket_res])
raw = traj.xyz * 10
def mind(a, b):
    return np.array([np.linalg.norm(raw[f][a][:, None, :] - raw[f][b][None, :, :], axis=2).min() for f in range(traj.n_frames)])
pickle.dump(dict(time=traj.time, d_pocket=mind(ain, pocket), d_ca=mind(ain, ca), d_prot=mind(ain, prot), d_caprot=mind(ca, prot),
                 cog=raw[:, ain, :].mean(1), edr=pyedr.edr_to_dict(os.path.join(DATA, 'aspirin-phospholipase.edr'))),
            open(os.path.join(WORK, 'asp_analysis.pkl'), 'wb'))
print('wrote', WORK)

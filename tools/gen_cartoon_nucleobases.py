#!/usr/bin/env python3
"""Generates src/cartoon_nucleobases.inl: the rings of the bases of the Martini nucleic acids, as weighted sums of their
beads, for the cartoon (cartoon_glyphs.cpp).

A Martini base is a rigid body of three to six beads in the plane of its rings, which holds no more of its geometry than
where it is. The cartoon draws it with the rings of an ideal base, placed where its beads are. Each ring atom is in the
plane of the beads, and so a fixed affine combination of three of them (weights that sum to 1): the one which gives the
ring atom from the beads of the ideal base mapped as the force field maps it. A rigid motion of the beads moves the rings
with them, and the combination is a derived point of md_gl (md_gl_mol_set_derived), evaluated on the GPU in every frame.
The weights are those of the real beads, the least squares ones where there are more than three (the four of a purine of
Martini 2, whose quadrangle of beads is not quite rigid): virtual sites follow from them, and a trajectory may lack them.

The ideal bases are the mean of those of an atomistic DNA (the nucleosome of nucl-dna.gro in mdlib's test data, B-DNA
of 87 A, 60 G, 60 C and 87 T, each within 0.03 A of the mean); uracil is the thymine without its methyl. The mappings:
  Martini 2 DNA and RNA   martinize-dna.py (Uusitalo et al. 2015) and martinize-nucleotide.py (2017), the bases of RNA
                          as those of DNA. The distances between the bead centres of the ideal bases are the constraint
                          lengths of the topologies (polyply's DNA_M2) within 0.01 nm.
  Martini 3 RNA           reForge's martinize_rna.py (Martini30RNA), virtual sites SC2 and SC5 of A, SC3 of C and U, SC2,
                          SC4 and SC6 of G.
The bases of a model are told apart by the number of their beads (what names a residue is its base: DA, A, RA, ...):
3 and 4 for the pyrimidines and purines of Martini 2, 4, 5 and 6 for C and U, A, G of Martini 3.

Checked against the frames of martini/martini2_dna, martini2_rna and martini3_rna of the test data (--check): the rings
of every base of them lie within 0.02 A (Martini 3) and 0.1 A (Martini 2, whose purines flex) of an ideal ring rigidly
placed (rmsd), and every ring atom of the atomistic DNA mapped to the beads of each model comes back within 0.09 A.

    python3 tools/gen_cartoon_nucleobases.py ext/mdlib/test_data [--check]
"""

import argparse, os, sys
import numpy as np

# The atoms of each base which the mappings use, and the ring atoms
BASE_ATOMS = {
    'A': "N9 C8 H8 N7 C5 C6 N6 H61 H62 N1 C2 H2 N3 C4".split(),
    'G': "N9 C8 H8 N7 C5 C6 O6 N1 H1 C2 N2 H21 H22 N3 C4".split(),
    'C': "N1 C6 H6 C5 H5 C4 N4 H41 H42 N3 C2 O2".split(),
    'T': "N1 C6 H6 C5 C7 H71 H72 H73 C4 O4 N3 H3 C2 O2".split(),
}
RING6 = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6']
PURINE_EXTRA = ['N7', 'C8', 'N9']   # the five ring C4 C5 N7 C8 N9

M2 = {
    'A': {'SC1': ['N9', 'C4'], 'SC2': ['C2', 'N3'], 'SC3': ['C6', 'N6', 'N1'], 'SC4': ['C8', 'N7', 'C5']},
    'G': {'SC1': ['N9', 'C4'], 'SC2': ['C2', 'N2', 'N3'], 'SC3': ['C6', 'O6', 'N1'], 'SC4': ['C8', 'N7', 'C5']},
    'C': {'SC1': ['N1', 'C6'], 'SC2': ['N3', 'C2', 'O2'], 'SC3': ['C5', 'C4', 'N4']},
    'T': {'SC1': ['N1', 'C6'], 'SC2': ['N3', 'C2', 'O2'], 'SC3': ['C5', 'C4', 'O4', 'C7']},
    'U': {'SC1': ['N1', 'C6'], 'SC2': ['N3', 'C2', 'O2'], 'SC3': ['C5', 'C4', 'O4']},
}
M3 = {
    'A': {'SC1': ['N9', 'C8', 'H8'], 'SC2': ['N3', 'C4'], 'SC3': ['N1', 'C2', 'H2'], 'SC4': ['N6', 'C6', 'H61', 'H62'], 'SC5': ['N7', 'C5']},
    'C': {'SC1': ['N1', 'C5', 'C6'], 'SC2': ['C2', 'O2'], 'SC3': ['N3'], 'SC4': ['N4', 'C4', 'H41', 'H42']},
    'G': {'SC1': ['C8', 'H8', 'N9'], 'SC2': ['C4', 'N3'], 'SC3': ['C2', 'N2', 'H21', 'H22'], 'SC4': ['N1'], 'SC5': ['C6', 'O6'], 'SC6': ['C5', 'N7']},
    'U': {'SC1': ['N1', 'C5', 'C6'], 'SC2': ['C2', 'O2'], 'SC3': ['N3'], 'SC4': ['C4', 'O4']},
}
M3_VSITES = {'A': {'SC2', 'SC5'}, 'C': {'SC3'}, 'G': {'SC2', 'SC4', 'SC6'}, 'U': {'SC3'}}

MODELS = [('Martini 2', M2, {}), ('Martini 3 RNA', M3, M3_VSITES)]


def read_gro(path):
    """[(resid, resname, {atom name: xyz in A})] in the order of the file"""
    lines = open(path).read().splitlines()
    n = int(lines[1])
    out = []
    for l in lines[2:2 + n]:
        key = (int(l[0:5]), l[5:10].strip())
        if not out or out[-1][0] != key:
            out.append((key, {}))
        out[-1][1][l[10:15].strip()] = np.array([float(l[20:28]), float(l[28:36]), float(l[36:44])]) * 10.0
    return [(k[0], k[1], d) for k, d in out]


def kabsch_fit(P, Q):
    """P rigidly moved onto Q (rows are points)"""
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    return (Vt.T @ D @ U.T @ (P - pc).T).T + qc


def rmsd(A, B):
    return float(np.sqrt(((A - B) ** 2).sum(1).mean()))


def ideal_bases(nucl_dna):
    sets = {b: [] for b in BASE_ATOMS}
    for _, rn, d in read_gro(nucl_dna):
        if rn in ('DA', 'DC', 'DG', 'DT') and all(a in d for a in BASE_ATOMS[rn[1]]):
            sets[rn[1]].append(np.array([d[a] for a in BASE_ATOMS[rn[1]]]))
    ideal = {}
    for b, s in sets.items():
        if not s:
            sys.exit(f'no complete {b} in {nucl_dna}')
        mean = s[0]
        for _ in range(10):
            mean = np.mean([kabsch_fit(x, mean) for x in s], 0)
        ideal[b] = dict(zip(BASE_ATOMS[b], mean))
    t = ideal['T']
    ideal['U'] = {k: v for k, v in t.items() if k not in ('C7', 'H71', 'H72', 'H73')}
    return ideal


def corners(b):
    return RING6 + PURINE_EXTRA if b in 'AG' else RING6


def weights(b, ideal_base, mapping, vsites):
    """The real beads, and the weights [corner][bead] which give each ring atom of base b from them"""
    names = sorted(mapping)
    centre = {k: np.mean([ideal_base[a] for a in mapping[k]], 0) for k in names}
    real = [k for k in names if k not in vsites]
    B = np.array([centre[k] for k in real])
    o = B.mean(0)
    _, _, vt = np.linalg.svd(np.array([centre[k] for k in names]) - np.array([centre[k] for k in names]).mean(0))
    e1, e2 = vt[0], vt[1]
    R = np.array([ideal_base[a] for a in corners(b)])
    A = np.vstack([(B - o) @ e1, (B - o) @ e2, np.ones(len(real))])
    X = np.vstack([(R - o) @ e1, (R - o) @ e2, np.ones(len(R))])
    W = (np.linalg.pinv(A) @ X).T
    return real, W


def table(ideal):
    out = []
    for model, mapping, vsites in MODELS:
        for b in sorted(mapping):
            real, W = weights(b, ideal[b], mapping[b], vsites.get(b, set()))
            out.append((model, b, len(mapping[b]), real, W))
    return out


def letter(rn):
    rn = rn.rstrip('35')
    if len(rn) == 1 and rn in 'ACGTU':
        return rn
    if len(rn) == 2 and rn[0] in 'DR' and rn[1] in 'ACGTU':
        return rn[1]
    return None


def check(test_data, ideal, entries):
    """The rings of the bases of the CG frames against an ideal ring placed rigidly; the atomistic DNA mapped and back"""
    ok = True
    by_key = {(b, n): (real, W) for _, b, n, real, W in entries}
    for name in ('martini2_dna', 'martini2_rna', 'martini3_rna'):
        worst = {}
        for _, rn, d in read_gro(os.path.join(test_data, 'martini', name + '.gro')):
            b = letter(rn)
            if not b:
                continue
            n = sum(1 for k in d if k.startswith('SC'))
            real, W = by_key[(b, n)]
            rec = W @ np.array([d[k] for k in real])
            ref = np.array([ideal[b][a] for a in corners(b)])
            worst[b] = max(worst.get(b, 0.0), rmsd(kabsch_fit(ref, rec), rec))
        print(f'{name}: the rings of each base within ' + ', '.join(f'{b} {v:.3f}' for b, v in sorted(worst.items())) + ' A of an ideal ring')
        ok &= all(v < 0.15 for v in worst.values())
    for model, mapping, _ in MODELS:
        worst = {}
        for _, rn, d in read_gro(os.path.join(test_data, 'nucl-dna.gro')):
            if rn not in ('DA', 'DC', 'DG', 'DT'):
                continue
            b = rn[1]
            if b not in mapping or not all(a in d for a in BASE_ATOMS[b]):
                continue
            real, W = by_key[(b, len(mapping[b]))]
            beads = np.array([np.mean([d[a] for a in mapping[b][k]], 0) for k in real])
            err = np.linalg.norm(W @ beads - np.array([d[a] for a in corners(b)]), axis=1).max()
            worst[b] = max(worst.get(b, 0.0), float(err))
        print(f'nucl-dna mapped to {model}: every ring atom back within ' + ', '.join(f'{b} {v:.3f}' for b, v in sorted(worst.items())) + ' A')
        ok &= all(v < 0.15 for v in worst.values())
    return ok


def write(path, entries):
    lines = [
        '// Generated by tools/gen_cartoon_nucleobases.py (see there). Do not edit.',
        '// The rings of the bases of the Martini nucleic acids: each ring atom a weighted sum of the real beads of its base,',
        '// N1 C2 N3 C4 C5 C6, then N7 C8 N9 for a purine.',
        '',
        'static const CartoonNucleobase cartoon_nucleobases[] = {',
    ]
    for model, b, n, real, W in entries:
        lines.append(f"    // {model} {b}: {', '.join(real)} of {n} beads")
        names = ', '.join(f'"{k}"' for k in real)
        rows = ', '.join('{' + ', '.join(f'{w:.4f}f' for w in row) + '}' for row in W)
        lines.append(f"    {{'{b}', {n}, {len(real)}, {{{names}}}, {len(W)}, {{{rows}}}}},")
    lines.append('};')
    lines.append('')
    open(path, 'w').write('\n'.join(lines))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('test_data', help="mdlib's test data (nucl-dna.gro, martini/)")
    ap.add_argument('--out', default=os.path.join(os.path.dirname(__file__), '..', 'src', 'cartoon_nucleobases.inl'))
    ap.add_argument('--check', action='store_true', help='check the weights against the frames of the test data')
    args = ap.parse_args()
    ideal = ideal_bases(os.path.join(args.test_data, 'nucl-dna.gro'))
    entries = table(ideal)
    if args.check and not check(args.test_data, ideal, entries):
        sys.exit('check failed')
    write(args.out, entries)
    print(f'wrote {len(entries)} bases to {os.path.normpath(args.out)}')


if __name__ == '__main__':
    main()

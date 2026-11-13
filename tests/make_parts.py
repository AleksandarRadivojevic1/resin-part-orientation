"""Synthetic parts with known answers, built by marching cubes on simple implicit shapes.

Used by run_tests.py. Every part is a few centimetres across; all lengths are in mm.
"""
import os, zipfile
import numpy as np
from skimage.measure import marching_cubes

DT = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])


def write_stl(path, V, F):
    T = V[F].astype(np.float32)
    n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    rec = np.zeros(len(F), DT); rec['n'] = n; rec['v'] = T
    with open(path, 'wb') as f:
        f.write(b'synthetic test part'.ljust(80) + np.uint32(len(F)).tobytes() + rec.tobytes())


def mesh(vol, h):
    V, F, _, _ = marching_cubes(np.pad(vol.astype(np.float32), 1), 0.5, spacing=(h, h, h))
    if np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum() < 0:
        F = F[:, ::-1]  # outward winding
    return V, F


def make_all(out):
    os.makedirs(out, exist_ok=True)
    h = 0.25
    X, Y, Z = np.mgrid[0:120, 0:120, 0:100] * h
    R = np.hypot(X - 15, Y - 15)

    # cup: cylinder r 10, h 20, with a pocket r 7 open at +z. Pocket area = pi * 7^2 = 153.9 mm2
    cup = mesh((R < 10) & (Z > 1) & (Z < 21) & ~((R < 7) & (Z > 6)), h); write_stl(f'{out}/cup.stl', *cup)

    # hollow sphere: r 10 with a sealed r 6 cavity. Largest cavity slice = pi * 6^2 = 113.1 mm2
    D = np.sqrt((X - 15) ** 2 + (Y - 15) ** 2 + (Z - 12) ** 2)
    write_stl(f'{out}/hollow.stl', *mesh((D < 10) & (D > 6), h))

    # slab on a pillar with three stalactites (2, 4, 8 mm) hanging toward -z: 3 islands when +z points away from the plate
    slab = (Z > 15) & (Z < 19) & (X > 3) & (X < 27) & (Y > 3) & (Y < 27)
    st = np.zeros_like(slab)
    for x, y, L in [(6, 6, 2), (24, 6, 4), (6, 24, 8)]:
        st |= (np.hypot(X - x, Y - y) < 1.2) & (Z > 15 - L) & (Z <= 15)
    stal = mesh(slab | ((R < 3) & (Z > 1) & (Z <= 15)) | st, h); write_stl(f'{out}/stalactites.stl', *stal)

    # low-poly box 30 x 20 x 10: 12 large triangles, binary and ASCII
    v = np.array([[x, y, z] for x in (0, 30) for y in (0, 20) for z in (0, 10)], float)
    f = np.array([[0, 1, 3], [0, 3, 2], [4, 6, 7], [4, 7, 5], [0, 4, 5], [0, 5, 1],
                  [2, 3, 7], [2, 7, 6], [0, 2, 6], [0, 6, 4], [1, 5, 7], [1, 7, 3]])
    write_stl(f'{out}/box.stl', v, f)
    with open(f'{out}/box_ascii.stl', 'w') as fh:
        fh.write('solid box\n')
        for t in v[f]:
            fh.write(' facet normal 0 0 0\n  outer loop\n' + ''.join(f'   vertex {p[0]} {p[1]} {p[2]}\n' for p in t)
                     + '  endloop\n endfacet\n')
        fh.write('endsolid box\n')

    # ball hanging on a 0.2 mm strap: no islands in truth, but the strap is thinner than the default 0.3 mm slice pixel
    hh = 0.05
    ax = lambda n: (np.arange(n) * hh).astype(np.float32)
    x, y, z = np.meshgrid(ax(200), ax(200), ax(300), indexing='ij', sparse=True)
    base = (z > 0.5) & (z < 2.5) & (x > 1) & (x < 9) & (y > 1) & (y < 9)
    strap = (np.abs((x - 5) + 0.07 * (z - 2.5)) < 0.1) & (np.abs(y - 5) < 1.5) & (z >= 2.5) & (z < 10)
    ball = np.sqrt((x - 4.5) ** 2 + (y - 5) ** 2 + (z - 11.5) ** 2) < 2
    write_stl(f'{out}/thin.stl', *mesh(base | strap | ball, hh))

    # a three-part "kit": the same parts stacked in assembled positions (kit/) and laid flat (flat/)
    os.makedirs(f'{out}/kit', exist_ok=True); os.makedirs(f'{out}/flat', exist_ok=True)
    z0 = 0.0
    for name, (V, F) in {'legs': stal, 'torso': cup, 'head': cup}.items():
        V = V - V.min(0); write_stl(f'{out}/flat/{name}.stl', V, F)
        write_stl(f'{out}/kit/{name}.stl', V + [0, 0, z0], F); z0 += V[:, 2].max() - 2.0  # 2 mm key overlap
    kit_height = z0 + 2.0

    # SatelLite-style 3MF holding the cup as a named mesh object, plus a component-only object
    V, F = cup
    xml = ['<?xml version="1.0" encoding="UTF-8"?>\n<model unit="millimeter"><resources>',
           '<object id="1" name="Cup_A" type="model"><mesh><vertices>']
    xml += [f'<vertex x="{p[0]:.6f}" y="{p[1]:.6f}" z="{p[2]:.6f}" />' for p in V]
    xml += ['</vertices><triangles>'] + [f'<triangle v1="{t[0]}" v2="{t[1]}" v3="{t[2]}" />' for t in F]
    xml += ['</triangles></mesh></object>',
            '<object id="2" name="Group" type="model"><components><component objectid="1"/></components></object>',
            '</resources><build><item objectid="1"/></build></model>']
    with zipfile.ZipFile(f'{out}/project.goo.3mf', 'w', zipfile.ZIP_DEFLATED) as zf:
        zf.writestr('3D/3dmodel.model', '\n'.join(xml))
    return dict(kit_height=float(kit_height))


if __name__ == '__main__':
    import sys
    print(make_all(sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), '_parts')))

#!/usr/bin/env python3
"""Tests for resin_orient.py on synthetic parts with known answers.

    pip install -r requirements.txt -r tests/requirements.txt
    python3 tests/run_tests.py

Plain asserts, no test framework. Exits non-zero on the first failure.
"""
import os, re, subprocess, sys, json, tempfile
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT); sys.path.insert(0, HERE)
import resin_orient as ro
from make_parts import make_all

TOOL = os.path.join(ROOT, 'resin_orient.py')
P = os.path.join(HERE, '_parts')
OUT = tempfile.mkdtemp(prefix='resin_orient_tests_')
Z = np.array([0, 0, 1.0])


def cli(*args, ok=True):
    r = subprocess.run([sys.executable, TOOL, *map(str, args)], capture_output=True, text=True)
    assert (r.returncode == 0) == ok, f'exit {r.returncode} for {args}\n{r.stdout}\n{r.stderr}'
    return r.stdout + r.stderr


def check(part, u, px=0.3, scale=1.0):
    V, F = ro.load(os.path.join(P, part))
    return ro.slice_check(V, F, np.asarray(u, float), scale, px, 0.05, 0.05)


def test_skill_and_tool_in_sync():
    """The tool embedded in SKILL.md is the same file as resin_orient.py."""
    md = open(os.path.join(ROOT, 'SKILL.md')).read()
    code = re.search(r'^```python\n(.*?)^```$', md, re.S | re.M).group(1)
    assert code == open(TOOL).read(), 'SKILL.md and resin_orient.py differ: re-extract the tool from SKILL.md'


def test_loaders():
    Vb, Fb = ro.load(f'{P}/box.stl'); Va, Fa = ro.load(f'{P}/box_ascii.stl')
    assert len(Fb) == len(Fa) == 12 and np.allclose(Vb, Va)
    V1, F1 = ro.load(f'{P}/cup.stl'); V2, F2 = ro.load(f'{P}/project.goo.3mf#Cup_A')
    assert len(F1) == len(F2) and np.allclose(V1.max(0) - V1.min(0), V2.max(0) - V2.min(0), atol=1e-4)
    assert 'no mesh of its own' in cli('info', f'{P}/project.goo.3mf#Group', ok=False)


def test_info():
    out = cli('info', f'{P}/hollow.stl')
    assert 'shells 2 (1 inward' in out and 'open edges 0' in out
    assert 'volume 6.0 ml' in cli('info', f'{P}/box.stl')


def test_suction_cup():
    up = check('cup.stl', Z)            # pocket faces the vat film: sealed
    assert abs(up['cup_max_area'] - np.pi * 49) < 0.05 * np.pi * 49, up['cup_max_area']
    assert up['cup_from'] is not None and 4.5 < up['cup_from'] < 5.5
    assert check('cup.stl', -Z)['cup_max_area'] < 1.0       # pocket faces the plate: open
    assert check('cup.stl', [1, 0, 0])['cup_max_area'] < 1.0  # on its side: open


def test_sealed_cavity():
    d = check('hollow.stl', Z)
    assert abs(d['cup_max_area'] - np.pi * 36) < 0.05 * np.pi * 36
    assert abs(d['cup_volume'] - 4 / 3 * np.pi * 216 / 1000) < 0.03


def test_islands_volume_and_surface():
    assert (check('stalactites.stl', Z)['islands2'], check('stalactites.stl', Z)['islands5']) == (3, 3)
    assert check('stalactites.stl', -Z)['islands2'] == 0
    V, F, fu, area, vn, E = ro.prep(f'{P}/stalactites.stl', 1.0)
    assert (ro.island_heights(V, E, vn, Z) > 0.1).sum() == 3
    assert (ro.island_heights(V, E, vn, -Z) > 0.1).sum() == 0


def test_thin_feature_needs_fine_pixels():
    """Known limit: a 0.2 mm strap is thinner than the 0.3 mm default pixel. --px 0.15 resolves it."""
    assert check('thin.stl', Z, px=0.15)['islands2'] == 0


def test_inverted_winding_is_handled():
    b = open(f'{P}/cup.stl', 'rb').read()
    dt = np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])
    R = np.frombuffer(b, dt, offset=84).copy(); R['v'] = R['v'][:, ::-1]
    open(f'{P}/cup_flipped.stl', 'wb').write(b[:84] + R.tobytes())
    assert check('cup_flipped.stl', Z)['cup_max_area'] == check('cup.stl', Z)['cup_max_area']


def test_rotation():
    rng = np.random.default_rng(0)
    for u in list(rng.normal(size=(50, 3))) + [Z, -Z, [1e-8, 0, -1], [0, 1e-7, 1]]:
        Rm = ro.rot_to_z(u); un = np.asarray(u, float) / np.linalg.norm(u)
        assert np.allclose(Rm @ Rm.T, np.eye(3), atol=1e-9) and np.isclose(np.linalg.det(Rm), 1)
        assert np.allclose(Rm @ un, Z, atol=1e-6)
    for bad in ([0, 0, 0], [0, 1]):
        try: ro.rot_to_z(bad); raise AssertionError('accepted a bad direction')
        except SystemExit: pass


def test_sweep_and_rank():
    npy = f'{OUT}/st.npy'
    cli('sweep', f'{P}/stalactites.stl', '--n', 60, '--out', npy)
    r = np.load(npy)
    assert r.shape == (61, 8) and np.allclose(r[-1, :3], Z) and tuple(r[-1, 3:5]) == (3, 3)
    best = r[np.argmin(ro.scores(r, [1, 3, 1, 0.5, 0.25]))]
    assert best[4] == 0 and best[2] < 0, 'best pose should flip the stalactites away from the plate'
    out = cli('check', f'{P}/stalactites.stl', '--from', npy, '--top', 3, '--out', f'{OUT}/c.png', '--json', f'{OUT}/c.json')
    res = json.load(open(f'{OUT}/c.json'))
    assert res[0]['label'] == 'as exported' and res[0]['islands5'] == 3 and all(d['islands5'] == 0 for d in res[1:])
    assert out.count('as exported') == 1


def test_render_fills_large_triangles():
    """A 12-triangle box must render solid, not as scattered dots."""
    V, F = ro.load(f'{P}/box.stl')
    img, _ = ro.render(V, F, np.tile([0.5, 0.5, 0.5], (len(F), 1)), 0, 89.9, px=0.1, pad=1.0)
    filled = (img < 0.999).any(2).sum() * 0.1 ** 2
    assert filled > 0.97 * 600, f'only {filled:.0f} of 600 mm2 drawn'
    assert 'islands 0' in cli('show', f'{P}/box.stl', '--u=0.3,0.2,0.93', '--out', f'{OUT}/s.png')


def test_export_round_trip():
    src, dst, s = f'{P}/stalactites.stl', f'{OUT}/exp.stl', 0.6
    cli('export', src, '--scale', s, '--u=0.18,-0.7,-0.7', '--rotz', 40, '--out', dst)
    V0, F0 = ro.load(src); V1, F1 = ro.load(dst)
    vol = lambda V, F: abs(np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum()) / 6
    assert len(F0) == len(F1) and abs(vol(V1, F1) / vol(V0, F0) - s ** 3) < 1e-4
    assert abs(V1[:, 2].min()) < 1e-4 and np.allclose(V1[:, :2].min(0), -V1[:, :2].max(0), atol=1e-3)
    assert check(os.path.relpath(dst, P), Z)['islands2'] == check('stalactites.stl', [0.18, -0.7, -0.7], scale=s)['islands2']


def test_pack():
    spec = [f'{P}/stalactites.stl@0.18,-0.7,-0.7@1', f'{P}/cup.stl@0,0,-1@1', f'{P}/box.stl@0.3,0.2,0.93@2']
    out = cli('pack', '--margin', 6, '--out', f'{OUT}/layout.json', '--png', f'{OUT}/layout.png', *spec)
    L = json.load(open(f'{OUT}/layout.json'))
    assert 'FITS' in out and len(L) == 3
    pts = []
    for p in L:  # rebuild each footprint where pack put it and measure the real gaps
        path, u, s = next(x for x in spec if x.startswith(p['part'] + '@')).split('@')
        V, F = ro.load(path); W = (V * float(s)) @ (ro.rotz(p['rotz']) @ ro.rot_to_z(list(map(float, u.split(','))))).T
        xy = W[:, :2] - (W[:, :2].min(0) + W[:, :2].max(0)) / 2 + p['center']
        T = xy[F]; t = np.linspace(0, 1, 60)[:, None, None]
        xy = np.concatenate([(T[:, i] * (1 - t) + T[:, (i + 1) % 3] * t).reshape(-1, 2) for i in range(3)])[::7]
        assert xy.min() > 6 and xy[:, 0].max() < 153.36 - 6 and xy[:, 1].max() < 77.76 - 6, 'part too close to the plate edge'
        pts.append(xy)
    from scipy.spatial import cKDTree
    for i in range(3):
        for j in range(i + 1, 3):
            assert cKDTree(pts[i]).query(pts[j])[0].min() > 6 * 1.4, 'parts closer than about 1.4 x margin'
    big = [f'{P}/box.stl@0,0,1@2'] * 5
    assert 'DOES NOT FIT' in cli('pack', '--out', f'{OUT}/none.json', *big, ok=False) and not os.path.exists(f'{OUT}/none.json')
    assert 'exceeds 165.0 mm build height' in cli('pack', '--plate', '400x400', '--out', f'{OUT}/tall.json', f'{P}/box.stl@1,0,0@5.6')


def test_size(kit_height):
    kit = [f'{P}/kit/{n}.stl' for n in ('legs', 'torso', 'head')]; flat = [f'{P}/flat/{n}.stl' for n in ('legs', 'torso', 'head')]
    out = cli('size', *kit)
    for name, mm in ro.SIZES.items():
        assert f'--scale {mm / kit_height:.4f}'[:-1] in out, (name, out)
    assert ro.SIZES == {'S': 150.0, 'M': 180.0, 'L': 200.0}
    assert 'NOT ASSEMBLED' in cli('size', *flat, ok=False)
    assert '(given)' in cli('size', *flat, '--native', kit_height)
    assert 'one part only' in cli('size', kit[0])


if __name__ == '__main__':
    print('generating synthetic parts ...', flush=True)
    meta = make_all(P)
    tests = [(k, v) for k, v in sorted(globals().items()) if k.startswith('test_')]
    for name, fn in tests:
        fn(meta['kit_height']) if name == 'test_size' else fn()
        print('ok  ', name, flush=True)
    print(f'{len(tests)} tests passed')

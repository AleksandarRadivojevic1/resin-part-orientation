#!/usr/bin/env python3
"""Resin print orientation helper: island sweep, peel/cup check, previews, plate packing, oriented STL export.

Subcommands (all lengths in mm, scale applied before analysis):
  info   PART [--scale S]                        size, volume, shells, manifold check
  size   PART [PART ...] [--up z] [--native MM]  assembled figure height and the --scale for sizes S/M/L
  sweep  PART --scale S [--n 300] --out F.npy     score build directions (surface islands, overhangs, support load)
  rank   F.npy [--top 12]                        print best directions
  check  PART --scale S (--from F.npy [--top 8] | --u x,y,z [--u ...]) --out F.png [--json F.json]
                                                 slice candidates: layer area (peel force), suction cups, volume islands
  show   PART --scale S --u x,y,z --out F.png [--views "az,el;..."]   preview: red = islands, orange = needs support
  pack   --plate 153.36x77.76 --margin 6 --out layout.json [--png layout.png] PART@x,y,z@S [...]
  export PART --scale S --u x,y,z [--rotz deg] --out F.stl
PART = binary/ASCII .stl, or a SatelLite/3MF file (first object, or FILE#name).
u = build direction: unit vector in the part's own coordinates that should point AWAY from the plate.
"""
import sys, os, re, json, zipfile, argparse
import numpy as np
import scipy.sparse as sp
from scipy import ndimage
from scipy.sparse.csgraph import connected_components
from scipy.signal import fftconvolve

try:
    import numba
    njit = numba.njit(cache=__name__ == '__main__')  # a cache written while imported as a module breaks later runs under -I
except ImportError:  # slower fallback
    njit = lambda f: f

Z_MAX = 165.0  # Mars 5 Ultra build height

# ---------------------------------------------------------------- loading
def load_stl(path):
    b = open(path, 'rb').read()
    n = int(np.frombuffer(b, '<u4', 1, 80)[0]) if len(b) >= 84 else 0
    is_ascii = b[:5].lower() == b'solid' and b'facet' in b[:4096]
    if not is_ascii and n > 0 and len(b) >= 84 + 50 * n:  # binary, tolerates trailing bytes
        rec = np.frombuffer(b, np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')]), n, 84)
        T = rec['v'].reshape(-1, 3)
    else:
        T = np.array(re.findall(rb'vertex\s+(\S+)\s+(\S+)\s+(\S+)', b), dtype=np.float32)
    if len(T) == 0:
        raise SystemExit(f'{path}: no triangles read (corrupt or unsupported STL)')
    V, inv = np.unique(T, axis=0, return_inverse=True)
    F = inv.reshape(-1, 3).astype(np.int32)
    F = F[(F[:, 0] != F[:, 1]) & (F[:, 1] != F[:, 2]) & (F[:, 0] != F[:, 2])]
    return V.astype(np.float64), F

def load_3mf(path, name=None):
    data = zipfile.ZipFile(path).read('3D/3dmodel.model')
    for m in re.finditer(rb'<object\b([^>]*)>', data):
        nm = re.search(rb'name="([^"]*)"', m.group(1))
        nm = nm.group(1).decode() if nm else ''
        if name and nm != name:
            continue
        end = data.find(b'</object>', m.end())
        seg = data[m.end():end if end >= 0 else len(data)]
        if b'<vertices>' not in seg:  # object built from <components>, not a mesh
            if name:
                raise SystemExit(f'object {name!r} in {path} has no mesh of its own; export it as STL from the slicer')
            continue
        vs = seg.index(b'<vertices>'); ve = seg.index(b'</vertices>', vs)
        ts = seg.index(b'<triangles>', ve); te = seg.index(b'</triangles>', ts)
        v = seg[vs + 10:ve]
        for a in (b'<vertex x="', b'" y="', b'" z="', b'" />', b'"/>'):
            v = v.replace(a, b' ')
        f = seg[ts + 11:te]
        for a in (b'<triangle v1="', b'" v2="', b'" v3="', b'" />', b'"/>'):
            f = f.replace(a, b' ')
        V = np.array(v.split(), dtype=np.float64).reshape(-1, 3)
        F = np.array(f.split(), dtype=np.int64).reshape(-1, 3).astype(np.int32)
        return V, F
    raise SystemExit(f'no mesh object {name!r} in {path}')

def load(part):
    path, _, name = part.partition('#')
    if path.lower().endswith('.3mf'):
        return load_3mf(path, name or None)
    return load_stl(path)

# ---------------------------------------------------------------- geometry
def rot_to_z(u):
    u = np.asarray(u, float); n = np.linalg.norm(u) if u.shape == (3,) else np.nan
    if not np.isfinite(n) or n < 1e-12:
        raise SystemExit(f'bad build direction {u.tolist()}: give x,y,z, not all zero')
    u = u / n; z = np.array([0, 0, 1.0])
    v = np.cross(u, z); c = u @ z
    if np.linalg.norm(v) < 1e-4:  # (anti)parallel to z: the general formula divides by 1 + c
        return np.eye(3) if c > 0 else np.diag([1, -1, -1.0])
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx / (1 + c)

def rotz(deg):
    a = np.radians(deg); c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])

def parse_u(s):
    u = np.array(list(map(float, s.split(',')))); return u / np.linalg.norm(u)

def face_normals(V, F):
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    if np.einsum('ij,ij->i', V[F[:, 0]], fn).sum() < 0:  # make outward
        fn = -fn
    area = np.linalg.norm(fn, axis=1) / 2
    return fn / (2 * area[:, None] + 1e-15), area

def outward(V, F):
    """Same global winding fix as face_normals, applied to F itself (needed by the slicer)."""
    fn = np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]])
    return F[:, [0, 2, 1]] if np.einsum('ij,ij->i', V[F[:, 0]], fn).sum() < 0 else F

def vertex_normals(V, F, fu, area):
    n = len(V); vn = np.zeros((n, 3)); fn = fu * area[:, None]
    for k in range(3):
        for c in range(3):
            vn[:, c] += np.bincount(F[:, k], weights=fn[:, c], minlength=n)
    return vn / (np.linalg.norm(vn, axis=1, keepdims=True) + 1e-15)

def edges(V, F):
    n = len(V)
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]).astype(np.int64); e.sort(1)
    e = np.unique(e[:, 0] * n + e[:, 1])
    return np.stack([e // n, e % n], 1)

def shells(V, F):
    """Connected surface shells: (count, per-shell signed volume in mm3, vertex labels)."""
    E = edges(V, F)
    G = sp.coo_matrix((np.ones(len(E)), (E[:, 0], E[:, 1])), shape=(len(V), len(V)))
    k, lab = connected_components(G, directed=False)
    sv = np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])) / 6
    vol = np.bincount(lab[F[:, 0]], weights=sv, minlength=k)
    used = np.bincount(lab[F[:, 0]], minlength=k) > 0
    return int(used.sum()), vol[used], lab

# ---------------------------------------------------------------- islands (surface graph, fast; one shell only)
@njit
def _persist(E, z, order, anc, n):
    parent = np.arange(n); bv = np.arange(n); a = anc.copy(); per = np.full(n, -1.0)
    for i in order:
        u = E[i, 0]; v = E[i, 1]; t = max(z[u], z[v])
        ru = u
        while parent[ru] != ru: ru = parent[ru]
        rv = v
        while parent[rv] != rv: rv = parent[rv]
        x = u
        while parent[x] != ru:
            nx = parent[x]; parent[x] = ru; x = nx
        x = v
        while parent[x] != rv:
            nx = parent[x]; parent[x] = rv; x = nx
        if ru == rv: continue
        if a[ru] and not a[rv]: old, young = ru, rv
        elif a[rv] and not a[ru]: old, young = rv, ru
        elif z[bv[ru]] <= z[bv[rv]]: old, young = ru, rv
        else: old, young = rv, ru
        if not a[young]: per[bv[young]] = t - z[bv[young]]
        parent[young] = old; a[old] = a[old] or a[young]
    return per

def island_heights(V, E, vn, u):
    """Per-vertex height (mm) an unsupported island floats before joining the part; -1 = not an island start.
    Up-facing local minima (cavity floors) and the global minimum count as anchored.
    Walks surface edges only, so on meshes with several overlapping shells it misses and invents islands:
    use `check` (volume slicing) for those."""
    n = len(V); z = V @ u
    low = np.full(n, np.inf)
    np.minimum.at(low, E[:, 0], z[E[:, 1]]); np.minimum.at(low, E[:, 1], z[E[:, 0]])
    lm = z < low; anc = (z <= low) & (vn @ u >= 0); anc[np.argmin(z)] = True  # <=: flat floors (key sockets) are plateaus
    order = np.argsort(np.maximum(z[E[:, 0]], z[E[:, 1]]))
    return _persist(E, z, order, anc, n)

def fib_dirs(n):
    i = np.arange(n) + 0.5; phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.vstack([np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1), [[0, 0, 1.0]]])

# ---------------------------------------------------------------- slicing (volume; handles overlapping shells)
def _crossings(W, F, px, dz, lo):
    """Signed Z-ray crossings at pixel centres: (I, J, K, delta). delta +1 entering solid, -1 leaving."""
    jit = np.array([0.000371, 0.000613]) * px  # keeps pixel centres off shared edges
    out = [[], [], [], []]
    for s in range(0, len(F), 400_000):
        T = W[F[s:s + 400_000]]
        x0 = (T[:, :, 0].min(1) - lo[0] - jit[0]) / px - 0.5; x1 = (T[:, :, 0].max(1) - lo[0] - jit[0]) / px - 0.5
        y0 = (T[:, :, 1].min(1) - lo[1] - jit[1]) / px - 0.5; y1 = (T[:, :, 1].max(1) - lo[1] - jit[1]) / px - 0.5
        i0 = np.ceil(x0).astype(np.int64); j0 = np.ceil(y0).astype(np.int64)
        ni = np.maximum(np.floor(x1).astype(np.int64) - i0 + 1, 0); nj = np.maximum(np.floor(y1).astype(np.int64) - j0 + 1, 0)
        cnt = ni * nj; tot = int(cnt.sum())
        if tot == 0: continue
        t = np.repeat(np.arange(len(T)), cnt); k = np.arange(tot) - np.repeat(np.cumsum(cnt) - cnt, cnt)
        I = i0[t] + k % ni[t]; J = j0[t] + k // ni[t]
        cx = lo[0] + (I + 0.5) * px + jit[0]; cy = lo[1] + (J + 0.5) * px + jit[1]
        A, B, C = T[t, 0], T[t, 1], T[t, 2]
        den = (B[:, 0] - A[:, 0]) * (C[:, 1] - A[:, 1]) - (B[:, 1] - A[:, 1]) * (C[:, 0] - A[:, 0])
        ok = np.abs(den) > 1e-12; den = np.where(ok, den, 1.0)
        w1 = ((cx - A[:, 0]) * (C[:, 1] - A[:, 1]) - (cy - A[:, 1]) * (C[:, 0] - A[:, 0])) / den
        w2 = ((B[:, 0] - A[:, 0]) * (cy - A[:, 1]) - (B[:, 1] - A[:, 1]) * (cx - A[:, 0])) / den
        w0 = 1 - w1 - w2; hit = ok & (w0 >= 0) & (w1 >= 0) & (w2 >= 0)
        z = w0 * A[:, 2] + w1 * B[:, 2] + w2 * C[:, 2]
        out[0].append(I[hit]); out[1].append(J[hit])
        out[2].append(np.ceil((z[hit] - lo[2]) / dz - 0.5).astype(np.int64))
        out[3].append(np.where(den[hit] > 0, -1, 1).astype(np.int16))  # up-facing = leaving
    return [np.concatenate(o) if o else np.zeros(0, np.int64) for o in out]

class _UF:
    def __init__(self): self.p = []; self.birth = []; self.flag = []
    def add(self, n, birth, flag):
        s = len(self.p); self.p += range(s, s + n); self.birth += [birth] * n; self.flag += [flag] * n
        return np.arange(s, s + n)
    def find(self, x):
        p = self.p
        while p[x] != x: p[x] = p[p[x]]; x = p[x]
        return x

def _pairs(a, b, nb):
    """Distinct (a, b) label pairs; 1-D unique on a packed key is much faster than unique(axis=0)."""
    key = np.unique(a.astype(np.int64) * (nb + 1) + b)
    return zip((key // (nb + 1)).tolist(), (key % (nb + 1)).tolist())

def slice_check(V, F, u, scale, px, dz, layer):
    """Slice the part (union of all shells) along u. Returns per-layer solid area, sealed (cup) area, and
    float heights of volume islands."""
    W = (V * scale) @ rot_to_z(u).T; W -= W.min(0); F = outward(W, F)
    lo = np.array([-px, -px, 0.0]); hi = W.max(0)
    nx = int(np.ceil((hi[0] + px) / px)) + 1; ny = int(np.ceil((hi[1] + px) / px)) + 1; nz = int(np.ceil(hi[2] / dz)) + 1
    I, J, K, D = _crossings(W, F, px, dz, lo)
    K = np.clip(K, 0, nz); o = np.argsort(K, kind='stable'); I, J, K, D = I[o], J[o], K[o], D[o]
    cut = np.searchsorted(K, np.arange(nz + 1))
    wind = np.zeros((nx, ny), np.int16)
    s8 = np.ones((3, 3), bool); px2 = px * px
    area = np.zeros(nz); sealed = np.zeros(nz)
    S = _UF(); E = _UF(); floats = []
    plab = pids = pelab = peids = None; started = False
    for k in range(nz):
        np.add.at(wind, (I[cut[k]:cut[k + 1]], J[cut[k]:cut[k + 1]]), D[cut[k]:cut[k + 1]])
        solid = wind > 0; area[k] = solid.sum() * px2
        if not started and area[k] == 0: continue
        first = not started; started = True
        # solid components: track islands
        lab, n = ndimage.label(solid, s8); ids = S.add(n, k, first)
        if plab is not None and n:
            m = (lab > 0) & (plab > 0)
            for a, b in _pairs(plab[m], lab[m], n):
                ra, rb = S.find(pids[a - 1]), S.find(ids[b - 1])
                if ra == rb: continue
                if S.flag[ra] != S.flag[rb]: old, young = (ra, rb) if S.flag[ra] else (rb, ra)
                else: old, young = (ra, rb) if S.birth[ra] <= S.birth[rb] else (rb, ra)
                if not S.flag[young] and k > S.birth[young]: floats.append((k - S.birth[young]) * dz)
                S.p[young] = old; S.flag[old] = S.flag[old] or S.flag[young]
        plab, pids = lab, ids
        # empty components: a pocket is sealed if no path through printed layers reaches open resin
        elab, en = ndimage.label(~solid)
        border = np.unique(np.concatenate([elab[0], elab[-1], elab[:, 0], elab[:, -1]]))
        eids = E.add(en, k, first)  # everything in the first layer opens to the space under the part
        for b in border[border > 0]: E.flag[eids[b - 1]] = True
        if pelab is not None and en:
            m = (elab > 0) & (pelab > 0)
            for a, b in _pairs(pelab[m], elab[m], en):
                ra, rb = E.find(peids[a - 1]), E.find(eids[b - 1])
                if ra != rb: E.p[rb] = ra; E.flag[ra] = E.flag[ra] or E.flag[rb]
        if en:
            closed = np.array([not E.flag[E.find(x)] for x in eids])
            if closed.any():
                sz = ndimage.sum_labels(np.ones_like(elab), elab, np.arange(1, en + 1))
                sealed[k] = sz[closed].sum() * px2
        pelab, peids = elab, eids
    for r in range(len(S.p)):  # bodies that never join the anchored part float to the top
        if S.p[r] == r and not S.flag[r]: floats.append((nz - S.birth[r]) * dz)
    floats = np.array(floats)
    h = hi[2]; zz = (np.arange(nz) + 0.5) * dz; top = zz >= 0.8 * h
    kmax = int(np.argmax(area))
    cup = sealed > 1.0
    return dict(u=[round(float(x), 3) for x in u], height=round(float(h), 1),
                islands2=int((floats > 2 * layer + 1e-9).sum()), islands5=int((floats > 5 * layer + 1e-9).sum()),
                max_area=round(area.max(), 1), max_area_at=round(zz[kmax] / h, 2),
                late_max_area=round(area[top].max() if top.any() else 0, 1),
                cup_max_area=round(sealed.max(), 1), cup_volume=round(sealed.sum() * dz / 1000, 3),
                cup_from=round(zz[cup].min(), 1) if cup.any() else None, cup_to=round(zz[cup].max(), 1) if cup.any() else None,
                _area=area, _sealed=sealed, _z=zz)

# ---------------------------------------------------------------- render (numpy z-buffer splats)
def _view(az, el):
    a, e = np.radians(az), np.radians(el)
    d = -np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    up = np.array([0, 0, 1.0]) if abs(d[2]) < 0.99 else np.array([0, 1.0, 0])
    r = np.cross(d, up); r /= np.linalg.norm(r); return np.stack([r, np.cross(r, d), -d])

def render(V, F, colors, az, el, px=0.08, pad=2.0, seed=0):
    R = _view(az, el); rng = np.random.default_rng(seed); P = V @ R.T
    ox, oy = P[:, 0].min() - pad, P[:, 1].min() - pad
    W = int((P[:, 0].max() + pad - ox) / px) + 1; H = int((P[:, 1].max() + pad - oy) / px) + 1
    zbuf = np.full(W * H, -np.inf); img = np.ones((W * H, 3))
    for s in range(0, len(F), 1_000_000):
        f = F[s:s + 1_000_000]; T = P[f]; c = colors[s:s + 1_000_000]
        nrm = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0]); nn = nrm / (np.linalg.norm(nrm, axis=1, keepdims=True) + 1e-12)
        shade = np.abs(nn[:, 2]) * 0.75 + 0.25 * np.clip(nn[:, 1] * 0.5 + 0.5, 0, 1)
        k = np.ceil(np.abs(nrm[:, 2]) / px ** 2 * 1.5).astype(np.int64) + 1  # about 3 samples per covered pixel
        k = np.minimum(np.where(k > 400, 2 * k, k), 4_000_000)  # no low cap: big flat faces (low-poly, decimated) must fill in
        cs = np.cumsum(k); cut = np.unique(np.r_[0, np.searchsorted(cs, np.arange(4_000_000, cs[-1], 4_000_000)), len(f)])
        for i0, i1 in zip(cut[:-1], cut[1:]):  # about 4M samples per pass, however large the triangles are
            idx = np.repeat(np.arange(i0, i1), k[i0:i1])
            r1 = np.sqrt(rng.random(len(idx))); r2 = rng.random(len(idx))
            b = np.stack([1 - r1, r1 * (1 - r2), r1 * r2], 1)
            pts = np.einsum('ij,ijk->ik', b, T[idx])
            ix = ((pts[:, 0] - ox) / px).astype(int); iy = ((pts[:, 1] - oy) / px).astype(int)
            pix = iy * W + ix; o = np.lexsort((-pts[:, 2], pix)); pix, zz, ii = pix[o], pts[o, 2], idx[o]
            keep = np.r_[True, pix[1:] != pix[:-1]]; pix, zz, ii = pix[keep], zz[keep], ii[keep]
            better = zz > zbuf[pix]
            zbuf[pix[better]] = zz[better]; img[pix[better]] = shade[ii[better], None] * c[ii[better]]
    return np.clip(img.reshape(H, W, 3)[::-1], 0, 1), (R, ox, oy, px, H, zbuf.reshape(H, W)[::-1])

def project(points, meta):
    R, ox, oy, px, H = meta[:5]; P = points @ R.T
    return (P[:, 0] - ox) / px, H - 1 - (P[:, 1] - oy) / px, P[:, 2]

def visible(x, y, d, meta, tol=0.5):
    zb = meta[5]; r = np.clip(np.round(y).astype(int), 0, zb.shape[0] - 1); c = np.clip(np.round(x).astype(int), 0, zb.shape[1] - 1)
    return d >= zb[r, c] - tol

# ---------------------------------------------------------------- commands
def prep(part, scale):
    V, F = load(part); V = V * scale
    fu, area = face_normals(V, F)
    return V, F, fu, area, vertex_normals(V, F, fu, area), edges(V, F)

def cmd_info(a):
    V, F = load(a.part); V = V * a.scale
    vol = abs(np.einsum('ij,ij->i', V[F[:, 0]], np.cross(V[F[:, 1]], V[F[:, 2]])).sum()) / 6
    e = np.concatenate([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]]).astype(np.int64); e.sort(1)
    _, c = np.unique(e[:, 0] * len(V) + e[:, 1], return_counts=True)
    ns, sv, _ = shells(V, F); sgn = 1 if sv.sum() >= 0 else -1
    size = V.max(0) - V.min(0)
    print(f'{a.part}: {len(F)} tris, size {np.round(size, 1)} mm, volume {vol / 1000:.1f} ml, '
          f'open edges {(c == 1).sum()}, non-manifold edges {(c > 2).sum()}, '
          f'shells {ns} ({(sgn * sv < 0).sum()} inward: cavities or flipped)')
    if ns > 1:
        print('  several shells: surface island counts from sweep/show are unreliable here, trust `check`')

SIZES = {'S': 150.0, 'M': 180.0, 'L': 200.0}  # height of the assembled figure, mm

def cmd_size(a):
    ax = 'xyz'.index(a.up); lo = np.full(3, np.inf); hi = -lo; z0 = []; tall = 0.0
    for part in a.parts:
        V, _ = load(part); l, h = V.min(0), V.max(0)
        lo = np.minimum(lo, l); hi = np.maximum(hi, h); z0.append(l[ax]); tall = max(tall, (h - l)[ax])
        print(f'{part}: {a.up} from {l[ax]:.1f} to {h[ax]:.1f} mm, size {np.round(h - l, 1)}')
    ext = hi - lo; native = a.native if a.native else ext[ax]
    print(f'all {len(a.parts)} parts together: {np.round(ext, 1)} mm (x y z), up axis {a.up}')
    if not a.native and len(a.parts) > 1:
        if np.ptp(z0) < 0.5:
            raise SystemExit('NOT ASSEMBLED: every part starts at the same height, so the parts are laid out for printing '
                             'and the figure height cannot be measured. Pass the real one with --native MM.')
        if ext[ax] < 1.02 * tall:
            print('WARNING: all parts together are no taller than the tallest part. If the parts are not in their '
                  'assembled positions this is not the figure height: pass the real one with --native MM.')
    if not a.native and len(a.parts) == 1:
        print('NOTE: one part only. This is the figure height only if the kit is a single piece.')
    print(f'figure height at 100%: {native:.1f} mm' + (' (given)' if a.native else ' (measured)'))
    for k, mm in SIZES.items():
        print(f'  {k} = {mm / 10:g} cm: --scale {mm / native:.4f}  ({100 * mm / native:.1f}%)')
    print('  other heights: scale = wanted height in mm / figure height at 100%')

def cmd_sweep(a):
    V, F, fu, area, vn, E = prep(a.part, a.scale)
    ns = shells(V, F)[0]
    if ns > 1: print(f'warning: {ns} shells, surface island counts are approximate; confirm finalists with `check`')
    C = V[F].mean(1); U = fib_dirs(a.n); res = []
    for k, u in enumerate(U):
        per = island_heights(V, E, vn, u)
        z = V @ u; down = fu @ u < -0.7071; zc = C @ u - z.min()
        res.append([*u, (per > 2 * a.layer).sum(), (per > 5 * a.layer).sum(), area[down].sum(),
                    (area[down] * (zc[down] + a.lift)).sum() / 1000, z.max() - z.min()])
        if k % 50 == 0: print(f'{k}/{len(U)}', flush=True)
    np.save(a.out, np.array(res)); print('saved', a.out, '(last row = as exported, u=0,0,1)')

def scores(r, weights):
    cols = r[:, 3:8]; med = np.maximum(np.median(cols, 0), [1, 1, 1, 0.1, 1])  # floor: medians are often 0 islands
    return (cols / med) @ np.asarray(weights, float)

WEIGHTS = '1,3,1,0.5,0.25'  # islands>2, islands>5, down-facing area, support load, height; each / its median

def cmd_rank(a):
    r = np.load(a.file); sc = scores(r, list(map(float, a.weights.split(','))))
    print('columns: u | islands >2 layers | >5 layers | down-facing area mm2 | support load | height mm | score (lower = better)')
    print('as exported:', np.round(r[-1, :3], 2), r[-1, 3:5].astype(int), r[-1, 5:].round(1), round(sc[-1], 2))
    for i in np.argsort(sc)[:a.top]:
        print(np.round(r[i, :3], 2), r[i, 3:5].astype(int), r[i, 5:].round(1), round(sc[i], 2))

def cmd_check(a):
    V, F = load(a.part)
    cands = []
    if a.sweep:
        r = np.load(a.sweep); sc = scores(r, list(map(float, a.weights.split(','))))
        cands += [('as exported', r[-1, :3])] + [(f'rank {j + 1}', r[i, :3]) for j, i in enumerate(np.argsort(sc)[:a.top])
                                                 if i != len(r) - 1]  # last row is the as-exported pose: slice it once
    cands += [('given', parse_u(s)) for s in (a.u or [])]
    if not cands: raise SystemExit('give --from F.npy or --u')
    res = []
    print('label | u | height | vol islands >2 / >5 | max layer area mm2 @ height frac | max area in top 20% | '
          'cup max area mm2, volume ml, z range')
    for lab, u in cands:
        d = slice_check(V, F, np.asarray(u, float), a.scale, a.px, a.dz, a.layer); d['label'] = lab; res.append(d)
        cz = f"{d['cup_from']}-{d['cup_to']} mm" if d['cup_from'] is not None else '-'
        print(f"{lab:12s} {np.round(u, 2)} {d['height']:6.1f} | {d['islands2']:4d} / {d['islands5']:4d} | "
              f"{d['max_area']:8.1f} @ {d['max_area_at']:.2f} | {d['late_max_area']:8.1f} | "
              f"{d['cup_max_area']:7.1f}, {d['cup_volume']:.3f}, {cz}", flush=True)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    cols = min(4, len(res)); rows = (len(res) + cols - 1) // cols
    fig, axs = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.2 * rows), squeeze=False)
    for ax, d in zip(axs.flat, res):
        ax.fill_between(d['_z'], d['_area'], color='#8a93b8', alpha=0.6, lw=0, label='layer area')
        if d['_sealed'].max() > 0: ax.fill_between(d['_z'], d['_sealed'], color='#d0312d', alpha=0.8, lw=0, label='sealed cup')
        ax.axvspan(0.8 * d['height'], d['height'], color='#f0a030', alpha=0.15, lw=0)
        ax.set_title(f"{d['label']}  u={d['u']}\nislands>5: {d['islands5']}  cup: {d['cup_max_area']} mm2", fontsize=9)
        ax.set_xlabel('height above plate, mm', fontsize=8); ax.set_ylabel('mm2', fontsize=8); ax.tick_params(labelsize=7)
    for ax in list(axs.flat)[len(res):]: ax.axis('off')
    axs.flat[0].legend(fontsize=7, loc='upper left')
    fig.suptitle(f'{os.path.basename(a.part)} at {a.scale:g}: layer area (peel force) and sealed cups; shaded = last 20% of height', fontsize=10)
    fig.tight_layout(); fig.savefig(a.out, dpi=80); print('saved', a.out)
    if a.json:
        json.dump([{k: v for k, v in d.items() if not k.startswith('_')} for d in res], open(a.json, 'w'), indent=1)

def cmd_show(a):
    V, F, fu, area, vn, E = prep(a.part, a.scale)
    u = parse_u(a.u)
    per = island_heights(V, E, vn, u); isl = np.nonzero(per > 2 * a.layer)[0]
    W = V @ rot_to_z(u).T; W -= W.min(0)
    fz, _ = face_normals(W, F)
    col = np.tile([0.70, 0.72, 0.82], (len(F), 1)); col[fz[:, 2] < -0.7071] = (0.98, 0.62, 0.25)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    views = [tuple(map(float, v.split(','))) for v in a.views.split(';')]
    fig, axs = plt.subplots(1, len(views), figsize=(6 * len(views), 6)); axs = np.atleast_1d(axs)
    for ax, (az, el) in zip(axs, views):
        img, meta = render(W, F, col, az, el, px=a.px); ax.imshow(img)
        x, y, d = project(W[isl], meta); front = visible(x, y, d, meta)
        ax.scatter(x[front], y[front], s=30, facecolors='none', edgecolors='red', linewidths=1.4)
        ax.set_title(f'az {az:.0f} el {el:.0f}'); ax.axis('off')
    fig.suptitle(f'{os.path.basename(a.part)}  u={np.round(u, 2)}  islands {len(isl)} (red)  '
                 f'orange = needs support  height {W[:, 2].max():.0f} mm  (plate at bottom)')
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(a.out, dpi=60); print('saved', a.out, '| islands', len(isl))

def footprint(V, F, u, scale, ang, margin, px):
    W = (V * scale) @ (rotz(ang) @ rot_to_z(u)).T; lo = W.min(0)
    ij = np.floor((W[:, :2] - lo[:2]) / px).astype(int)
    m = np.zeros(ij.max(0) + 1, bool); m[ij[:, 0], ij[:, 1]] = True
    I, J, _, _ = _crossings(W, F, px, 1.0, lo)  # fill triangle interiors too: low-poly parts have few vertices
    ok = (I >= 0) & (J >= 0) & (I < m.shape[0]) & (J < m.shape[1]); m[I[ok], J[ok]] = True
    m = ndimage.binary_fill_holes(ndimage.binary_closing(m, iterations=3))
    r = int(round(margin / px)); m = np.pad(m, r)
    return ndimage.binary_dilation(m, ndimage.generate_binary_structure(2, 1), iterations=r), W[:, 2].max() - W[:, 2].min()

def cmd_pack(a):
    px = 0.5; PW, PH = map(float, a.plate.lower().split('x')); GW, GH = int(PW / px), int(PH / px); e = 2
    parts = []
    for spec in a.parts:
        path, u, s = spec.split('@'); V, F = load(path); u = np.array(list(map(float, u.split(','))))
        masks = [(ang, *footprint(V, F, u, float(s), ang, a.margin, px)) for ang in range(0, 360, 20)]
        parts.append((path, masks))
    parts.sort(key=lambda p: -max(m.sum() for _, m, _ in p[1]))
    def place(grid, m):
        mw, mh = m.shape
        if mw > GW - 2 * e or mh > GH - 2 * e: return None
        occ = fftconvolve(grid.astype(float), m[::-1, ::-1].astype(float), mode='valid') > 0.5
        ok = np.argwhere(~occ[e:GW - e - mw + 1, e:GH - e - mh + 1])
        return None if len(ok) == 0 else tuple(ok[np.lexsort((ok[:, 1], ok[:, 0]))[0]] + e)
    states = [(0, np.zeros((GW, GH), bool), [])]  # beam search: keeps the `beam` tightest partial layouts
    for path, masks in parts:
        nxt = []
        for _, grid, placed in states:
            for ang, m, h in masks:
                pos = place(grid, m)
                if pos is None: continue
                g = grid.copy(); g[pos[0]:pos[0] + m.shape[0], pos[1]:pos[1] + m.shape[1]] |= m
                cx, cy = (pos[0] + m.shape[0] / 2) * px, (pos[1] + m.shape[1] / 2) * px
                nxt.append((int(np.nonzero(g.any(1))[0].max()) + 1, g,
                            placed + [dict(part=path, rotz=ang, x=float(pos[0] * px), y=float(pos[1] * px),
                                           center=[round(float(cx), 1), round(float(cy), 1)], height=round(float(h), 1))]))
        if not nxt: raise SystemExit(f'DOES NOT FIT: no room for {path} on one plate with margin {a.margin:g} mm; no layout written')
        nxt.sort(key=lambda s: s[0]); states = nxt[:a.beam]
    span, grid, placed = states[0]
    print(f'FITS: uses {span * px:.0f} of {PW} mm'); [print(p) for p in placed]
    tall = max(p['height'] for p in placed)
    if tall + a.lift > Z_MAX: print(f'WARNING: tallest part {tall} mm + {a.lift} mm lift exceeds {Z_MAX} mm build height')
    json.dump(placed, open(a.out, 'w'), indent=1)
    if a.png:
        import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(PW / 15, PH / 15))
        ax.imshow(grid.T, origin='lower', extent=(0, PW, 0, PH), cmap='Blues', vmin=0, vmax=1.6)
        for p in placed:
            ax.text(*p['center'], f"{os.path.basename(p['part'])}\nrotz {p['rotz']}  h {p['height']}",
                    ha='center', va='center', fontsize=7)
        ax.set_xlim(0, PW); ax.set_ylim(0, PH); ax.set_title(f'plate {PW} x {PH} mm, margin {a.margin} mm (shaded = footprint + margin)', fontsize=8)
        fig.tight_layout(); fig.savefig(a.png, dpi=90); print('saved', a.png)

def cmd_export(a):
    V, F = load(a.part); u = parse_u(a.u)
    W = (V * a.scale) @ (rotz(a.rotz) @ rot_to_z(u)).T
    W[:, :2] -= (W[:, :2].min(0) + W[:, :2].max(0)) / 2; W[:, 2] -= W[:, 2].min()
    T = W[F].astype(np.float32); n = np.cross(T[:, 1] - T[:, 0], T[:, 2] - T[:, 0])
    n /= np.linalg.norm(n, axis=1, keepdims=True) + 1e-12
    rec = np.zeros(len(F), np.dtype([('n', '<f4', 3), ('v', '<f4', (3, 3)), ('a', '<u2')])); rec['n'] = n; rec['v'] = T
    with open(a.out, 'wb') as f:
        f.write(f'scale {a.scale} pre-oriented'.encode()[:80].ljust(80)); f.write(np.uint32(len(F)).tobytes()); f.write(rec.tobytes())
    print('saved', a.out, 'size', np.round(W.max(0) - W.min(0), 1), 'mm')

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    def common(q, need_u=False):
        q.add_argument('part'); q.add_argument('--scale', type=float, default=1.0)
        q.add_argument('--layer', type=float, default=0.05); q.add_argument('--lift', type=float, default=5.0)
        if need_u: q.add_argument('--u', required=True)
    q = sub.add_parser('info'); common(q)
    q = sub.add_parser('size'); q.add_argument('parts', nargs='+'); q.add_argument('--up', choices='xyz', default='z')
    q.add_argument('--native', type=float, help='figure height at 100%% in mm, when the parts are not assembled')
    q = sub.add_parser('sweep'); common(q); q.add_argument('--n', type=int, default=300); q.add_argument('--out', required=True)
    q = sub.add_parser('rank'); q.add_argument('file'); q.add_argument('--top', type=int, default=12); q.add_argument('--weights', default=WEIGHTS)
    q = sub.add_parser('check'); common(q); q.add_argument('--from', dest='sweep'); q.add_argument('--top', type=int, default=8)
    q.add_argument('--weights', default=WEIGHTS); q.add_argument('--u', action='append')
    q.add_argument('--px', type=float, default=0.3); q.add_argument('--dz', type=float, default=None)
    q.add_argument('--out', required=True); q.add_argument('--json')
    q = sub.add_parser('show'); common(q, True); q.add_argument('--out', required=True)
    q.add_argument('--views', default='-90,10;0,10;90,10;0,-45'); q.add_argument('--px', type=float, default=0.08)
    q = sub.add_parser('pack'); q.add_argument('--plate', default='153.36x77.76'); q.add_argument('--margin', type=float, default=6.0)
    q.add_argument('--beam', type=int, default=24); q.add_argument('--lift', type=float, default=5.0)
    q.add_argument('--out', required=True); q.add_argument('--png'); q.add_argument('parts', nargs='+')
    q = sub.add_parser('export'); common(q, True); q.add_argument('--rotz', type=float, default=0.0); q.add_argument('--out', required=True)
    a = p.parse_args()
    if a.cmd == 'check' and a.dz is None: a.dz = a.layer
    {'info': cmd_info, 'size': cmd_size, 'sweep': cmd_sweep, 'rank': cmd_rank, 'check': cmd_check, 'show': cmd_show,
     'pack': cmd_pack, 'export': cmd_export}[a.cmd](a)

if __name__ == '__main__':
    main()

import json, math, sys, zlib, random

src, extra_src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
random.seed(14)
els = json.load(open(src, encoding='utf8'))['elements']
try:
    els += json.load(open(extra_src, encoding='utf8'))['elements']
except Exception as e:
    print('no extra', e)
nodes = {}
for e in els:
    if e['type'] == 'node' and 'lat' in e: nodes[e['id']] = e
ways = {e['id']: e for e in els if e['type'] == 'way'}
rels = {e['id']: e for e in els if e['type'] == 'relation'}

LAT0, LON0 = 49.1840, -0.3630
KX = math.cos(math.radians(LAT0)) * 111320
KZ = 110540
HX, HZ = 1200, 900  # demi-dimensions de la maquette (m)


def xy(nid):
    n = nodes[nid]
    return (round((n['lon'] - LON0) * KX, 1), round(-(n['lat'] - LAT0) * KZ, 1))


def pts_of(w):
    return [xy(n) for n in w['nodes'] if n in nodes]


def ring_of(w):
    p = pts_of(w)
    if len(p) > 1 and p[0] == p[-1]: p = p[:-1]
    return p


def inside_box(p, m=0):
    return abs(p[0]) <= HX - m and abs(p[1]) <= HZ - m


def clip_poly(poly):
    def edge(pts, inside, inter):
        out = []
        for i in range(len(pts)):
            a, b = pts[i - 1], pts[i]
            ia, ib = inside(a), inside(b)
            if ib:
                if not ia: out.append(inter(a, b))
                out.append(b)
            elif ia: out.append(inter(a, b))
        return out
    ix = lambda v: (lambda a, b: (v, round(a[1] + (b[1] - a[1]) * (v - a[0]) / (b[0] - a[0]), 1)))
    iz = lambda v: (lambda a, b: (round(a[0] + (b[0] - a[0]) * (v - a[1]) / (b[1] - a[1]), 1), v))
    p = poly
    for ins, inter in [(lambda q: q[0] >= -HX, ix(-HX)), (lambda q: q[0] <= HX, ix(HX)),
                       (lambda q: q[1] >= -HZ, iz(-HZ)), (lambda q: q[1] <= HZ, iz(HZ))]:
        if not p: break
        p = edge(p, ins, inter)
    return p


def clip_line(pts):
    """Liang-Barsky par segment ; renvoie des morceaux de polyligne dans la boîte."""
    pieces, cur = [], []
    for i in range(len(pts) - 1):
        (x0, z0), (x1, z1) = pts[i], pts[i + 1]
        dx, dz = x1 - x0, z1 - z0
        t0, t1, ok = 0.0, 1.0, True
        for p, q in ((-dx, x0 + HX), (dx, HX - x0), (-dz, z0 + HZ), (dz, HZ - z0)):
            if p == 0:
                if q < 0: ok = False; break
            else:
                t = q / p
                if p < 0: t0 = max(t0, t)
                else: t1 = min(t1, t)
        if not ok or t0 > t1:
            if len(cur) > 1: pieces.append(cur)
            cur = []; continue
        a = (round(x0 + dx * t0, 1), round(z0 + dz * t0, 1))
        b = (round(x0 + dx * t1, 1), round(z0 + dz * t1, 1))
        if cur and cur[-1] == a: cur.append(b)
        else:
            if len(cur) > 1: pieces.append(cur)
            cur = [a, b]
        if t1 < 1:
            pieces.append(cur); cur = []
    if len(cur) > 1: pieces.append(cur)
    return pieces


def join_rings(way_ids):
    segs = [pts_of(ways[i]) for i in way_ids if i in ways]
    segs = [s for s in segs if len(s) > 1]
    rings = []
    while segs:
        cur = segs.pop(0)
        changed = True
        while cur[0] != cur[-1] and changed:
            changed = False
            for k, s in enumerate(segs):
                if s[0] == cur[-1]: cur += s[1:]
                elif s[-1] == cur[-1]: cur += s[::-1][1:]
                elif s[-1] == cur[0]: cur = s + cur[1:]
                elif s[0] == cur[0]: cur = s[::-1] + cur[1:]
                else: continue
                segs.pop(k); changed = True; break
        if cur[0] == cur[-1]: cur = cur[:-1]
        if len(cur) >= 3: rings.append(cur)
    return rings


def centroid(r): return (sum(p[0] for p in r) / len(r), sum(p[1] for p in r) / len(r))
def area(r): return abs(sum(r[i - 1][0] * r[i][1] - r[i][0] * r[i - 1][1] for i in range(len(r)))) / 2
def flat(r):
    o = []
    for p in r: o += [p[0], p[1]]
    return o


def pip(pt, ring):
    x, z = pt; c = False
    for i in range(len(ring)):
        a, b = ring[i - 1], ring[i]
        if (a[1] > z) != (b[1] > z) and x < (b[0] - a[0]) * (z - a[1]) / (b[1] - a[1]) + a[0]: c = not c
    return c


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3: return pts
    def cr(o, a, b): return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cr(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cr(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def obb(r):
    h = hull(r)
    best = None
    for i in range(len(h)):
        a, b = h[i - 1], h[i]
        ang = math.atan2(b[1] - a[1], b[0] - a[0])
        c, s = math.cos(ang), math.sin(ang)
        us = [p[0] * c + p[1] * s for p in h]
        vs = [-p[0] * s + p[1] * c for p in h]
        A = (max(us) - min(us)) * (max(vs) - min(vs))
        if best is None or A < best[0]:
            cu, cv = (max(us) + min(us)) / 2, (max(vs) + min(vs)) / 2
            best = (A, cu * c - cv * s, cu * s + cv * c, ang, max(us) - min(us), max(vs) - min(vs))
    A, cx, cz, ang, L, W = best
    if W > L: ang += math.pi / 2; L, W = W, L
    return A, [round(cx, 1), round(cz, 1), round(ang, 4), round(L, 1), round(W, 1)]


# --- château : polygone du tertre ---
castle = []
for w in ways.values():
    if w.get('tags', {}).get('historic') == 'castle':
        r = ring_of(w)
        if len(r) >= 3 and inside_box(centroid(r)): castle.append(r)
castle.sort(key=area, reverse=True)
castle = castle[:1]
MOUND = 7.0
def on_mound(p): return bool(castle) and pip(p, castle[0])


def height_of(t, oid):
    h = t.get('height', '').replace('m', '').strip()
    try:
        v = float(h)
        if v > 0: return v
    except ValueError: pass
    lv = t.get('building:levels')
    try:
        if lv is not None and float(lv) > 0: return float(lv) * 3.2 + 1.0
    except ValueError: pass
    b = t.get('building', 'yes')
    rnd = (zlib.crc32(str(oid).encode()) % 1000) / 1000
    if b in ('church', 'cathedral'): return 22
    if b == 'chapel': return 11
    if b == 'monastery': return 15
    if b in ('garage', 'garages', 'shed', 'roof', 'kiosk', 'hut', 'carport'): return 3.2
    if b in ('house', 'detached', 'semidetached_house', 'terrace'): return 6.5 + rnd * 3
    return 11 + rnd * 7


MODERN = {'commercial', 'retail', 'industrial', 'office', 'university', 'hospital', 'school', 'college',
          'public', 'civic', 'parking', 'warehouse', 'sports_hall', 'government', 'train_station', 'hotel'}


def kind_of(t, h):
    b = t.get('building', 'yes')
    if b in ('church', 'cathedral', 'chapel', 'monastery') or t.get('amenity') == 'place_of_worship': return 1
    if t.get('historic') or 'heritage' in t or 'ref:mhs' in t or 'wikidata' in t and t.get('name'): return 2
    if b in ('roof', 'garage', 'garages', 'shed', 'carport', 'hut'): return 3
    sd = t.get('start_date', '')[:4]
    if sd.isdigit():
        y = int(sd)
        if y < 1900: return 2
        if 1944 <= y <= 1965: return 0
        if y > 1965: return 4
    if b in MODERN or h > 22: return 4
    if b in ('house', 'detached', 'semidetached_house', 'terrace'): return 5
    return 0


ROOF = {'flat': 0, 'hipped': 1, 'gabled': 2, 'mansard': 3, 'pyramidal': 1, 'half-hipped': 1, 'gambrel': 3}

out = {'b': [], 'r': [], 'w': [], 'p': [], 'wood': [], 'walls': [], 'hedges': [], 'tram': [], 'labels': [],
       'spires': [], 'castle': [flat(c) for c in castle], 'mound': MOUND, 'trees': [], 'lamps': [], 'benches': [],
       'boats': [], 'piers': []}
named = {}
bpolys = []  # pour éviter de planter des arbres dans les bâtiments


def add_building(outer, holes, t, oid):
    if area(outer) < 4: return
    if not all(inside_box(p, 1) for p in outer): return
    h = height_of(t, oid)
    mh = 0.0
    try: mh = float(t.get('min_height', '0').replace('m', ''))
    except ValueError: pass
    b = t.get('building', 'yes')
    if b == 'roof' and mh == 0: mh = max(0, h - 1.2)
    k = kind_of(t, h)
    c = centroid(outer)
    base = MOUND if on_mound(c) else 0
    # toit
    A, box = obb(outer)
    rect = area(outer) / A if A > 0 else 0
    rs = t.get('roof:shape')
    roof = 0
    rnd = zlib.crc32(str(oid).encode()) % 100
    if rs in ROOF: roof = ROOF[rs]
    elif rect > 0.86 and mh == 0 and box[4] < 26 and box[4] > 4:
        if k == 1: roof = 2
        elif k == 5: roof = 2 if rnd < 55 else 1
        elif k == 2: roof = 2 if rnd < 60 else 1
        elif k == 0: roof = 3 if rnd < 55 else 1
        elif k == 4: roof = 0
    if roof and (rect < 0.75 or box[4] > 40): roof = 0
    out['b'].append([round(h, 1), k, round(mh, 1), base, roof, box if roof else 0, flat(outer)] + [flat(hh) for hh in holes])
    bpolys.append(outer)
    if t.get('name'): named[t['name']] = (outer, h, t)


ROADW = {'motorway': 14, 'trunk': 12, 'primary': 11, 'secondary': 9, 'tertiary': 8, 'unclassified': 6,
         'residential': 6, 'living_street': 5, 'service': 4, 'pedestrian': 6, 'footway': 2.5, 'path': 2,
         'cycleway': 2, 'steps': 2.5, 'primary_link': 7, 'secondary_link': 6, 'tertiary_link': 6, 'trunk_link': 8}
road_segs = []
tram_ways = []
park_polys, wood_polys, garden_polys = [], [], []
tree_rows = []

for w in ways.values():
    t = w.get('tags', {})
    if not t: continue
    if 'building' in t:
        r = ring_of(w)
        if len(r) >= 3: add_building(r, [], t, w['id'])
    elif 'highway' in t:
        hw = t['highway']
        if hw not in ROADW: continue
        cls = 0 if hw in ('footway', 'path', 'cycleway', 'steps') else (2 if hw in ('pedestrian', 'living_street') else 1)
        if t.get('area') == 'yes' and hw == 'pedestrian':
            r = clip_poly(ring_of(w))
            if len(r) >= 3: out['p'].append([2, flat(r)])
            continue
        for piece in clip_line(pts_of(w)):
            out['r'].append([ROADW[hw], cls, flat(piece)])
            for i in range(len(piece) - 1): road_segs.append((piece[i], piece[i + 1], ROADW[hw] / 2))
    elif t.get('railway') == 'tram':
        tram_ways.append(w)
    elif t.get('natural') == 'water' or t.get('waterway') == 'riverbank':
        r = clip_poly(ring_of(w))
        if len(r) >= 3: out['w'].append([t.get('name') or t.get('water', ''), flat(r)])
    elif t.get('leisure') in ('park', 'garden') or t.get('landuse') in ('grass', 'recreation_ground', 'meadow', 'village_green'):
        r = clip_poly(ring_of(w))
        if len(r) >= 3:
            out['p'].append([0, flat(r)])
            if t.get('leisure') == 'park' or t.get('landuse') in ('recreation_ground', 'village_green'): park_polys.append(r)
            elif t.get('leisure') == 'garden': garden_polys.append(r)
    elif t.get('natural') in ('wood', 'scrub') or t.get('landuse') == 'forest':
        r = clip_poly(ring_of(w))
        if len(r) >= 3: out['p'].append([1, flat(r)]); wood_polys.append(r)
    elif t.get('natural') == 'tree_row':
        tree_rows.append(pts_of(w))
    elif t.get('barrier') == 'wall':
        for piece in clip_line(pts_of(w)):
            out['walls'].append(flat(piece))
    elif t.get('barrier') == 'hedge':
        for piece in clip_line(pts_of(w)):
            out['hedges'].append(flat(piece))
    elif t.get('man_made') == 'pier':
        for piece in clip_line(pts_of(w)):
            out['piers'].append(flat(piece))

for rel in rels.values():
    t = rel.get('tags', {})
    outers = [m['ref'] for m in rel['members'] if m['type'] == 'way' and m.get('role') in ('outer', '')]
    inners = [m['ref'] for m in rel['members'] if m['type'] == 'way' and m.get('role') == 'inner']
    if 'building' in t:
        ors, irs = join_rings(outers), join_rings(inners)
        for o in ors: add_building(o, irs, t, rel['id'])
    elif t.get('natural') == 'water' or t.get('waterway') == 'riverbank':
        ors, irs = join_rings(outers), join_rings(inners)
        for o in ors:
            c = clip_poly(o)
            if len(c) >= 3:
                out['w'].append([t.get('name') or t.get('water', ''), flat(c)] + [flat(clip_poly(i)) for i in irs if len(clip_poly(i)) >= 3])

# --- tram : chaînage des voies en lignes continues ---
segs = [pts_of(w) for w in tram_ways]
segs = [s for s in segs if len(s) > 1]
for s in segs:
    for piece in clip_line(s): out['tram'].append(flat(piece))
chains = []
pool = [list(s) for s in segs]
while pool:
    cur = pool.pop(0)
    grown = True
    while grown:
        grown = False
        for k, s in enumerate(pool):
            if s[0] == cur[-1]: cur += s[1:]
            elif s[-1] == cur[-1]: cur += s[::-1][1:]
            elif s[-1] == cur[0]: cur = s + cur[1:]
            elif s[0] == cur[0]: cur = s[::-1] + cur[1:]
            else: continue
            pool.pop(k); grown = True; break
    chains.append(cur)
def plen(p): return sum(math.dist(p[i], p[i + 1]) for i in range(len(p) - 1))
runs = []
for c in sorted(chains, key=plen, reverse=True)[:6]:
    for piece in clip_line(c):
        if plen(piece) > 400: runs.append(flat(piece))
out['tramRuns'] = runs[:3]

# --- grille spatiale bâtiments / routes ---
CELL = 25
grid = {}
for i, poly in enumerate(bpolys):
    xs = [p[0] for p in poly]; zs = [p[1] for p in poly]
    for gx in range(int(min(xs) // CELL), int(max(xs) // CELL) + 1):
        for gz in range(int(min(zs) // CELL), int(max(zs) // CELL) + 1):
            grid.setdefault((gx, gz), []).append(i)
rgrid = {}
for i, (a, b, hw) in enumerate(road_segs):
    for gx in range(int(min(a[0], b[0]) // CELL) - 1, int(max(a[0], b[0]) // CELL) + 2):
        for gz in range(int(min(a[1], b[1]) // CELL) - 1, int(max(a[1], b[1]) // CELL) + 2):
            rgrid.setdefault((gx, gz), []).append(i)


def seg_dist(p, a, b):
    dx, dz = b[0] - a[0], b[1] - a[1]
    L = dx * dx + dz * dz
    t = 0 if L == 0 else max(0, min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dz) / L))
    return math.hypot(p[0] - a[0] - t * dx, p[1] - a[1] - t * dz)


def free(p, margin=1.5):
    if not inside_box(p, 3): return False
    key = (int(p[0] // CELL), int(p[1] // CELL))
    for i in grid.get(key, []):
        if pip(p, bpolys[i]): return False
    for i in rgrid.get(key, []):
        a, b, hw = road_segs[i]
        if seg_dist(p, a, b) < hw + margin: return False
    return True


trees = []
def add_tree(p, r, kind=0):
    base = MOUND if on_mound(p) else 0
    trees.append([round(p[0], 1), round(p[1], 1), round(r, 1), base, kind])

for e in nodes.values():
    t = e.get('tags', {})
    if t.get('natural') == 'tree':
        p = xy(e['id'])
        if inside_box(p, 2): add_tree(p, 2.6 + random.random() * 2.2, 0)
    elif t.get('highway') == 'street_lamp':
        p = xy(e['id'])
        if inside_box(p, 1): out['lamps'].append([round(p[0], 1), round(p[1], 1)])
    elif t.get('amenity') == 'bench':
        p = xy(e['id'])
        if inside_box(p, 1): out['benches'].append([round(p[0], 1), round(p[1], 1), round(random.random() * 3.14, 2)])
for row in tree_rows:
    for i in range(len(row) - 1):
        a, b = row[i], row[i + 1]
        n = max(1, int(math.dist(a, b) // 8))
        for k in range(n):
            p = (a[0] + (b[0] - a[0]) * k / n, a[1] + (b[1] - a[1]) * k / n)
            if inside_box(p, 2): add_tree(p, 2.8 + random.random() * 1.5, 0)


def scatter(polys, per_m2, rmin, rmax, kind):
    for poly in polys:
        A = area(poly)
        n = int(A * per_m2)
        xs = [p[0] for p in poly]; zs = [p[1] for p in poly]
        tries = 0
        while n > 0 and tries < n * 6 + 20:
            tries += 1
            p = (random.uniform(min(xs), max(xs)), random.uniform(min(zs), max(zs)))
            if pip(p, poly) and free(p):
                add_tree(p, random.uniform(rmin, rmax), kind); n -= 1

scatter(wood_polys, 1 / 45, 2.5, 4.5, 1)
scatter(park_polys, 1 / 260, 3, 5.5, 0)
scatter(garden_polys, 1 / 380, 2.2, 4, 0)
existing = len(trees)
out['trees'] = trees
print('trees', len(trees))

# lampadaires : générés le long des axes si OSM en a peu
if len(out['lamps']) < 300:
    for w in ways.values():
        t = w.get('tags', {})
        if t.get('highway') not in ('primary', 'secondary', 'tertiary', 'pedestrian', 'residential', 'living_street'): continue
        hw = ROADW[t['highway']] / 2
        for piece in clip_line(pts_of(w)):
            acc = 0
            for i in range(len(piece) - 1):
                a, b = piece[i], piece[i + 1]
                L = math.dist(a, b)
                if L < 0.1: continue
                nx, nz = -(b[1] - a[1]) / L, (b[0] - a[0]) / L
                d = 14 - acc
                while d < L:
                    side = 1 if int(d + acc) % 2 == 0 else -1
                    p = (a[0] + (b[0] - a[0]) * d / L + nx * (hw + 0.8) * side, a[1] + (b[1] - a[1]) * d / L + nz * (hw + 0.8) * side)
                    if inside_box(p, 1): out['lamps'].append([round(p[0], 1), round(p[1], 1)])
                    d += 28
                acc = (acc + L) % 28
print('lamps', len(out['lamps']), 'benches', len(out['benches']))

# bateaux dans le bassin Saint-Pierre
for wv in out['w']:
    if 'assin' in wv[0] or wv[0] == 'canal':
        ring = list(zip(wv[1][::2], wv[1][1::2]))
        n, tries = 0, 0
        while n < 45 and tries < 5000 and len(out['boats']) < 110:
            tries += 1
            xs = [p[0] for p in ring]; zs = [p[1] for p in ring]
            p = (random.uniform(min(xs), max(xs)), random.uniform(min(zs), max(zs)))
            if not pip(p, ring): continue
            best = None
            for i in range(len(ring)):
                d = seg_dist(p, ring[i - 1], ring[i])
                if best is None or d < best[0]: best = (d, ring[i - 1], ring[i])
            if not (4 < best[0] < 11): continue
            if any(math.dist(p, (b[0], b[1])) < 4.5 for b in out['boats']): continue
            a, b = best[1], best[2]
            ang = math.atan2(b[1] - a[1], b[0] - a[0]) + math.pi / 2
            out['boats'].append([round(p[0], 1), round(p[1], 1), round(ang, 2), round(random.uniform(7, 12), 1), random.random() < 0.6])
            n += 1
print('boats', len(out['boats']))

# --- monuments ---
LAND = {
    'Abbaye aux Dames': 'Abbaye aux Dames', 'Abbaye aux Hommes': 'Abbaye aux Hommes · Hôtel de Ville',
    'Église Saint-Étienne': 'Abbatiale Saint-Étienne', 'Église Saint-Pierre': 'Église Saint-Pierre',
    'Église Saint-Jean': 'Église Saint-Jean', 'Église Saint-Sauveur': 'Église Saint-Sauveur',
    'Église Notre-Dame de la Gloriette': 'Notre-Dame de la Gloriette', 'Lycée Malherbe': 'Lycée Malherbe',
    'Palais de Justice': 'Palais de Justice', 'Maison des Quatrans': 'Maison des Quatrans',
    "Hôtel d'Escoville": "Hôtel d'Escoville", 'Église du Vieux Saint-Sauveur': 'Vieux Saint-Sauveur',
}
for name, label in LAND.items():
    if name in named:
        o, h, t = named[name]
        c = centroid(o)
        out['labels'].append([label, round(c[0], 1), round(c[1], 1), h])


def west_towers(name, n, h, r):
    if name not in named: return
    o = named[name][0]
    xs = [p[0] for p in o]; x0, x1 = min(xs), max(xs)
    west = [p for p in o if p[0] < x0 + (x1 - x0) * 0.12]
    zs = [p[1] for p in west]; z0, z1 = min(zs), max(zs)
    x = x0 + r + 1
    if n == 1: out['spires'].append([round(x, 1), round((z0 + z1) / 2, 1), h, r, 8])
    else:
        for zz in (z0 + r + 1, z1 - r - 1): out['spires'].append([round(x, 1), round(zz, 1), h, r, 8])

west_towers('Église Saint-Pierre', 1, 72, 4.5)
west_towers('Église Saint-Étienne', 2, 80, 4.2)
west_towers('Église Saint-Jean', 1, 40, 4.5)
west_towers('Église Notre-Dame de la Gloriette', 1, 34, 3)
west_towers('Église du Vieux Saint-Sauveur', 1, 40, 3)
if 'Église abbatiale de la Trinité' in named:
    o = named['Église abbatiale de la Trinité'][0]
    xs = [p[0] for p in o]; x0 = min(xs); x1 = max(xs)
    west = [p for p in o if p[0] < x0 + (x1 - x0) * 0.12]
    zs = [p[1] for p in west]
    for zz in (min(zs) + 5, max(zs) - 5): out['spires'].append([round(x0 + 5, 1), round(zz, 1), 38, 4.5, 0])
if castle:
    c = centroid(castle[0])
    out['labels'].append(['Château de Caen', round(c[0], 1), round(c[1], 1), 18 + MOUND])

json.dump(out, open(dst, 'w', encoding='utf8'), ensure_ascii=False, separators=(',', ':'))
print({k: len(v) for k, v in out.items() if isinstance(v, list)})

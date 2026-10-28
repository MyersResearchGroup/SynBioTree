"""Build a nested similarity tree over all parts and join part text onto it.

Levels:  part leaves -> 90% cluster -> 70% supercluster -> single-linkage tree over
superclusters, built from the MinHash neighbour graph via a minimum spanning tree
(exact single linkage, no full distance matrix -- 20,648 nodes would need 1.7 GB).

Reads  : seqtree/parts_meta.tsv, clu90_cluster.tsv, clu70_cluster.tsv,
         pairs_masked.tsv  (or pairs.tsv for the unmasked version)
Writes : seqtree/tree_all_parts.nwk       54,597 leaves
         seqtree/tree_superclusters.nwk   20,648 leaves, labelled
         seqtree/leaf_annotations.tsv     per-leaf displayId / name / description / role
         seqtree/cluster_text.tsv         merged text per cluster

CAVEAT: single linkage chains. With the unmasked pairs, 66% of parts collapsed into one
component of diameter 12 hops and edge density 1.9% -- the raw dendrogram topology is
not meaningful below 70% identity. The viewer uses distance cuts instead (see the
artifact). Masking generic k-mers (03_sketch_masked.py) reduces but does not remove this.
"""
import csv, os, re, sys, collections, json, threading
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, connected_components

csv.field_size_limit(10_000_000)
D = '/home/sophia/synbiohub_rdf_export/seqtree'
PAIRS = os.path.join(D, 'pairs_masked.tsv')
H90, H70, HROOT = 0.10, 0.30, 1.0
STOP = set('a an the of and or to in for is are was were with that this from by on at as it its '
           'be been has have had not no also which can may using used use we our they their '
           'part parts sequence sequences dna bba igem'.split())

# ---- inputs (manual TSV parse: descriptions contain quotes that break csv) ----
meta = {}
with open(os.path.join(D, 'parts_meta.tsv'), encoding='utf-8') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        f = line.rstrip('\n').split('\t')
        if len(f) == len(hdr):
            meta[f[0]] = dict(zip(hdr, f))
print('parts with sequence : %d' % len(meta))

rep90, rep70 = {}, {}
for line in open(os.path.join(D, 'clu90_cluster.tsv')):
    r, m = line.rstrip('\n').split('\t')[:2]; rep90[m] = r
for line in open(os.path.join(D, 'clu70_cluster.tsv')):
    r, m = line.rstrip('\n').split('\t')[:2]; rep70[m] = r
members90 = collections.defaultdict(list)
for mem, r in rep90.items():
    members90[r].append(mem)
members70 = collections.defaultdict(list)
for mem, r in rep70.items():
    members70[r].append(mem)

sf = sorted(members70)
idx = {s: i for i, s in enumerate(sf)}
n = len(sf)

best = {}
for line in open(PAIRS):
    f = line.rstrip('\n').split('\t')
    if len(f) < 3:
        continue
    qi, ti = idx.get(f[0]), idx.get(f[1])
    if qi is None or ti is None or qi == ti:
        continue
    d = float(f[2])
    key = (qi, ti) if qi < ti else (ti, qi)
    if key not in best or d < best[key]:
        best[key] = d
print('distinct pairs      : %d' % len(best))

rows = np.fromiter((k[0] for k in best), np.int32, len(best))
cols = np.fromiter((k[1] for k in best), np.int32, len(best))
vals = np.maximum(np.fromiter(best.values(), np.float64, len(best)), 1e-6)
g = coo_matrix((vals, (rows, cols)), shape=(n, n)).tocsr()
ncomp, _ = connected_components(g, directed=False)
mst = minimum_spanning_tree(g).tocoo()
edges = sorted(zip(mst.data, mst.row, mst.col))
print('components %d, MST edges %d' % (ncomp, len(edges)))

# ---- single linkage from the MST (Kruskal + union-find) ----
parent = list(range(n)); node_of = list(range(n))
children, height = {}, {}
nid = n


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]; x = parent[x]
    return x


for w, a, b in edges:
    ra, rb = find(int(a)), find(int(b))
    if ra == rb:
        continue
    children[nid] = (node_of[ra], node_of[rb])
    height[nid] = max(H70, float(w))
    parent[rb] = ra; node_of[ra] = nid; nid += 1

roots = sorted({node_of[find(i)] for i in range(n)})
print('disconnected families: %d' % len(roots))
level = list(roots)                      # balanced join -> shallow root, not a 13k chain
while len(level) > 1:
    nxt = []
    for i in range(0, len(level) - 1, 2):
        children[nid] = (level[i], level[i + 1]); height[nid] = HROOT
        nxt.append(nid); nid += 1
    if len(level) % 2:
        nxt.append(level[-1])
    level = nxt
root = level[0]

WORD = re.compile(r"[A-Za-z][A-Za-z'-]{2,}")


def label_for(ids):
    titles, wc = collections.Counter(), collections.Counter()
    for s in ids:
        m = meta.get(s)
        if not m:
            continue
        t = (m['title'] or m['displayId']).strip()
        if t:
            titles[t] += 1
        wc.update(w for w in WORD.findall((m['title'] + ' ' + m['description']).lower())
                  if w not in STOP)
    return (titles.most_common(1)[0][0] if titles else 'unnamed',
            '|'.join(w for w, _ in wc.most_common(3)))


def sanitize(s):
    return (re.sub(r"[\s(),:;\[\]'\"]+", '_', s)[:60]).strip('_') or 'x'


def nw_parts(rep, parent_h):
    mem = members90[rep]
    leaves = []
    for s in mem:
        m = meta.get(s)
        if m is None:
            continue
        leaves.append('%s:%.4f' % (sanitize('%s_%s' % (s, m['displayId'] or m['title'])), H90))
    if not leaves:
        return '%s:%.4f' % (sanitize(rep), parent_h)
    if len(leaves) == 1:
        return leaves[0].rsplit(':', 1)[0] + ':%.4f' % parent_h
    lab, _ = label_for(mem)
    return '(%s)%s:%.4f' % (','.join(leaves), sanitize('c90_%s_n%d' % (lab, len(mem))),
                            max(parent_h - H90, 0.0))


def nw_super(sfrep, parent_h):
    kids = members70[sfrep]
    if len(kids) == 1:
        return nw_parts(kids[0], parent_h)
    inner = ','.join(nw_parts(k, H70) for k in kids)
    allmem = [s for k in kids for s in members90[k]]
    lab, _ = label_for(allmem)
    return '(%s)%s:%.4f' % (inner, sanitize('c70_%s_n%d' % (lab, len(allmem))),
                            max(parent_h - H70, 0.0))


def nw_node(node, parent_h):
    if node < n:
        return nw_super(sf[node], parent_h)
    h = height[node]
    a, b = children[node]
    return '(%s,%s):%.4f' % (nw_node(a, h), nw_node(b, h), max(parent_h - h, 0.0))


def nw_node_sf(node, parent_h):
    if node < n:
        allmem = [s for k in members70[sf[node]] for s in members90[k]]
        lab, _ = label_for(allmem)
        return '%s:%.4f' % (sanitize('%s_%s_n%d' % (sf[node], lab, len(allmem))),
                            max(parent_h - H70, 0.0))
    h = height[node]
    a, b = children[node]
    return '(%s,%s):%.4f' % (nw_node_sf(a, h), nw_node_sf(b, h), max(parent_h - h, 0.0))


sys.setrecursionlimit(3000000)
threading.stack_size(512 * 1024 * 1024)       # deep recursion: chains can be thousands long
for fn, out in ((nw_node, 'tree_all_parts.nwk'), (nw_node_sf, 'tree_superclusters.nwk')):
    def emit(fn=fn, out=out):
        with open(os.path.join(D, out), 'w') as fh:
            fh.write(fn(root, HROOT) + ';\n')
    t = threading.Thread(target=emit); t.start(); t.join()

with open(os.path.join(D, 'leaf_annotations.tsv'), 'w', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['leaf', 'seqid', 'displayId', 'name', 'description', 'roles', 'types',
                'graph', 'length', 'cluster90', 'cluster90_size', 'supercluster70',
                'supercluster70_size', 'cd_uri'])
    for s, m in meta.items():
        r90 = rep90.get(s, s); r70 = rep70.get(r90, r90)
        w.writerow([sanitize('%s_%s' % (s, m['displayId'] or m['title'])), s, m['displayId'],
                    m['title'], m['description'], m['roles'], m['types'], m['graph'],
                    m['length'], r90, len(members90.get(r90, [])), r70,
                    sum(len(members90[k]) for k in members70.get(r70, [])), m['cd_uri']])

with open(os.path.join(D, 'cluster_text.tsv'), 'w', encoding='utf-8') as fh:
    w = csv.writer(fh, delimiter='\t')
    w.writerow(['level', 'cluster', 'n_parts', 'label', 'keywords', 'roles',
                'display_ids', 'merged_description'])
    for rep, mem in members90.items():
        lab, kw = label_for(mem)
        roles = collections.Counter(x for s in mem if s in meta
                                    for x in (meta[s]['roles'].split(';') if meta[s]['roles'] else []))
        descs, seen = [], set()
        for s in mem:
            if s not in meta:
                continue
            d = meta[s]['description'].strip()
            if d and d.lower() not in seen:
                seen.add(d.lower()); descs.append(d)
        w.writerow(['90', rep, len(mem), lab, kw,
                    ';'.join('%s:%d' % kv for kv in roles.most_common(5)),
                    ';'.join(sorted({meta[s]['displayId'] for s in mem
                                     if s in meta and meta[s]['displayId']}))[:500],
                    ' || '.join(descs)[:2000]])

json.dump({'parts': len(meta), 'clusters90': len(members90), 'superclusters70': n,
           'pairs': len(best), 'components': int(ncomp)},
          open(os.path.join(D, 'tree_stats.json'), 'w'), indent=1)
print('done')

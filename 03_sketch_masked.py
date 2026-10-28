"""MinHash sketches with generic-backbone k-mers masked out, then candidate pairs.

A k-mer carried by more than DF_MAX representatives is treated as shared scaffolding
(B0034-type RBS, B0010/B0012/B0015 terminators, BioBrick prefix/suffix scars, assembly
linkers) and removed before sketching, so family membership reflects the distinctive
part of each sequence instead of the glue between parts.

Reads  : seqtree/clu70_rep_seq.fasta   (output of 02_cluster.sh)
Writes : seqtree/pairs_masked.tsv      repA, repB, 1 - jaccard

Masking cut the candidate pairs from 1,391,278 to 211,637 and the largest family
from 13,477 parts to 3,704 (which is a single 70% supercluster, i.e. real homology).
Drop the masking step (use the raw sketch) by setting DF_MAX very high.
"""
import numpy as np, os, collections, time

D = '/home/sophia/synbiohub_rdf_export/seqtree'
FA = os.path.join(D, 'clu70_rep_seq.fasta')
K = 21             # canonical k-mer size
SKETCH = 256       # bottom-k sketch size
BAND = 48          # sequences sharing one of their 48 smallest hashes become candidates
MAXBUCKET = 400    # skip degenerate buckets (low-complexity k-mers shared by everything)
DF_MAX = 100       # k-mer in >100 of 20,648 reps (0.5%) counts as generic glue

CODE = np.full(256, 255, dtype=np.uint8)
for i, c in enumerate('ACGT'):
    CODE[ord(c)] = i
    CODE[ord(c.lower())] = i


def kmers(seq):
    """Canonical k-mer hashes (uint64, splitmix64-finalised) for one sequence."""
    a = CODE[np.frombuffer(seq.encode(), dtype=np.uint8)]
    if a.size < K:
        return np.empty(0, dtype=np.uint64)
    bad = np.flatnonzero(a == 255)
    segs, start = [], 0
    for b in bad:
        if b - start >= K:
            segs.append(a[start:b])
        start = b + 1
    if a.size - start >= K:
        segs.append(a[start:])
    out = []
    shifts = np.arange(K - 1, -1, -1, dtype=np.uint64) * np.uint64(2)
    for s in segs:
        s = s.astype(np.uint64)
        idx = np.arange(s.size - K + 1)[:, None] + np.arange(K)[None, :]
        w = s[idx]
        fwd = (w << shifts).sum(axis=1, dtype=np.uint64)
        rev = ((np.uint64(3) - w)[:, ::-1] << shifts).sum(axis=1, dtype=np.uint64)
        out.append(np.minimum(fwd, rev))
    if not out:
        return np.empty(0, dtype=np.uint64)
    x = np.concatenate(out) * np.uint64(0x9E3779B97F4A7C15)
    x ^= x >> np.uint64(30); x *= np.uint64(0xBF58476D1CE4E5B9)
    x ^= x >> np.uint64(27); x *= np.uint64(0x94D049BB133111EB)
    x ^= x >> np.uint64(31)
    return np.unique(x)


names, seqs, cur = [], [], []
for line in open(FA):
    if line[0] == '>':
        if cur:
            seqs.append(''.join(cur)); cur = []
        names.append(line[1:].strip().split()[0])
    else:
        cur.append(line.strip())
if cur:
    seqs.append(''.join(cur))
N = len(names)
print('reps: %d' % N, flush=True)

t0 = time.time()
km = [kmers(s) for s in seqs]
print('distinct k-mers per seq summed: %.1fM (%.0fs)'
      % (sum(h.size for h in km) / 1e6, time.time() - t0), flush=True)

allk = np.concatenate([h for h in km if h.size])
uniq, cnt = np.unique(allk, return_counts=True)
del allk
mask = np.sort(uniq[cnt > DF_MAX])
print('distinct k-mers: %.1fM | generic (DF>%d): %d (%.3f%%)'
      % (uniq.size / 1e6, DF_MAX, mask.size, 100 * mask.size / uniq.size), flush=True)
del uniq, cnt

sk, emptied = [], 0
for h in km:
    if h.size and mask.size:
        pos = np.searchsorted(mask, h)
        pos[pos >= mask.size] = mask.size - 1
        keep = h[mask[pos] != h]
    else:
        keep = h
    if keep.size == 0 and h.size:
        emptied += 1
    sk.append(keep[:SKETCH])
print('reps left with no distinctive k-mer: %d' % emptied, flush=True)

inv = collections.defaultdict(list)
for i, h in enumerate(sk):
    for v in h[:BAND]:
        inv[int(v)].append(i)
cand = set()
for v, ids in inv.items():
    if len(ids) < 2 or len(ids) > MAXBUCKET:
        continue
    for a in range(len(ids)):
        for b in range(a + 1, len(ids)):
            cand.add((ids[a], ids[b]) if ids[a] < ids[b] else (ids[b], ids[a]))
print('candidate pairs: %d' % len(cand), flush=True)

sets = [set(int(x) for x in h) for h in sk]
out = open(os.path.join(D, 'pairs_masked.tsv'), 'w')
kept = 0
for a, b in cand:
    A, B = sets[a], sets[b]
    if not A or not B:
        continue
    union = sorted(A | B)[:SKETCH]
    us = set(union)
    inter = sum(1 for v in union if v in A and v in B)
    j = inter / len(us) if us else 0.0
    if j <= 0:
        continue
    out.write('%s\t%s\t%.6f\n' % (names[a], names[b], max(1e-6, 1.0 - j)))
    kept += 1
out.close()
print('pairs written: %d -> pairs_masked.tsv' % kept)

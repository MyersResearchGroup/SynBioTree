"""Second-level grouping: merge the 95%-identity clusters into gene families via the
80%-identity clustering, so sequence variants of the same gene (with or without tags,
alternative start codons, truncations) end up under one name.

Reads  : seqtree/canonical_names.tsv, seqtree/pclu80_cluster.tsv, seqtree/parts_meta.tsv
Writes : seqtree/gene_families.tsv     one row per family
         seqtree/canonical_names.tsv   rewritten with gene_family / family_name /
                                       family_role appended
         seqtree/alias_lookup.tsv      rewritten with family_name as a second target

Result: 24,867 clusters -> 22,940 families, 1,100 of which merge more than one cluster.
Biggest: tetR (22 clusters, 66 parts, 64 names), mCherry (9 clusters, 105 names),
lacI (14 clusters, 34 parts, 26 names).
Four separate LacI CDS families remain (1,080 / 1,083 / 1,086-1,193 / 1,094-1,319 bp);
80% identity at 80% bidirectional coverage does not merge them, which is correct --
they are distinguishable sequence versions.
"""
import os, re, csv, collections, json

D = '/home/sophia/synbiohub_rdf_export/seqtree'
csv.field_size_limit(10_000_000)
VAGUE = {'sequence_feature', 'engineered_region', 'sequence', 'gene_fragment', 'unannotated'}
IDLIKE = re.compile(r'^(BBa_[A-Za-z]?\d+[\w.]*|BO_\d+|ACS_\d+|UU\d+|p\d{5,}|'
                    r'[A-Za-z]{0,3}\d{4,}|CDS|RBS|part\d*|Part_\w+|ComponentDefinition\w*)$', re.I)
norm = lambda s: re.sub(r'[^a-z0-9]', '', s.lower())

clusters = []
with open(os.path.join(D, 'canonical_names.tsv'), encoding='utf-8') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        v = line.rstrip('\n').split('\t')
        if len(v) == len(hdr):
            clusters.append(dict(zip(hdr, v)))

c80 = {}
for line in open(os.path.join(D, 'pclu80_cluster.tsv')):
    r, m = line.rstrip('\n').split('\t')[:2]
    c80[m] = r

lens = {}
with open(os.path.join(D, 'parts_meta.tsv'), encoding='utf-8') as fh:
    h = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        v = line.rstrip('\n').split('\t')
        if len(v) == len(h):
            lens[v[0]] = int(v[h.index('length')] or 0)

fam_of = {}
for c in clusters:
    votes = collections.Counter(c80.get(s) for s in c['seqids'].split() if c80.get(s))
    fam_of[c['cluster']] = votes.most_common(1)[0][0] if votes else c['cluster']

byfam = collections.defaultdict(list)
for c in clusters:
    byfam[fam_of[c['cluster']]].append(c)


def family_name(members):
    """Most used descriptive name across the family; shortest surface form breaks ties."""
    votes = collections.defaultdict(collections.Counter)
    for c in members:
        n = int(c['n_parts'])
        for nm in [c['canonical_name']] + [x.strip() for x in c['aliases'].split(';') if x.strip()]:
            if nm and not IDLIKE.match(nm):
                votes[norm(nm)][nm] += n if nm == c['canonical_name'] else 1
    if not votes:
        return members[0]['canonical_name']
    best = max(votes.values(), key=lambda g: (sum(g.values()), -min(len(x) for x in g)))
    return sorted(best.items(), key=lambda kv: (-kv[1], len(kv[0])))[0][0]


fams = []
for fid, mem in byfam.items():
    parts = sum(int(c['n_parts']) for c in mem)
    rolec = collections.Counter()
    for c in mem:
        if c['consensus_role'] not in VAGUE:
            rolec[c['consensus_role']] += int(c['n_parts'])
    if not rolec:
        rolec = collections.Counter({c['consensus_role']: int(c['n_parts']) for c in mem})
    role, rn = rolec.most_common(1)[0]
    seqs = [s for c in mem for s in c['seqids'].split()]
    L = sorted(lens.get(s, 0) for s in seqs if lens.get(s))
    fname = family_name(mem)
    seen, names = {norm(fname)}, []
    for c in mem:
        for nm in ([c['canonical_name']]
                   + [x.strip() for x in c['aliases'].split(';') if x.strip()]
                   + [x.strip() for x in c['case_variants'].split(';') if x.strip()]):
            if nm and norm(nm) not in seen:
                seen.add(norm(nm)); names.append(nm)
    flags = []
    if len(rolec) > 1:
        flags.append('mixed_roles')
    if L and L[-1] > 2 * L[0]:
        flags.append('length_range_%d-%d' % (L[0], L[-1]))
    fams.append({
        'gene_family': fid, 'family_name': fname, 'family_role': role,
        'role_support': round(rn / max(1, parts), 2),
        'n_clusters': len(mem), 'n_parts': parts,
        'bp_min': L[0] if L else 0, 'bp_median': L[len(L) // 2] if L else 0,
        'bp_max': L[-1] if L else 0,
        'n_names': len(names) + 1, 'other_names': '; '.join(names[:40]),
        'member_clusters': ' '.join(c['cluster'] for c in mem[:40]),
        'flags': ';'.join(flags),
    })

fams.sort(key=lambda f: (-f['n_parts'], -f['n_clusters']))
cols = ['gene_family', 'family_name', 'family_role', 'role_support', 'n_clusters', 'n_parts',
        'bp_min', 'bp_median', 'bp_max', 'n_names', 'other_names', 'member_clusters', 'flags']
with open(os.path.join(D, 'gene_families.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('\t'.join(cols) + '\n')
    for f in fams:
        fh.write('\t'.join(str(f[c]).replace('\t', ' ') for c in cols) + '\n')

lookup = {f['gene_family']: f for f in fams}
newcols = hdr + ['gene_family', 'family_name', 'family_role']
with open(os.path.join(D, 'canonical_names.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('\t'.join(newcols) + '\n')
    for c in clusters:
        f = lookup[fam_of[c['cluster']]]
        fh.write('\t'.join([c[k] for k in hdr]
                           + [f['gene_family'], f['family_name'], f['family_role']]) + '\n')

rows = 0
with open(os.path.join(D, 'alias_lookup.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('alias\tcanonical_name\tfamily_name\trole\tcluster\tgene_family\tkind\n')
    for c in clusters:
        f = lookup[fam_of[c['cluster']]]
        base = (c['canonical_name'], f['family_name'], c['consensus_role'],
                c['cluster'], f['gene_family'])
        fh.write('%s\t%s\t%s\t%s\t%s\t%s\tcanonical\n' % ((c['canonical_name'],) + base))
        rows += 1
        for kind, blob in (('alias', c['aliases']), ('case_variant', c['case_variants'])):
            for a in [x.strip() for x in blob.split(';') if x.strip()]:
                fh.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ((a,) + base + (kind,)))
                rows += 1

print('95%% clusters           : %d' % len(clusters))
print('gene families (80%% id) : %d' % len(fams))
print('  merging >1 cluster   : %d' % sum(1 for f in fams if f['n_clusters'] > 1))
print('  flagged mixed_roles  : %d' % sum(1 for f in fams if 'mixed_roles' in f['flags']))
print('alias lookup rows      : %d' % rows)
for f in fams[:10]:
    print('  %-20s %-12s clusters=%-3d parts=%-4d bp %d-%d names=%d'
          % (f['family_name'][:20], f['family_role'], f['n_clusters'], f['n_parts'],
             f['bp_min'], f['bp_max'], f['n_names']))
json.dump({'clusters': len(clusters), 'families': len(fams), 'lookup_rows': rows},
          open(os.path.join(D, 'family_stats.json'), 'w'))

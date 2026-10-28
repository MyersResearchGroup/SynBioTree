"""Sequence-first identity table for atomic parts.

The sequence cluster is the unit of identity. Names and roles are both derived from it:
  * canonical name  = most frequent descriptive title in the cluster (identifiers such as
                      BBa_K123456 are demoted to aliases)
  * consensus role  = majority SO role in the cluster, falling back to the 80%-identity
                      group when the 95% cluster carries no role of its own
A member whose annotated role disagrees with its cluster is reported as a suspected
mislabel rather than being allowed to split the cluster -- sequence is objective,
the role annotation is not.

Names come ONLY from the parts' own dcterms:title and sbol:displayId, never from
surrounding composite text (that is context, not identity).

Writes : seqtree/canonical_names.tsv   one row per 95% cluster
         seqtree/alias_lookup.tsv      flat alias -> canonical mapping for a search index
         seqtree/role_corrections.tsv  parts whose role contradicts their sequence group

Result: 24,867 clusters, 16,312 with at least one alias, 786 suspected role problems
(298 directly contradicting the cluster consensus).
"""
import os, re, csv, json, collections

D = '/home/sophia/synbiohub_rdf_export/seqtree'
csv.field_size_limit(10_000_000)
SONAME = {
    'SO:0000167': 'promoter', 'SO:0000316': 'CDS', 'SO:0000141': 'terminator',
    'SO:0000139': 'RBS', 'SO:0000057': 'operator', 'SO:0000112': 'primer',
    'SO:0000296': 'origin_of_replication', 'SO:0000417': 'protein_domain',
    'SO:0000031': 'aptamer', 'SO:0001687': 'restriction_site', 'SO:0000804': 'engineered_region',
    'SO:0000110': 'sequence_feature', 'SO:0000143': 'assembly_region', 'SO:0000318': 'start_codon',
    'SO:0000553': 'polyA_site', 'SO:0000997': 'gene_fragment', 'SO:0000834': 'mature_transcript',
    'SO:0000155': 'plasmid', 'SO:0000001': 'sequence', 'SO:0000324': 'tag', 'SO:0000035': 'riboswitch',
}
VAGUE = {'sequence_feature', 'engineered_region', 'sequence', 'gene_fragment'}
IDLIKE = re.compile(r'^(BBa_[A-Za-z]?\d+[\w.]*|BO_\d+|ACS_\d+|UU\d+|p\d{5,}|'
                    r'[A-Za-z]{0,3}\d{4,}|CDS|RBS|part\d*)$', re.I)

meta = {}
with open(os.path.join(D, 'parts_meta.tsv'), encoding='utf-8') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        v = line.rstrip('\n').split('\t')
        if len(v) == len(hdr):
            meta[v[0]] = dict(zip(hdr, v))

tiers, roles_of = {}, {}
with open(os.path.join(D, 'tiers.tsv'), encoding='utf-8') as fh:
    fh.readline()
    for line in fh:
        v = line.rstrip('\n').split('\t')
        tiers[v[0]] = v[2]
        roles_of[v[0]] = [SONAME.get(x, x) for x in (v[7] if len(v) > 7 else '').split(';')
                          if x.startswith('SO:')]


def load_clusters(path):
    m, of = collections.defaultdict(list), {}
    for line in open(path):
        r, x = line.rstrip('\n').split('\t')[:2]
        m[r].append(x); of[x] = r
    return m, of


m95, c95 = load_clusters(os.path.join(D, 'pclu95_cluster.tsv'))
m80, c80 = load_clusters(os.path.join(D, 'pclu80_cluster.tsv'))

host_count = collections.Counter()
with open(os.path.join(D, 'part_text.tsv'), encoding='utf-8') as fh:
    h = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        v = line.rstrip('\n').split('\t')
        if len(v) == len(h):
            host_count[v[0]] = int(v[4] or 0)


def vote(members):
    """Majority SO role over a set of parts; specific roles beat vague ones."""
    c = collections.Counter()
    for m in members:
        for r in roles_of.get(m, ()):
            c[r] += 1
    if not c:
        return None, 0, 0.0, c
    specific = collections.Counter({k: v for k, v in c.items() if k not in VAGUE})
    use = specific if specific else c
    role, n = use.most_common(1)[0]
    return role, n, n / sum(use.values()), c


norm = lambda s: re.sub(r'[^a-z0-9]', '', s.lower())
rows, corrections = [], []

for cl, mem_all in m95.items():
    mem = [m for m in mem_all if m in meta and tiers.get(m) == 'part']
    if not mem:
        continue
    desc, ids, lens = collections.Counter(), collections.Counter(), []
    for m in mem:
        for s in ((meta[m]['title'] or '').strip(), (meta[m]['displayId'] or '').strip()):
            if s:
                (ids if IDLIKE.match(s) else desc)[s] += 1
        lens.append(int(meta[m]['length'] or 0))

    if desc:
        groups = collections.defaultdict(collections.Counter)
        for n, c in desc.items():
            groups[norm(n)][n] += c
        best = max(groups.values(), key=lambda g: (sum(g.values()), -min(len(x) for x in g)))
        canon = sorted(best.items(), key=lambda kv: (-kv[1], len(kv[0])))[0][0]
    elif ids:
        canon = sorted(ids.items(), key=lambda kv: (-kv[1], len(kv[0])))[0][0]
    else:
        canon = cl

    seen, aliases, cases = {norm(canon)}, [], []
    for n in list(desc) + list(ids):
        if n == canon:
            continue
        if norm(n) in seen:
            if n not in cases:
                cases.append(n)
            continue
        seen.add(norm(n)); aliases.append(n)

    role, nvote, frac, allc = vote(mem)
    level = '95'
    if role is None:
        grp = {x for m in mem for x in m80.get(c80.get(m, ''), [])}
        role, nvote, frac, _ = vote(grp)
        level = '80' if role else '-'
    if role is None:
        role, frac, level = 'unannotated', 0.0, '-'

    for m in mem:
        own = roles_of.get(m, [])
        if not own:
            if role != 'unannotated':
                corrections.append((m, meta[m]['displayId'], '(none)', role, level,
                                    len(mem), round(frac, 2), 'missing_role'))
        elif role != 'unannotated' and role not in own:
            spec = [r for r in own if r not in VAGUE]
            corrections.append((m, meta[m]['displayId'], ';'.join(own), role, level,
                                len(mem), round(frac, 2),
                                'contradicts_cluster' if spec else 'vague_role'))

    rows.append({
        'cluster': cl, 'canonical_name': canon, 'consensus_role': role,
        'role_support': round(frac, 2), 'role_level': level,
        'roles_seen': ';'.join('%s:%d' % kv for kv in allc.most_common(4)),
        'n_parts': len(mem), 'median_bp': sorted(lens)[len(lens) // 2] if lens else 0,
        'in_composites': max((host_count.get(m, 0) for m in mem), default=0),
        'n_aliases': len(aliases), 'aliases': '; '.join(aliases[:30]),
        'case_variants': '; '.join(cases[:12]), 'seqids': ' '.join(mem[:40]),
    })

rows.sort(key=lambda r: (-r['n_parts'], -r['in_composites']))
cols = ['cluster', 'canonical_name', 'consensus_role', 'role_support', 'role_level',
        'roles_seen', 'n_parts', 'median_bp', 'in_composites', 'n_aliases', 'aliases',
        'case_variants', 'seqids']
with open(os.path.join(D, 'canonical_names.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('\t'.join(cols) + '\n')
    for r in rows:
        fh.write('\t'.join(str(r[c]).replace('\t', ' ') for c in cols) + '\n')

pairs = 0
with open(os.path.join(D, 'alias_lookup.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('alias\tcanonical_name\tconsensus_role\tcluster\tkind\n')
    for r in rows:
        fh.write('%s\t%s\t%s\t%s\tcanonical\n'
                 % (r['canonical_name'], r['canonical_name'], r['consensus_role'], r['cluster']))
        pairs += 1
        for kind, blob in (('alias', r['aliases']), ('case_variant', r['case_variants'])):
            for a in [x.strip() for x in blob.split(';') if x.strip()]:
                fh.write('%s\t%s\t%s\t%s\t%s\n'
                         % (a, r['canonical_name'], r['consensus_role'], r['cluster'], kind))
                pairs += 1

with open(os.path.join(D, 'role_corrections.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('seqid\tdisplayId\tannotated_roles\tconsensus_role\tevidence_level\t'
             'cluster_size\tconsensus_support\tkind\n')
    for c in corrections:
        fh.write('\t'.join(str(x) for x in c) + '\n')

kinds = collections.Counter(c[7] for c in corrections)
print('clusters                 : %d' % len(rows))
print('  with >=1 alias         : %d' % sum(1 for r in rows if r['n_aliases']))
print('  role from own cluster  : %d' % sum(1 for r in rows if r['role_level'] == '95'))
print('  role borrowed from 80%% : %d' % sum(1 for r in rows if r['role_level'] == '80'))
print('  still unannotated      : %d' % sum(1 for r in rows if r['consensus_role'] == 'unannotated'))
print('alias lookup rows        : %d' % pairs)
print('suspected role problems  : %d  %s' % (len(corrections), dict(kinds)))
json.dump({'clusters': len(rows), 'lookup_rows': pairs, 'corrections': len(corrections),
           'kinds': dict(kinds)}, open(os.path.join(D, 'canonical_stats.json'), 'w'))

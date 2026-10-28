"""Classify every sequence-bearing ComponentDefinition into plasmid / device / construct /
part, then compute the pan-structure (core / shell / cloud part content) of the composites.

Writes : seqtree/tiers.tsv        seqid, cd_uri, tier, n_children, n_annotations, length,
                                  displayId, roles
         seqtree/pan_device.tsv   part family -> how many devices contain it, class
         seqtree/pan_plasmid.tsv  same for plasmids
         seqtree/pan_*.json       accumulation curves + Heaps exponent

Result on the 2026-09-27 dump: part 30,885 / construct 6,067 / device 16,457 /
plasmid 1,180 / large 8.  Pan-structure of devices: 2,443 part families, core 0,
shell 6, cloud 2,437, Heaps gamma 0.562 (completely open).
"""
import gzip, glob, re, os, json, collections, random
import numpy as np

DUMP = sorted(glob.glob('/home/sophia/synbiohub_rdf_export/output*.nq.gz'))
D = '/home/sophia/synbiohub_rdf_export/seqtree'
P = 'http://sbols.org/v2#'
SO = 'http://identifiers.org/so/'
IGEM = 'http://wiki.synbiohub.org/wiki/Terms/igem#partType/'
PLASMID = {SO + 'SO:0000155', SO + 'SO:0000637', SO + 'SO:0000755',
           IGEM + 'Plasmid', IGEM + 'Plasmid_Backbone'}
CIRCULAR = SO + 'SO:0000988'
ENGREG = SO + 'SO:0000804'
DEVICEROLE = {IGEM + x for x in ('Device', 'Generator', 'Composite', 'Inverter',
                                 'Translational_Unit', 'Reporter', 'Signalling',
                                 'Measurement', 'Project', 'Intermediate')}
ATOMIC = {SO + x for x in ('SO:0000167', 'SO:0000316', 'SO:0000141', 'SO:0000139',
                           'SO:0000057', 'SO:0000112', 'SO:0000296', 'SO:0000417',
                           'SO:0000031', 'SO:0001687', 'SO:0000318', 'SO:0000553',
                           'SO:0000143', 'SO:0000997')}
SONAME = {'SO:0000167': 'promoter', 'SO:0000316': 'CDS', 'SO:0000141': 'terminator',
          'SO:0000139': 'RBS', 'SO:0000057': 'operator', 'SO:0000112': 'primer',
          'SO:0000296': 'origin', 'SO:0000417': 'protein domain', 'SO:0000031': 'aptamer',
          'SO:0001687': 'restriction site', 'SO:0000143': 'assembly region',
          'SO:0000804': 'engineered region', 'SO:0000155': 'plasmid',
          'SO:0000755': 'plasmid vector', 'SO:0000110': 'sequence feature'}

isCD = set()
comp_def, cd_components = {}, collections.defaultdict(list)
sa_count = collections.Counter()
roles, types = collections.defaultdict(set), collections.defaultdict(set)

for path in DUMP:
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 4:
                continue
            s, p, o = f[0].strip('<>'), f[1].strip('<>'), f[2]
            if p == 'http://www.w3.org/1999/02/22-rdf-syntax-ns#type':
                if o.strip('<>') == P + 'ComponentDefinition':
                    isCD.add(s)
            elif p == P + 'component':
                cd_components[s].append(o.strip('<>'))
            elif p == P + 'definition':
                comp_def[s] = o.strip('<>')
            elif p == P + 'sequenceAnnotation':
                sa_count[s] += 1
            elif p == P + 'role':
                roles[s].add(o.strip('<>'))
            elif p == P + 'type' and o.startswith('<'):
                types[s].add(o.strip('<>'))

meta = {}
with open(os.path.join(D, 'parts_meta.tsv'), encoding='utf-8') as fh:
    hdr = fh.readline().rstrip('\n').split('\t')
    for line in fh:
        v = line.rstrip('\n').split('\t')
        if len(v) == len(hdr):
            r = dict(zip(hdr, v))
            meta[r['cd_uri']] = r

children = {cd: [comp_def[c] for c in cd_components.get(cd, []) if c in comp_def] for cd in meta}


def tier(cd):
    r, t = roles.get(cd, set()), types.get(cd, set())
    L, nk, nsa = int(meta[cd]['length'] or 0), len(children[cd]), sa_count.get(cd, 0)
    if r & PLASMID or CIRCULAR in t:
        return 'plasmid'
    if nk >= 2:
        return 'device' if (r & DEVICEROLE or ENGREG in r) else 'construct'
    if nk == 1 or (nsa >= 2 and (ENGREG in r or r & DEVICEROLE)):
        return 'construct'
    if L >= 20000:
        return 'large'
    if r & ATOMIC or ENGREG not in r:
        return 'part'
    return 'construct' if L > 500 else 'part'


tiers = {cd: tier(cd) for cd in meta}
cnt = collections.Counter(tiers.values())
print('%-10s %7s %12s %10s' % ('tier', 'count', 'bp', 'median bp'))
for t in ('part', 'construct', 'device', 'plasmid', 'large'):
    ls = sorted(int(meta[cd]['length'] or 0) for cd in meta if tiers[cd] == t)
    print('%-10s %7d %12d %10d' % (t, cnt[t], sum(ls), ls[len(ls) // 2] if ls else 0))

with open(os.path.join(D, 'tiers.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('seqid\tcd_uri\ttier\tn_children\tn_annotations\tlength\tdisplayId\troles\n')
    for cd, m in meta.items():
        fh.write('%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\n' % (
            m['seqid'], cd, tiers[cd], len(children[cd]), sa_count.get(cd, 0),
            m['length'], m['displayId'],
            ';'.join(sorted(x.replace(SO, '').replace(IGEM, 'igem:') for x in roles.get(cd, ())))))

# ---------- pan-structure ----------
rep90 = {}
for line in open(os.path.join(D, 'clu90_cluster.tsv')):
    r, m = line.rstrip('\n').split('\t')[:2]
    rep90[m] = r
seq_of_cd = {m['cd_uri']: m['seqid'] for m in meta.values()}
fam_of_cd = {cd: rep90.get(seq_of_cd.get(cd, ''), seq_of_cd.get(cd, '')) for cd in meta}
label_of_fam = {}
for cd, f in fam_of_cd.items():
    label_of_fam.setdefault(f, meta[cd]['displayId'] or meta[cd]['title'])


def flatten(cd, depth=0, seen=None):
    """Atomic part families contained in a composite, recursively."""
    if seen is None:
        seen = set()
    if cd in seen or depth > 6:
        return []
    seen.add(cd)
    kids = children.get(cd, [])
    if not kids:
        return [fam_of_cd[cd]] if cd in meta and tiers.get(cd) in ('part', 'construct') else []
    out = []
    for k in kids:
        if k in meta and not children.get(k):
            if tiers.get(k) == 'part':
                out.append(fam_of_cd[k])
        else:
            out.extend(flatten(k, depth + 1, seen))
    return out


for group in ('device', 'plasmid'):
    hosts = [cd for cd in meta if tiers[cd] == group and children[cd]]
    content = {cd: set(flatten(cd)) for cd in hosts}
    content = {k: v for k, v in content.items() if v}
    if not content:
        continue
    df = collections.Counter()
    for fams in content.values():
        df.update(fams)
    n, pan = len(content), len(df)
    core = sum(1 for v in df.values() if v >= 0.95 * n)
    soft = sum(1 for v in df.values() if 0.75 * n <= v < 0.95 * n)
    shell = sum(1 for v in df.values() if 0.15 * n <= v < 0.75 * n)
    cloud = sum(1 for v in df.values() if v < 0.15 * n)
    uniq = sum(1 for v in df.values() if v == 1)
    print('\n=== pan-structure: %s ===' % group)
    print('%d %ss, pan %d | core %d soft %d shell %d cloud %d (singletons %d)'
          % (n, group, pan, core, soft, shell, cloud, uniq))
    for fam, c in df.most_common(10):
        print('    %-22s %5d  %5.1f%%' % (label_of_fam.get(fam, fam)[:22], c, 100 * c / n))

    keys = list(content)
    cps, ccs = [], []
    for _ in range(20):
        random.shuffle(keys)
        seen, inter, cp, cc = set(), None, [], []
        for k in keys:
            seen |= content[k]
            inter = set(content[k]) if inter is None else (inter & content[k])
            cp.append(len(seen)); cc.append(len(inter))
        cps.append(cp); ccs.append(cc)
    cp, cc = np.mean(cps, axis=0), np.mean(ccs, axis=0)
    x = np.arange(1, n + 1); k = x >= 10
    gamma = np.polyfit(np.log(x[k]), np.log(np.maximum(cp[k], 1)), 1)[0]
    print('  Heaps gamma = %.3f (<1 means open pan-structure)' % gamma)

    with open(os.path.join(D, 'pan_%s.tsv' % group), 'w', encoding='utf-8') as fh:
        fh.write('part_family\tlabel\tin_%ss\tfraction\tclass\n' % group)
        for fam, c in df.most_common():
            frac = c / n
            cls = ('core' if frac >= 0.95 else 'soft-core' if frac >= 0.75
                   else 'shell' if frac >= 0.15 else 'cloud')
            fh.write('%s\t%s\t%d\t%.4f\t%s\n' % (fam, label_of_fam.get(fam, fam), c, frac, cls))
    json.dump({'hosts': n, 'pan': pan, 'core': core, 'soft': soft, 'shell': shell,
               'cloud': cloud, 'singletons': uniq, 'gamma': float(gamma),
               'curve_pan': [float(v) for v in cp[::max(1, n // 200)]],
               'curve_core': [float(v) for v in cc[::max(1, n // 200)]]},
              open(os.path.join(D, 'pan_%s.json' % group), 'w'))

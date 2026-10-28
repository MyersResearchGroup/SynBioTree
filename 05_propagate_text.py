"""Push device / plasmid / construct descriptions down onto the atomic parts they contain,
and mine the context words each part sequence is surrounded by.

Writes : seqtree/part_text.tsv   one row per atomic part: own text + text inherited
                                 from every ancestor, prefixed with the ancestor's name
         seqtree/aliases.tsv     per strict cluster: names the parts carry themselves,
                                 plus gene-like tokens mined from the host composites

Only 7,948 of the 30,885 atomic parts are inside any composite, so most rows have no
inherited text.

IMPORTANT: the mined column is CONTEXT, not names. It is scored for specificity
(token rate among this cluster's hosts vs. the corpus-wide rate, >= 3x) but tokens like
"LacI" still appear next to parts that are not LacI. Use own_names + display_ids for
identity; see 06_canonical_names.py.
"""
import gzip, glob, re, os, html, json, collections, math

DUMP = sorted(glob.glob('/home/sophia/synbiohub_rdf_export/output*.nq.gz'))
D = '/home/sophia/synbiohub_rdf_export/seqtree'
P = 'http://sbols.org/v2#'
TEXTP = {
    'http://purl.org/dc/terms/title': 'title',
    'http://purl.org/dc/terms/description': 'description',
    'http://wiki.synbiohub.org/wiki/Terms/synbiohub#mutableDescription': 'longdesc',
    'http://wiki.synbiohub.org/wiki/Terms/synbiohub#mutableNotes': 'notes',
    P + 'displayId': 'displayId',
}
TAG = re.compile(r'<[^>]{1,200}>')
ESC = {'n': ' ', 't': ' ', 'r': ' ', '"': '"', '\\': '\\'}


def unesc(s):
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            nx = s[i + 1]
            if nx in ESC:
                out.append(ESC[nx]); i += 2; continue
            if nx == 'u' and i + 5 < len(s):
                try:
                    out.append(chr(int(s[i + 2:i + 6], 16))); i += 6; continue
                except ValueError:
                    pass
        out.append(c); i += 1
    return ''.join(out)


def clean(s):
    return re.sub(r'\s+', ' ', html.unescape(TAG.sub(' ', unesc(s)))).strip()


text = collections.defaultdict(dict)
cd_comp = collections.defaultdict(list)
comp_def = {}
for path in DUMP:
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            f = line.split('\t')
            if len(f) < 4:
                continue
            s, p, o = f[0].strip('<>'), f[1].strip('<>'), f[2]
            if p == P + 'component':
                cd_comp[s].append(o.strip('<>'))
            elif p == P + 'definition':
                comp_def[s] = o.strip('<>')
            elif p in TEXTP:
                e = o.rfind('"')
                if o.startswith('"') and e > 0:
                    text[s][TEXTP[p]] = clean(o[1:e])

tier, seqid_of, cd_of_seq = {}, {}, {}
with open(os.path.join(D, 'tiers.tsv'), encoding='utf-8') as fh:
    fh.readline()
    for line in fh:
        v = line.rstrip('\n').split('\t')
        tier[v[1]] = v[2]; seqid_of[v[1]] = v[0]; cd_of_seq[v[0]] = v[1]

rep = {}
for line in open(os.path.join(D, 'pclu95_cluster.tsv')):   # strict: 95% id, 90% both-way cov
    r, m = line.rstrip('\n').split('\t')[:2]
    rep[m] = r

children = {cd: [comp_def[c] for c in cd_comp.get(cd, []) if c in comp_def] for cd in tier}


def atomic_descendants(cd, depth=0, seen=None):
    if seen is None:
        seen = set()
    if cd in seen or depth > 6:
        return set()
    seen.add(cd)
    kids = children.get(cd, [])
    if not kids:
        return {cd} if tier.get(cd) == 'part' else set()
    out = set()
    for k in kids:
        if k in tier and not children.get(k):
            if tier[k] == 'part':
                out.add(k)
        else:
            out |= atomic_descendants(k, depth + 1, seen)
    return out


parents_of = collections.defaultdict(set)
composites = [cd for cd in tier if tier[cd] in ('device', 'plasmid', 'construct') and children.get(cd)]
for cd in composites:
    for p in atomic_descendants(cd):
        parents_of[p].add(cd)
print('composites %d, atomic parts with a parent %d' % (len(composites), len(parents_of)))


def alltext(cd):
    t = text.get(cd, {})
    return ' '.join(x for x in (t.get('title'), t.get('description'),
                                t.get('longdesc'), t.get('notes')) if x)


atomic = [cd for cd in tier if tier[cd] == 'part']
with open(os.path.join(D, 'part_text.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('seqid\tdisplayId\tname\tcluster95\tn_parents\town_text\tinherited_text\tparent_names\n')
    for cd in atomic:
        t = text.get(cd, {})
        ps = sorted(parents_of.get(cd, ()))
        inh, pn, seen = [], [], set()
        for p in ps:
            pt = text.get(p, {})
            nm = pt.get('title') or pt.get('displayId') or ''
            if nm:
                pn.append(nm)
            s = ' '.join(x for x in (pt.get('description'), pt.get('longdesc')) if x)
            k = s.lower()[:160]
            if s and k not in seen:
                seen.add(k); inh.append('[%s] %s' % (nm or pt.get('displayId', '?'), s))
        fh.write('\t'.join([
            seqid_of[cd], t.get('displayId', ''), t.get('title', ''),
            rep.get(seqid_of[cd], seqid_of[cd]), str(len(ps)),
            (' '.join(x for x in (t.get('description'), t.get('longdesc')) if x))[:2000].replace('\t', ' '),
            (' || '.join(inh))[:6000].replace('\t', ' '),
            (';'.join(pn[:40]))[:800].replace('\t', ' '),
        ]) + '\n')
print('part_text.tsv written (%d atomic parts)' % len(atomic))

# ---------- context mining, scored for specificity ----------
GENEISH = re.compile(r'\b(?:[a-z]{2,5}[A-Z][A-Za-z0-9]{0,4}|[A-Z][a-z]{2,4}[A-Z][A-Za-z0-9]{0,4}|'
                     r'p[A-Z][A-Za-z0-9]{2,8}|[A-Z]{2,6}\d{0,3})\b')
GENERIC = set('''DNA RNA PCR RBS CDS ORF UTR GFP RFP YFP CFP BFP SBOL BioBrick IGEM RFC PRO NOT AND OR
NB TODO NA OK USA UK MIT BBa SO EC ATP GTP NADH IPTG AHL LB OD IU PDF URL HTML XML ID IDs
E Coli DH FACS FRET ELISA SDS PAGE RT qPCR KO WT MCS SNP'''.split())

host_tokens, bg = {}, collections.Counter()
for p in composites:
    toks = {w for w in GENEISH.findall(alltext(p)) if w not in GENERIC and len(w) >= 3}
    if toks:
        host_tokens[p] = toks
        bg.update(toks)
NH = max(1, len(host_tokens))

members = collections.defaultdict(list)
for seq, r in rep.items():
    members[r].append(seq)

rows = []
for cl, mem in members.items():
    cds = [cd_of_seq[s] for s in mem if s in cd_of_seq]
    cds = [c for c in cds if tier.get(c) == 'part']
    if not cds:
        continue
    own, dids = collections.Counter(), set()
    for c in cds:
        t = text.get(c, {})
        if t.get('title'):
            own[t['title'].strip()] += 1
        if t.get('displayId'):
            dids.add(t['displayId'])
    hostset = {p for c in cds for p in parents_of.get(c, ())}
    mined = collections.Counter()
    for p in hostset:
        mined.update(host_tokens.get(p, ()))
    nh = max(1, len(hostset))
    scored = []
    for w, c in mined.items():
        if c < 2:
            continue
        enr = (c / nh) / max(bg[w] / NH, 1e-9)
        if enr >= 3.0:
            scored.append((c * math.log(enr), w, c))
    scored.sort(reverse=True)
    top = [(w, c) for _, w, c in scored[:12]]
    if not own and not top:
        continue
    rows.append({'cluster': cl, 'n_parts': len(cds), 'n_hosts': len(hostset),
                 'own': own.most_common(), 'dids': sorted(dids)[:20], 'mined': top})

rows.sort(key=lambda r: -(len(r['own']) + len(r['mined'])))
with open(os.path.join(D, 'aliases.tsv'), 'w', encoding='utf-8') as fh:
    fh.write('cluster95\tn_parts\tn_host_composites\tn_distinct_own_names\town_names\t'
             'display_ids\tmined_from_parent_text\n')
    for r in rows:
        fh.write('\t'.join([
            r['cluster'], str(r['n_parts']), str(r['n_hosts']), str(len(r['own'])),
            '; '.join('%s(%d)' % (n, c) for n, c in r['own'])[:1200],
            ' '.join(r['dids']),
            '; '.join('%s(%d)' % (w, c) for w, c in r['mined'])[:600],
        ]) + '\n')
print('aliases.tsv written (%d clusters)' % len(rows))

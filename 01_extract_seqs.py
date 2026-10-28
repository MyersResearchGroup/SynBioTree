"""Extract DNA sequences + their part metadata from the SynBioHub N-Quads dump.

Reads  : /home/sophia/synbiohub_rdf_export/output*.nq.gz
Writes : seqtree/parts.fasta        one record per ComponentDefinition that has a sequence
         seqtree/parts_meta.tsv     seqid, cd_uri, seq_uri, displayId, title, description,
                                    roles, types, graph, length
"""
import gzip, glob, re, os, html, json, collections

DUMP = sorted(glob.glob('/home/sophia/synbiohub_rdf_export/output*.nq.gz'))
OUT = '/home/sophia/synbiohub_rdf_export/seqtree'
os.makedirs(OUT, exist_ok=True)

P = 'http://sbols.org/v2#'
WANT = {
    P + 'elements': 'elements',
    P + 'sequence': 'sequence',
    P + 'displayId': 'displayId',
    P + 'role': 'role',
    P + 'type': 'type',
    'http://purl.org/dc/terms/title': 'title',
    'http://purl.org/dc/terms/description': 'description',
    'http://www.w3.org/1999/02/22-rdf-syntax-ns#type': 'rdftype',
}
TAG = re.compile(r'<[^>]{1,200}>')
ESC = {'n': ' ', 't': ' ', 'r': ' ', '"': '"', '\\': '\\'}


def unescape(s):
    out, i = [], 0
    while i < len(s):
        c = s[i]
        if c == '\\' and i + 1 < len(s):
            n = s[i + 1]
            if n in ESC:
                out.append(ESC[n]); i += 2; continue
            if n == 'u' and i + 5 < len(s):
                try:
                    out.append(chr(int(s[i + 2:i + 6], 16))); i += 6; continue
                except ValueError:
                    pass
        out.append(c); i += 1
    return ''.join(out)


def clean_text(s):
    return re.sub(r'\s+', ' ', html.unescape(TAG.sub(' ', unescape(s)))).strip()


elements = {}
cd_seq = collections.defaultdict(list)
meta = collections.defaultdict(dict)
roles = collections.defaultdict(set)
types = collections.defaultdict(set)
rdftypes = {}
graph_of = {}

for path in DUMP:
    with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
        for line in fh:
            parts = line.split('\t')
            if len(parts) < 4:
                continue
            pred = parts[1].strip('<>')
            key = WANT.get(pred)
            if key is None:
                continue
            subj = parts[0].strip('<>')
            obj = parts[2]
            graph = parts[3].strip().rstrip(' .').strip('<>')
            if key == 'elements':
                end = obj.rfind('"')
                if obj.startswith('"') and end > 0:
                    seq = re.sub(r'[^A-Za-z]', '', obj[1:end]).upper()
                    if seq:
                        elements[subj] = seq
            elif key == 'sequence':
                cd_seq[subj].append(obj.strip('<>'))
            elif key == 'role':
                roles[subj].add(obj.strip('<>'))
            elif key == 'type':
                if obj.startswith('<'):
                    types[subj].add(obj.strip('<>'))
            elif key == 'rdftype':
                rdftypes[subj] = obj.strip('<>')
                graph_of[subj] = graph
            else:
                end = obj.rfind('"')
                if obj.startswith('"') and end > 0:
                    meta[subj][key] = clean_text(obj[1:end])
                    graph_of.setdefault(subj, graph)

print('sequences with elements : %d' % len(elements))
print('CDs referencing sequence: %d' % len(cd_seq))

SO = 'http://identifiers.org/so/'
rows = []
for cd, seqs in cd_seq.items():
    if rdftypes.get(cd) != P + 'ComponentDefinition':
        continue
    for s in seqs:
        dna = elements.get(s)
        if not dna:
            continue
        m = meta.get(cd, {})
        so_roles = sorted(r.replace(SO, '') for r in roles.get(cd, ()) if r.startswith(SO))
        rows.append({
            'cd_uri': cd, 'seq_uri': s,
            'displayId': m.get('displayId', ''), 'title': m.get('title', ''),
            'description': m.get('description', ''),
            'roles': ';'.join(so_roles),
            'types': ';'.join(sorted(t.split('#')[-1] for t in types.get(cd, ()))),
            'graph': graph_of.get(cd, ''), 'length': len(dna), 'dna': dna,
        })

rows.sort(key=lambda r: r['cd_uri'])
print('rows (CD with sequence)  : %d' % len(rows))

fa = open(os.path.join(OUT, 'parts.fasta'), 'w')
tsv = open(os.path.join(OUT, 'parts_meta.tsv'), 'w', encoding='utf-8')
cols = ['seqid', 'cd_uri', 'seq_uri', 'displayId', 'title', 'description',
        'roles', 'types', 'graph', 'length']
tsv.write('\t'.join(cols) + '\n')
for i, r in enumerate(rows, 1):
    sid = 'p%06d' % i
    fa.write('>%s\n' % sid)
    dna = r['dna']
    for j in range(0, len(dna), 80):
        fa.write(dna[j:j + 80] + '\n')
    tsv.write('\t'.join([sid] + [str(r[c]).replace('\t', ' ') for c in cols[1:]]) + '\n')
fa.close(); tsv.close()

lens = sorted(r['length'] for r in rows)
uniq = len({r['dna'] for r in rows})
print('unique DNA strings       : %d (%.0f%% redundancy)' % (uniq, 100 * (1 - uniq / len(rows))))
print('length  min %d  median %d  max %d  total %.1f Mbp'
      % (lens[0], lens[len(lens) // 2], lens[-1], sum(lens) / 1e6))
json.dump({'rows': len(rows), 'unique': uniq}, open(os.path.join(OUT, 'extract_stats.json'), 'w'))

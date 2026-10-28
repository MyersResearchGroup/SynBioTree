#!/bin/bash
# mmseqs2 clustering steps. Run after 01_extract_seqs.py.
#
# Two families of clustering are used for different purposes:
#
#   clu90 / clu70   cov-mode 1 (coverage of the shorter sequence only).
#                   Used ONLY to pick representatives for the similarity tree.
#                   NOT suitable for identity: a short part contained inside a long
#                   one lands in the same cluster, which put 2,310 unrelated parts
#                   (PhlF, YFP, AmpR ...) into one group.
#
#   pclu95 / pclu80 cov-mode 0 (bidirectional coverage), atomic parts only.
#                   These are the clusters used for canonical names, aliases and
#                   role consensus. Largest cluster is 157 members.
#
set -euo pipefail
D=/home/sophia/synbiohub_rdf_export/seqtree
TMP=$(mktemp -d)
cd "$D"

# ---- representatives for the tree (all 54,597 sequences) ----
mmseqs easy-linclust parts.fasta clu90 "$TMP/t1" \
    --min-seq-id 0.9 -c 0.8 --cov-mode 1 --threads 8        # -> 23,827 clusters
mmseqs easy-linclust clu90_rep_seq.fasta clu70 "$TMP/t2" \
    --min-seq-id 0.7 -c 0.5 --cov-mode 1 --threads 8        # -> 20,648 superclusters

# ---- strict clusters for identity (atomic parts only; needs tiers.tsv) ----
# Produce parts_atomic.fasta first:
#   awk -F'\t' 'NR>1 && $3=="part" {print $1}' tiers.tsv > atomic_ids.txt
#   (then filter parts.fasta by those ids -- see README)
if [ -f parts_atomic.fasta ]; then
    mmseqs easy-linclust parts_atomic.fasta pclu95 "$TMP/t3" \
        --min-seq-id 0.95 -c 0.9 --cov-mode 0 --threads 8   # -> 24,867 clusters
    mmseqs easy-linclust parts_atomic.fasta pclu80 "$TMP/t4" \
        --min-seq-id 0.8  -c 0.8 --cov-mode 0 --threads 8   # -> 22,940 gene families
fi

rm -rf "$TMP"

# NOTE: mmseqs easy-search (all-vs-all alignment) was tried for the tree distances
# and abandoned -- it ran 8.6 hours without finishing the prefilter on this machine
# (2 GB free RAM, one 4.65 Mbp genome in the input). 03_sketch_masked.py replaces it
# and finishes in ~30 seconds.

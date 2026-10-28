# SynBioTree

Sequence clustering and similarity tree for the parts in a SynBioHub
instance, built from a full N-Quads dump of its Virtuoso store.

The scripts run in order from this directory. Each step reads the previous
step's output from the data directory.

| step | script | produces |
|---|---|---|
| 1 | `01_extract_seqs.py` | `parts.fasta`, `parts_meta.tsv` |
| 2a | `02_cluster.sh` | `clu90_*`, `clu70_*` (for the tree) |
| 3 | `03_sketch_masked.py` | `pairs_masked.tsv` |
| 4 | `04_classify_and_pan.py` | `tiers.tsv`, `pan_device.*`, `pan_plasmid.*` |
| 2b | `02_cluster.sh` on atomic parts | `pclu95_*`, `pclu80_*` |
| 5 | `05_propagate_text.py` | `part_text.tsv`, `aliases.tsv` |
| 6 | `06_canonical_names.py` | `canonical_names.tsv`, `alias_lookup.tsv`, `role_corrections.tsv` |
| 7 | `07_gene_families.py` | `gene_families.tsv` |
| 8 | `08_build_tree.py` | `tree_*.nwk`, `leaf_annotations.tsv`, `cluster_text.tsv` |

Step 4 has to run before the second half of step 2, because the atomic
parts list comes from `tiers.tsv`.

## Inputs

These scripts read `output*.nq.gz` (an N-Quads dump of the SynBioHub Virtuoso
store). The dump and every derived file are not in this repository.

## Known caveats

- `clu90` and `clu70` are run with `--cov-mode 1`, which only requires the
  shorter sequence to be covered. A short part can fall into the cluster of a
  longer sequence that contains it. Use these clusters to pick tree
  representatives, not to judge identity.
- `pclu95` and `pclu80` use `--cov-mode 0` (two-way coverage) and are the
  identity clusters.
- Distances are k-mer MinHash Jaccard estimates, not alignment identities.
  `tree_*.nwk` are hierarchical clustering dendrograms, not phylogenies.
- Roles are a majority vote inside a sequence cluster. Disagreeing parts are
  listed in `role_corrections.tsv`; clusters are not split on them.

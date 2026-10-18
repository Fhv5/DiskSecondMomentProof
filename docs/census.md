# Graph census

Inputs are included in `data/planar_graphs/`; the proof classes are in
`data/connected_graphs/`. Commands start at the repository root.

| n | Planar classes | Outerplanar classes |
|---:|---:|---:|
| 2 | 1 | 1 |
| 3 | 2 | 2 |
| 4 | 6 | 5 |
| 5 | 20 | 13 |
| 6 | 99 | 46 |

The filter validates syntax, connectedness and absence of isomorphic
duplicates. It tests circular vertex orders for an embedding without
alternating endpoints of independent edges, then compares the accepted
classes with the stored proof list up to graph isomorphism. The proof
launcher performs every selected comparison. To check an input directly:

```powershell
julia --startup-file=no --project=. scripts/filter_census.jl --n 6 --input data/planar_graphs/planar6.g6
```

Use `--n N --input data/planar_graphs/planarN.g6` for `N=2,3,4,5`.
The two-vertex class is derived directly. For four vertices the checker
also enumerates all 64 labeled graphs and verifies the K4 Gram obstruction.
The [mathematical argument](proof.md) establishes the scope of the census.

## Regeneration

With Bash, plantri, nauty's `shortg`, `mktemp` and `sha256sum` installed:

```sh
bash scripts/generate_planar_census.sh output/census
```

This runs `plantri -p -c1 -g N | shortg` for `N=3,4,5,6`, checks the counts
and writes executable/output digests to `output/census/toolchain.txt`.
Filter the generated files using the following pattern:

```powershell
julia --startup-file=no --project=. scripts/filter_census.jl --n 6 --input output/census/planar6.g6 --output output/census/graphs6.txt
```

Repeat for `N=3,4,5` with the corresponding filenames. The scripts refuse
to overwrite existing result files. Generated data stay under `output/`;
the versioned inputs remain in `data/`. Their commands, counts and SHA-256
digests are listed in [the input inventory](../data/planar_graphs/provenance.txt).

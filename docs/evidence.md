# Evidence formats

All source paths start at the repository root. `<case>` denotes
`output/releases/<run-id>/nN/` for a release or `output/nN/` for a default
single-case run. The full launcher enables evidence recording; the Julia
entry point enables it with `--record-evidence`.

## Inventory

`output/releases/<run-id>/release.json` identifies the commit, selected disk
counts, thread policy and SHA-256 digests. Its `source.zip` is `git archive`
of that commit. `<case>/run.json` records parameters, environment, graph
counts (`interval_graphs`, `leaf_graphs`, `leaf_excluded`, `isolated_boxes`),
conclusions, source digests and result-file digests. Source inventories
include documentation, tests and data. Text identities also record SHA-256
after replacing CRLF with LF; all other bytes are preserved. The run inventory
uses schema 3 and the search-conventions record uses schema 2.

## Coordinates and transcripts

`<case>/search/conventions.json` records each graph's route, anchor, oriented
tree, angle and multiplier order, reflection coordinate, rational initial
box, symmetry orders, permutations and transcript path. Indices are one-based;
angles are in radians along parent-to-child edges. Permutations map old
labels to new labels: `p[old]=new`, so `c'_(p(i))=c_i`.

`<case>/evidence/angular/` contains one TSV transcript per searched graph,
including graphs with leaves. `<case>/search/angle_boxes/nN-summary.tsv`
records the corresponding counters and completion status.
Filenames encode graph6 bytes in hexadecimal. Each line gives box depth and
decision in depth-first order. The rational initial box, widest-coordinate
bisection rule and lower-child-first order reconstruct every tested box.
The [proof](proof.md) derives the necessary conditions behind the decisions.

| Decision | Verified condition |
|---|---|
| `contact` | Required contact interval excludes zero. |
| `overlap` | A prescribed nonedge has no strictly positive gap. |
| `symmetry` | Bounds violate a justified angle order. |
| `leaf_equilibrium` | Leaf balance or force bound fails. |
| `bridge_equilibrium` | Bridge force balance fails. |
| `bridge_rotation` | Necessary rotational derivative condition fails. |
| `multiplier` | Multiplier enclosure misses required force bounds. |
| `stationarity` | A balance-residual component excludes zero. |
| `moment_balance` | `2 sum(lambda)-M` excludes zero. |
| `bridge_krawczyk` | A square subsystem of bridge-rotation equations excludes a root. |
| `krawczyk` | Residual or Krawczyk image excludes a root. |
| `boundary_root` | Uniqueness identifies an exact additional-contact witness. |
| `isolated` | Strict Krawczyk inclusion and strict noncontacts hold. |

`split` replaces a box by its two closed children. Uncertain rank, failed
verification or an operational limit cannot establish exclusion; unresolved
boxes prevent closure. Transcripts describe the computation; rerunning the
source evaluates its interval inequalities.

## Candidate analysis

`<case>/refinement/candidates.jsonl` stores root enclosures with rational
endpoints and their oriented tree: lower numerator/denominator and upper numerator/denominator,
encoded as decimal integer strings. Serialization and restoration verify
outward containment. Reported decimal displays are for inspection.

`<case>/analysis/nN-interval-decisions.jsonl` records the objective, its margin
above the target, multiplier and rank statuses, and verified negative intrinsic
values with their tangent directions. Objective exclusion requires
`inf(M-U_n)>0`; a negative-direction exclusion requires `sup(Q)<0`.
Targets and continuous families have dedicated exact-analysis records.
Each isolated box is analysed, including boxes assigned to the same geometric
class. Full per-object reports are in `<case>/analysis/evidence/`.

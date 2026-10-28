# DiskSecondMomentProof

Computational verification for *Second Moment Penny Packings via Rolling
Spaces*, by José Ayala Hoffmann and Fabián Henry Vilaxa.

For packings of unit disks, the program verifies the minimum of
`M(c) = (1/2) sum_i ||c_i-c_bar||^2` and its equality cases:

| Disks | Minimum |
|---:|---:|
| 2 | 1 |
| 3 | 2 |
| 4 | 4 |
| 5 | 34/5 |
| 6 | 29/3 |

## Setup

Use Julia **1.12.7**, Git and PowerShell. All commands start at the
repository root. Install the dependency versions fixed by `Manifest.toml`:

```powershell
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Run

From a clean Git checkout, compute all five cases at 512-bit precision:

```powershell
& .\scripts\run_complete_proof.ps1
```

The launcher checks the included graph censuses, fixes Julia, GC, BLAS and
FLINT to one thread, and creates `output/releases/<run-id>/`. This contains
`source.zip`, `release.json`, logs and the case directories `n2/` through
`n6/`. The ZIP contains exactly the versioned repository files at the
recorded commit. A complete run requires zero unresolved boxes and verified
terminal conclusions for every graph, candidate and equality object.

Select cases with `-DiskCounts 5,6`. For one case without release packaging:

```powershell
julia --startup-file=no --threads=1,0 --gcthreads=1 --project=. --load=scripts/single_threaded.jl scripts/run_proof.jl --n 6 --precision 512
```

This writes `output/n6/`. Use `--n 2`, `3`, `4` or `5` for the other cases;
`--output output/check` changes the case directory to `output/check/nN/`.
Search limits produce an unresolved outcome and prevent proof closure.

## Tests

```powershell
julia --project=. test/runtests.jl
& .\test\native_logging.ps1
```

## Contents

- `src/`: exact arithmetic, a common interval search, refinement and analysis.
- `data/`: graph censuses and algebraic configurations and families.
- `scripts/`: proof launchers and census generation.
- `test/`: mathematical and integration checks.
- `output/`: generated files, ignored by Git.

[Mathematical verification](docs/proof.md) ·
[Census](docs/census.md) · [Evidence formats](docs/evidence.md)

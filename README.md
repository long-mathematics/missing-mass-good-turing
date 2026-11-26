# Missing-mass extremizers and canonical Good–Turing risk

Two mathematical manuscripts by Christopher D. Long, revised November 26, 2026:

- [Four Missing-Mass Extremizer Problems](paper/main.pdf) ([TeX](paper/main.tex)): the exact one-transition/two-transition diagram for the Poissonized positive missing-mass MGF; a fixed-sample active-level bound; and exact tail saturation with a bounded-support maximizer.
- [Canonical Poissonized Good–Turing Risk: One-Exceptional Extremizers](good-turing/main.pdf) ([TeX](good-turing/main.tex)): every worst-case law for the canonical unbiased order-r Poissonized estimator is uniform on its positive support or has one exceptional positive atom. For ordinary missing mass, support is full and the exceptional atom is heavy.

The November audit recovered the original July manuscript, checked its proof dependencies, corrected four supporting statements, and added independent symbolic, exact rational, and numerical checks. See [AUDIT.md](AUDIT.md) for the claim ledger and limitations, and [PROVENANCE.md](PROVENANCE.md) for source history and hashes.

The Good–Turing companion is a structural theorem for a fixed estimator. It does **not** establish a complete risk phase diagram, a general fixed-count extremizer theorem, or an exact minimax estimator. Historical manuscripts in [archive/](archive/README.md) contain claims outside the audited scope and are retained only as research provenance.

## Reproduce

Python 3.12 or newer and a LaTeX distribution with `pdflatex`, `amsmath`, `microtype`, `hyperref`, and `cleveref` are sufficient.

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
make check
make pdf
```

Alternatively, `uv venv .venv` and `uv pip install --python .venv/bin/python -r requirements.txt` install the Python dependencies.

`make check` writes JSON results to `checks/`. Finite searches are stress tests, not proofs of global optimality or absence of counterexamples. The mathematical arguments are in the manuscripts. `make pdf` compiles only the two current manuscripts, never the historical archive.

## Current scope

| Result | Status |
|---|---|
| Poissonized positive-MGF one-heavy reduction, duality, exact zero sets | Complete proof, audited and corrected |
| Continuity and strict increase of the late boundary | Complete proof, including first-return exclusion of pinch zeros |
| Fixed-sample MGF at most t(t+1) positive levels | Complete proof; not a sharp classification |
| Fixed-tail saturation and existence of bounded-support maximizer | Complete proof; nonsaturated optimizer classification open |
| Canonical order-r Poissonized Good–Turing one-exceptional structure | Complete proof in separate companion |
| Ordinary Good–Turing complete phase expansion from older draft | Not certified by this audit |
| General fixed-count Good–Turing and exact minimax claims | Not established here |

A targeted literature check verified the principal antecedents, but this repository does not claim an exhaustive priority search or external peer review. No journal or preprint submission is implied by repository publication.

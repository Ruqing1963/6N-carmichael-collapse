# 6N-carmichael-collapse

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20610382.svg)](https://doi.org/10.5281/zenodo.20610382)

**Part XXXVI** of *Arithmetic Geodynamics on the 6N Skeleton* — Ruqing Chen, 2026.

> **Topological Collapse and Asymmetric Divergence of Carmichael Pseudoprimes on the 6N Skeleton**

This repository contains the full reproducer, data, figure, and paper for the empirical study of
Carmichael numbers read through the 6N skeleton (the wings `6N+1` / `6N-1`).

## What the paper shows

On the 6N skeleton every integer coprime to 6 sits on a **wing**: the *right* wing `6N+1`
(≡ 1 mod 6) or the *left* wing `6N-1` (≡ 5 mod 6). Korselt's criterion (`p-1 | n-1` for every
prime `p | n`) forces three elementary but rigid structural laws, which we state and machine-verify:

1. **Semipermeable membrane (Prop. 1).** Every Carmichael number on the **left** wing has *all* of
   its prime factors on the left wing, and is not divisible by 3.
2. **Parity law (Prop. 2).** A Carmichael number coprime to 6 lies on the **left** wing iff it has an
   **odd** number of left-wing (`6N-1`) prime factors; on the **right** wing iff that number is even.
3. **Unified barrier (Prop. 3).** Whenever `3 ∤ (n-1)`, no `6N+1` prime can divide `n`. This fuses the
   left wing with the off-skeleton residue-3 numbers: a residue-3 Carmichael number is exactly a 3
   glued to a pure-left core (e.g. `561 = 3 × 11 × 17`). The single gate `3 | (n-1)` governs all three
   classes.

Enumerating **all 1547 Carmichael numbers up to 10¹⁰**, all three laws hold with **zero violations**.

The geometric consequence is a violent population asymmetry. The right-to-left ratio
`R(X) = #{right} / #{left}` rises monotonically across seven decades, and the same monochromatic
suppression acts *inside* the right wing, where the pure-right family is overtaken by the mixed family
(`M(X) = #{pure-right} / #{mixed}` crosses 1 between 10⁷ and 10⁸):

| log₁₀X | C(X) | right | left | R(X) | pure-right | mixed | M(X) | mean #factors |
|-------:|-----:|------:|-----:|-----:|-----------:|------:|-----:|--------------:|
| 4  | 7    | 5    | 1  | 5.00  | 3   | 2   | 1.50 | 3.00 |
| 5  | 16   | 13   | 1  | 13.00 | 8   | 5   | 1.60 | 3.25 |
| 6  | 43   | 38   | 2  | 19.00 | 23  | 15  | 1.53 | 3.49 |
| 7  | 105  | 98   | 4  | 24.50 | 53  | 45  | 1.18 | 3.58 |
| 8  | 255  | 243  | 6  | 40.50 | 114 | 129 | 0.88 | 3.78 |
| 9  | 646  | 623  | 12 | 51.92 | 254 | 369 | 0.69 | 4.00 |
| 10 | 1547 | 1510 | 17 | 88.82 | 587 | 923 | 0.64 | 4.23 |

**Mechanism.** A single-colour (monochromatic) core requires all factors in one wing, while the mean
number of prime factors grows with X, so single-colour states become exponentially scarce — the left
wing collapses outright and pure-right is overtaken by mixed.

**Honest ledger.** The two single-colour basins decay in the same *direction* but at materially
different *rates* (`M ~ (log X)^-1.1` vs `R ~ (log X)^+2.9`, a ~2.6× gap traced to the smoothness of
`p-1` for right primes vs its roughness for left primes), so **no single thermodynamic law** governs
them; and neither asymptotic rate is identifiable (Poisson sparsity + deep pre-asymptotic drift). What
is unified is **algebraic** — the single barrier `3 ∤ (n-1)` — not thermodynamic.

## Reproduce

```bash
python3 carmichael_6N_collapse.py 1e10     # ~20 s, single core; default bound is 1e10
python3 carmichael_6N_collapse.py 1e9      # ~3 s
```

Dependencies: `numpy`, `matplotlib`.

The script enumerates every Carmichael number ≤ X by the Korselt divisor method (recovering each
factorisation), **asserts** the count against the known value, runs the three structural-law verifiers,
prints the wing census, and writes the data files and the figure. The three `assert` conditions in
`verify_laws` **are** the machine-verification of Propositions 1, 2, and 3: if they pass, the membrane,
parity, and unified-barrier laws hold with zero violations on the entire enumerated set.

## Files

| File | Description |
|------|-------------|
| `carmichael_6N_collapse.py` | Self-contained reproducer (enumeration + 3-law verifiers + M(X) + CSV/figure output) |
| `Chen_6N_Paper36_v2.tex` | Paper source (v2: adds Prop. 3 and the M(X) section) |
| `Chen_6N_Paper36_v2.pdf` | Compiled paper (6 pp.) |
| `carmichael_Rx.pdf` / `.png` | Figure 1 (asymmetric divergence + mechanism) |
| `data/carmichael_le1e10.csv` | All 1547 Carmichael numbers ≤ 10¹⁰: `n, factorization, num_factors, n_mod6, wing, num_left_factors` |
| `data/wing_census.csv` | Cumulative wing census, `R(X)`, pure-right/mixed/`M(X)` by decade |

## Scope

The three structural laws are **elementary consequences of Korselt's criterion**, stated here in wing
coordinates and verified empirically; we do not claim them as new theorems. The framework says nothing
about the infinitude or the true counting function of Carmichael numbers beyond Alford–Granville–
Pomerance (1994), Erdős (1956), and Pinch.

## Citation

R. Chen, *Topological Collapse and Asymmetric Divergence of Carmichael Pseudoprimes on the 6N
Skeleton* (Part XXXVI), Zenodo, 2026. [doi:10.5281/zenodo.20610382](https://doi.org/10.5281/zenodo.20610382).
See also the series review (Parts I–XIX), [doi:10.5281/zenodo.20585301](https://doi.org/10.5281/zenodo.20585301).

## License

MIT — see [LICENSE](LICENSE).

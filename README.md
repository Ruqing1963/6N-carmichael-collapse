# 6N-carmichael-collapse

**Part XXXVI** of *Arithmetic Geodynamics on the 6N Skeleton* — Ruqing Chen, 2026.

> **Topological Collapse and Asymmetric Divergence of Carmichael Pseudoprimes on the 6N Skeleton**

This repository contains the full reproducer, data, figure, and paper for the empirical study of
Carmichael numbers read through the 6N skeleton (the wings `6N+1` / `6N-1`).

## What the paper shows

On the 6N skeleton every integer coprime to 6 sits on a **wing**: the *right* wing `6N+1`
(≡ 1 mod 6) or the *left* wing `6N-1` (≡ 5 mod 6). Korselt's criterion (`p-1 | n-1` for every
prime `p | n`) forces two elementary but rigid structural laws, which we state and machine-verify:

1. **Semipermeable membrane (Prop. 1).** Every Carmichael number on the **left** wing has *all* of
   its prime factors on the left wing, and is not divisible by 3. No `6N+1` factor, and no 3, can
   enter a left-wing Carmichael number.
2. **Parity law (Prop. 2).** A Carmichael number coprime to 6 lies on the **left** wing iff it has an
   **odd** number of left-wing (`6N-1`) prime factors; on the **right** wing iff that number is even.

Enumerating **all 1547 Carmichael numbers up to 10¹⁰**, both laws hold with **zero violations**.

The geometric consequence is a violent population asymmetry. The right-to-left ratio
`R(X) = #{right} / #{left}` rises monotonically across seven decades:

| log₁₀X | C(X) | right | left | R(X) | mean #factors |
|-------:|-----:|------:|-----:|-----:|--------------:|
| 4  | 7    | 5    | 1  | 5.00  | 3.00 |
| 5  | 16   | 13   | 1  | 13.00 | 3.25 |
| 6  | 43   | 38   | 2  | 19.00 | 3.49 |
| 7  | 105  | 98   | 4  | 24.50 | 3.58 |
| 8  | 255  | 243  | 6  | 40.50 | 3.78 |
| 9  | 646  | 623  | 12 | 51.92 | 4.00 |
| 10 | 1547 | 1510 | 17 | 88.82 | 4.23 |

**Mechanism.** The left wing admits only *monochromatic* products (all factors ≡ 5 mod 6) of *odd*
length, while the mean number of prime factors of Carmichael numbers grows with X, so the
monochromatic configuration becomes exponentially scarce.

**Honest ledger.** The *divergence* of `R(X)` is firm and structurally explained, but its **rate is
not identifiable**: a power of `log X`, a small power of `X`, and a stretched exponential all fit the
seven points equally well, and the left-wing counts (1, 1, 2, 4, 6, 12, 17) carry large Poisson
noise. No reachable computational window resolves the asymptotic law. We claim the structure and the
phenomenon, not an asymptotic divergence law.

## Reproduce

```bash
python3 carmichael_6N_collapse.py 1e10     # ~20 s, single core; default bound is 1e10
python3 carmichael_6N_collapse.py 1e9      # ~3 s
```

Dependencies: `numpy`, `matplotlib`.

The script enumerates every Carmichael number ≤ X by the Korselt divisor method (recovering each
factorisation), **asserts** the count against the known value, runs the two structural-law verifiers,
prints the wing census, and writes the data files and the figure. The two `assert` statements in
`verify_laws` **are** the machine-verification of Propositions 1 and 2: if they pass, the membrane and
parity laws hold with zero violations on the entire enumerated set.

## Files

| File | Description |
|------|-------------|
| `carmichael_6N_collapse.py` | Self-contained reproducer (enumeration + verifiers + CSV/figure output) |
| `Chen_6N_Paper36.tex` | Paper source |
| `Chen_6N_Paper36.pdf` | Compiled paper (5 pp.) |
| `carmichael_Rx.pdf` / `.png` | Figure 1 (asymmetric divergence + mechanism) |
| `data/carmichael_le1e10.csv` | All 1547 Carmichael numbers ≤ 10¹⁰: `n, factorization, num_factors, n_mod6, wing, num_left_factors` |
| `data/wing_census.csv` | Cumulative wing census and `R(X)` by decade |

## Scope

The two structural laws are **elementary consequences of Korselt's criterion**, stated here in wing
coordinates and verified empirically; we do not claim them as new theorems. The framework says nothing
about the infinitude or the true counting function of Carmichael numbers beyond Alford–Granville–
Pomerance (1994), Erdős (1956), and Pinch.

## Citation

R. Chen, *Topological Collapse and Asymmetric Divergence of Carmichael Pseudoprimes on the 6N
Skeleton* (Part XXXVI), 2026. See the series review (Parts I–XIX), Zenodo
[doi:10.5281/zenodo.20585301](https://doi.org/10.5281/zenodo.20585301).

## License

MIT — see [LICENSE](LICENSE).

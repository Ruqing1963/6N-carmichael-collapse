#!/usr/bin/env python3
"""
carmichael_6N_collapse.py  --  reproducer for
"Topological Collapse and Asymmetric Divergence of Carmichael Pseudoprimes
 on the 6N Skeleton" (Part XXXVI), R. Chen, 2026.

What it does
------------
  1. Enumerates ALL Carmichael numbers up to X (default 1e10) by the Korselt
     divisor method, recovering each number's full prime factorisation.
  2. Verifies the count against the known value C(1e10)=1547 (Pinch).
  3. Tabulates the cumulative wing census and the right/left ratio R(X) by decade.
  4. Runs the two structural-law verifiers with a hard zero-violation assertion:
        - membrane: every left-wing (n=5 mod6) Carmichael has ALL factors =5 mod6,
                    and is not divisible by 3;
        - parity  : every coprime-to-6 Carmichael lies on the left wing iff it has
                    an odd number of factors =5 mod6.
  5. Reports the mean prime-factor count kbar(X) (the divergence driver) and writes
     the R(X) figure.

Dependencies: numpy, matplotlib.  Single-core; ~20 s to 1e10, ~3 s to 1e9.

Usage:   python3 carmichael_6N_collapse.py            # X = 1e10
         python3 carmichael_6N_collapse.py 1e9         # custom bound
"""
import sys, time
import numpy as np

# ----------------------------------------------------------------------
# smallest-prime-factor sieve (factors n-1 targets up to SPF_N in O(log n))
# ----------------------------------------------------------------------
SPF_N = 2 * 10**7
def build_spf(n):
    spf = np.zeros(n + 1, dtype=np.int32)
    for i in range(2, int(n**0.5) + 1):
        if spf[i] == 0:                 # i is prime
            seg = spf[i * i::i]
            seg[seg == 0] = i
    return spf

# small primes used as inner factor candidates and as trial divisors
def build_small_primes(lim):
    s = np.ones(lim + 1, dtype=bool); s[:2] = False
    for i in range(2, int(lim**0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    p = np.nonzero(s)[0]
    return p[p >= 3]                     # odd primes from 3 (Carmichaels are odd)

# ----------------------------------------------------------------------
# deterministic Miller-Rabin (correct for all n < 3.3e24)
# ----------------------------------------------------------------------
_MR_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
def is_prime(n):
    if n < 2:
        return False
    for p in _MR_BASES:
        if n % p == 0:
            return n == p
    d = n - 1; r = 0
    while d % 2 == 0:
        d //= 2; r += 1
    for a in _MR_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(r - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

# ----------------------------------------------------------------------
# Carmichael enumeration by the Korselt divisor method
#
# A squarefree n = p_1 ... p_k (k>=3) is Carmichael iff (p_i - 1) | (n - 1) for all i.
# Build the ordered prefix of inner primes; the final (largest) prime r then satisfies
# n - 1 = M r - 1 == M - 1  (mod r-1), so (r-1) | (M-1).  Inner factors are checked
# explicitly.  Small M: factor M-1 via the SPF sieve and enumerate divisors.  Large M:
# only divisors d <= X//M matter (the final prime must fit under X), found by a short
# bounded scan.
# ----------------------------------------------------------------------
def divisors_from_spf(spf, t):
    f = {}
    while t > 1:
        p = int(spf[t]); p = p if p else t
        while t % p == 0:
            t //= p; f[p] = f.get(p, 0) + 1
    divs = [1]
    for p, c in f.items():
        divs = [d * p**e for d in divs for e in range(c + 1)]
    return divs

def enumerate_carmichael(X, spf, small_primes):
    SP = small_primes.tolist()
    out = []
    def rec(primes, M, start):
        last = primes[-1]
        if len(primes) >= 2:                 # try to close with a final prime r > last
            cap = X // M                      # r <= cap
            if cap > last:
                if M - 1 <= SPF_N:
                    for d in divisors_from_spf(spf, M - 1):
                        r = d + 1
                        if last < r <= cap and is_prime(r):
                            n = M * r
                            if all((n - 1) % (p - 1) == 0 for p in primes):
                                out.append((n, tuple(primes) + (r,)))
                else:                          # M large -> only small divisors fit under X
                    Mm = M - 1
                    for d in range(last, cap):
                        if Mm % d == 0:
                            r = d + 1
                            if is_prime(r):
                                n = M * r
                                if all((n - 1) % (p - 1) == 0 for p in primes):
                                    out.append((n, tuple(primes) + (r,)))
        for idx in range(start, len(SP)):      # extend with another inner prime
            p = SP[idx]
            if p <= last:
                continue
            newM = M * p
            if newM * p > X:                   # need a final prime > p as well
                break
            rec(primes + [p], newM, idx + 1)
    for i, p in enumerate(SP):
        if p * p * p > X:
            break
        rec([p], p, i + 1)
    return sorted(set(out))

# ----------------------------------------------------------------------
# structural-law verifiers (Propositions 1 and 2): must return 0 violations
# ----------------------------------------------------------------------
def verify_laws(carms):
    membrane_violations = 0   # Prop 1: left wing -> all factors =5 mod6
    parity_violations = 0     # Prop 2: right wing -> #left factors even
    barrier_violations = 0    # Prop 3: 3|n (residue 3) -> all non-3 factors =5 mod6
    for n, F in carms:
        res = n % 6
        if res == 5:                                   # left wing
            if 3 in F or any(p % 6 != 5 for p in F):   # membrane: all factors must be 5 mod6
                membrane_violations += 1
        elif res == 1:                                 # right wing
            n_left = sum(1 for p in F if p % 6 == 5)   # parity: #left must be even
            if n_left % 2 != 0:
                parity_violations += 1
        elif res == 3:                                 # ghost (3 | n)
            if any(p % 6 == 1 for p in F):             # unified barrier: no right-wing factor
                barrier_violations += 1
    return membrane_violations, parity_violations, barrier_violations

# ----------------------------------------------------------------------
def main():
    X = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10**10
    t0 = time.time()
    spf = build_spf(SPF_N)
    # inner-prime candidates must reach ~ sqrt(X/3); trial divisors must reach ~ sqrt(X^{2/3})
    plim = max(1_200_000, int((X / 3) ** 0.5) + 1)
    small_primes = build_small_primes(plim)
    carms = enumerate_carmichael(X, spf, small_primes)
    elapsed = time.time() - t0

    ns = np.array([n for n, _ in carms])
    res = ns % 6
    ks = np.array([len(F) for _, F in carms])
    nl = np.array([sum(1 for p in F if p % 6 == 5) for _, F in carms])  # #left factors

    print(f"Enumerated {len(carms)} Carmichael numbers up to {X:.0e} in {elapsed:.1f}s")
    known = {10**9: 646, 10**10: 1547, 10**11: 3522, 10**12: 8241}
    if X in known:
        assert len(carms) == known[X], f"count {len(carms)} != known {known[X]}"
        print(f"  count matches known C({X:.0e}) = {known[X]}  [OK]")

    mv, pv, bv = verify_laws(carms)
    print(f"\nStructural-law verification (must be zero):")
    print(f"  Prop 1 membrane violations (left-wing impurity)     : {mv}")
    print(f"  Prop 2 parity   violations (right-wing odd #left)    : {pv}")
    print(f"  Prop 3 barrier  violations (residue-3 with =1 factor): {bv}")
    assert mv == 0 and pv == 0 and bv == 0, "STRUCTURAL LAW VIOLATED"
    print("  -> all three laws hold with zero violations  [OK]")

    print(f"\n{'log10X':>6} {'C(X)':>6} {'right':>6} {'left':>5} {'div3':>5} "
          f"{'R(X)':>8} {'pureR':>6} {'mixed':>6} {'M(X)':>7} {'kbar':>6}")
    decades = []
    for e in range(4, int(np.log10(X)) + 1):
        Xd = 10**e; m = ns <= Xd
        R = int(((res == 1) & m).sum())
        L = int(((res == 5) & m).sum())
        D3 = int(((res == 3) & m).sum())
        pureR = int(((res == 1) & (nl == 0) & m).sum())
        mixed = int(((res == 1) & (nl > 0) & m).sum())
        kbar = ks[m].mean() if m.any() else float('nan')
        rr = R / L if L else float('nan')
        MM = pureR / mixed if mixed else float('nan')
        decades.append((e, int(m.sum()), R, L, D3, rr, pureR, mixed, MM, kbar))
        print(f"{e:6d} {int(m.sum()):6d} {R:6d} {L:5d} {D3:5d} {rr:8.2f} "
              f"{pureR:6d} {mixed:6d} {MM:7.3f} {kbar:6.2f}")

    # ---- write CSV data files -------------------------------------------------
    import csv
    wing = {1: "right", 5: "left", 3: "div3"}
    fn_full = f"carmichael_le{X:.0e}.csv".replace("e+0", "e").replace("e+", "e")
    with open(fn_full, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["n", "factorization", "num_factors", "n_mod6", "wing", "num_left_factors"])
        for n, F in carms:
            w.writerow([n, "*".join(map(str, F)), len(F), n % 6,
                        wing.get(n % 6, "?"), sum(1 for p in F if p % 6 == 5)])
    with open("wing_census.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["log10_X", "C_X", "right", "left", "div3", "R_X_right_over_left",
                    "pure_right", "mixed", "M_X_pureright_over_mixed", "mean_num_factors"])
        for e, c_x, R, L, D3, rr, pureR, mixed, MM, kbar in decades:
            w.writerow([e, c_x, R, L, D3, f"{rr:.4f}", pureR, mixed, f"{MM:.4f}", f"{kbar:.4f}"])
    print(f"\nWrote {fn_full} (all {len(carms)} numbers) and wing_census.csv")

    # figure (optional; skipped if matplotlib unavailable)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        D = np.array(decades, dtype=float)
        E, Rv, kbar = D[:, 0], D[:, 5], D[:, 9]
        good = np.isfinite(Rv)
        lnX = E[good] * np.log(10)
        A = np.polyfit(np.log(lnX), np.log(Rv[good]), 1)
        fig, ax = plt.subplots(1, 2, figsize=(13, 5))
        ax[0].plot(E[good], Rv[good], 'o-', color="#c0392b", lw=2, ms=7,
                   label=r"measured $R(X)=\#$right$/\#$left")
        xx = np.linspace(E[good].min(), E[good].max() + 1, 100); lnxx = xx * np.log(10)
        ax[0].plot(xx, np.exp(np.polyval(A, np.log(lnxx))), '--', color="#1f5fbf", lw=1.3,
                   label=fr"$\sim(\log X)^{{{A[0]:.1f}}}$ (reference; rate not pinnable)")
        ax[0].set_xlabel(r"$\log_{10} X$"); ax[0].set_ylabel(r"$R(X)$")
        ax[0].set_title("Asymmetric divergence"); ax[0].legend(fontsize=9); ax[0].grid(alpha=.3)
        ax[1].plot(E[good], Rv[good], 'o-', color="#c0392b", lw=2, ms=6, label=r"$R(X)$")
        ax[1].plot(E[good], 2 ** (kbar[good] - 1), 's--', color="#2e7d32", lw=1.5, ms=6,
                   label=r"$2^{\bar k(X)-1}$ (direction only)")
        ax[1].set_xlabel(r"$\log_{10} X$")
        ax[1].set_title(r"Mechanism: $R\sim2^{\bar k}$, $\bar k$ grows")
        ax[1].legend(fontsize=9); ax[1].grid(alpha=.3)
        fig.suptitle(f"Carmichael wing asymmetry on the 6N skeleton (to {X:.0e})")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig("carmichael_Rx.pdf")
        fig.savefig("carmichael_Rx.png", dpi=150)
        print("Wrote carmichael_Rx.pdf and carmichael_Rx.png")
    except ImportError:
        print("\n(matplotlib not available; skipped figure)")

if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Verifier r05 (independent; does NOT import donutloc).

Noise-free population bias of the naive MLE under TCP misalignment.
Model (as stated in paper Sec. VIII B): LG donut Eq. S17 (fwhm=300 nm, peak 1), TCP of diameter
L=100 nm with peripheral zeros at angles pi/2 + 2 pi k/3 and a central zero; every one of the four
zeros displaced by delta in an independent uniform random direction; Eq. S30 background with
SBR=10 in both the true and the nominal model; naive MLE = argmax sum_i q_i ln p_nom_i(r) with
q = p_true(r0) (expected counts; N cancels).

Method: vectorized Fisher scoring started at r0 (converges to the local maximum near r0), plus a
brute-force grid check over the disk of radius 0.75 L on a subset of patterns.  Also the
linear-response (delta -> 0) limit of |bias|/delta from the exact first-order map.
"""
import sys
import numpy as np

FWHM = 300.0
A = 4 * np.log(2) / FWHM ** 2
L = 100.0
SBR = 10.0
S = SBR / (SBR + 1)


def ideal(Lv=L):
    ang = np.pi / 2 + 2 * np.pi * np.arange(3) / 3
    b = np.zeros((4, 2))
    b[:3, 0] = Lv / 2 * np.cos(ang)
    b[:3, 1] = Lv / 2 * np.sin(ang)
    return b


def lg(u):  # u = rho^2
    return np.e * A * u * np.exp(-A * u)


def dlg(u):  # d lg / du
    return np.e * A * np.exp(-A * u) * (1 - A * u)


def probs(r, b, grad=False):
    """r (..., 2), b (..., 4, 2) broadcastable -> p (..., 4) [and dp (..., 4, 2)]."""
    d = r[..., None, :] - b                    # (...,4,2)
    u = np.sum(d * d, axis=-1)
    I = lg(u)
    Ssum = I.sum(-1, keepdims=True)
    p0 = I / Ssum
    p = S * p0 + (1 - S) / 4
    if not grad:
        return p
    dI = (dlg(u)[..., None] * 2 * d)           # (...,4,2)
    dS = dI.sum(-2, keepdims=True)
    dp0 = dI / Ssum[..., None] - I[..., None] * dS / Ssum[..., None] ** 2
    return p, S * dp0


def naive_mle_scoring(q, r0, bnom, iters=60):
    """Fisher scoring for argmax sum q ln p_nom(r); q (P,4), r0 (P,2)."""
    r = r0.copy()
    for _ in range(iters):
        p, dp = probs(r, bnom, grad=True)
        score = np.einsum("pi,pik->pk", q / p, dp)
        F = np.einsum("pik,pil->pkl", dp / p[..., None], dp)
        step = np.linalg.solve(F, score[..., None])[..., 0]
        r = r + step
        if np.max(np.abs(step)) < 1e-10:
            break
    return r


def draw(P, delta, rng):
    phi = rng.uniform(0, 2 * np.pi, size=(P, 4))
    disp = delta * np.stack([np.cos(phi), np.sin(phi)], -1)
    return ideal()[None] + disp


def pop_bias(P, delta, pos, rng_seed, grid_check=0):
    rng = np.random.default_rng(rng_seed)
    btrue = draw(P, delta, rng)
    r0 = np.tile(np.asarray(pos, float), (P, 1))
    q = probs(r0, btrue)
    bnom = ideal()[None]
    rh = naive_mle_scoring(q, r0, bnom)
    bias = np.hypot(*(rh - r0).T)
    bad = 0
    if grid_check:
        # brute force over the disk of radius 0.75 L (step 0.5 nm) for the first grid_check patterns
        g = np.arange(-75, 75.01, 0.5)
        X, Y = np.meshgrid(g, g)
        m = X ** 2 + Y ** 2 <= 75 ** 2
        G = np.stack([X[m], Y[m]], -1)
        pg = probs(G, ideal()[None])            # (G,4)
        lpg = np.log(pg)
        for j in range(grid_check):
            ll = lpg @ q[j]
            k = np.argmax(ll)
            if np.hypot(*(G[k] - rh[j])) > 1.0:
                bad += 1
    return bias, bad


def linear_ratio(pos, M=400000, seed=7):
    """delta->0 limit of E|bias|/delta via the exact first-order map (finite differences in the
    zero positions, h=1e-4 nm)."""
    b0 = ideal()
    r0 = np.asarray(pos, float)
    p_nom, dp = probs(r0[None], b0[None], grad=True)
    p_nom, dp = p_nom[0], dp[0]
    F = np.einsum("ik,il->kl", dp / p_nom[:, None], dp)
    # first-order change of q with each zero coordinate
    h = 1e-4
    D = np.zeros((4, 8))
    for c in range(8):
        bb = b0.copy().reshape(-1)
        bb[c] += h
        pp = probs(r0[None], bb.reshape(4, 2)[None])[0]
        bb[c] -= 2 * h
        pm = probs(r0[None], bb.reshape(4, 2)[None])[0]
        D[:, c] = (pp - pm) / (2 * h)
    # KL-minimizer: F dr = sum_i dq_i dp_i / p_i
    Mmap = np.linalg.solve(F, (dp / p_nom[:, None]).T @ D)   # (2,8)
    rng = np.random.default_rng(seed)
    phi = rng.uniform(0, 2 * np.pi, size=(M, 4))
    d = np.stack([np.cos(phi), np.sin(phi)], -1).reshape(M, 8)
    b = d @ Mmap.T
    nb = np.hypot(*b.T)
    return nb.mean(), nb.std(ddof=1) / np.sqrt(M)


if __name__ == "__main__":
    out = {}
    positions = {"center": (0.0, 0.0), "Lq": (L / 4, 0.0)}
    deltas = [2.0, 5.0, 10.0]
    for tag, pos in positions.items():
        lr, lse = linear_ratio(pos)
        print("%s linear-response limit |bias|/delta = %.4f +- %.4f" % (tag, lr, lse))
    P = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    for tag, pos in positions.items():
        for seed in (101, 202, 303):
            ratios = []
            per = []
            for d in deltas:
                bias, bad = pop_bias(P, d, pos, seed, grid_check=200 if seed == 101 else 0)
                ratios.append((bias.mean() / d, bias.std(ddof=1) / np.sqrt(P) / d, bad))
                per.append(bias)
            dl = np.array(deltas)
            perp = (dl[:, None] * np.array(per)).sum(0) / np.sum(dl ** 2)
            print("%s seed %d P=%d: ratios %s ; LSQ slope %.4f +- %.4f ; grid mismatches %s"
                  % (tag, seed, P, " ".join("%.4f(%.4f)" % r[:2] for r in ratios),
                     perp.mean(), perp.std(ddof=1) / np.sqrt(P), [r[2] for r in ratios]))
            # running means for convergence
            if seed == 101:
                for n in (400, 1000, 4000, P):
                    print("   first %6d patterns: slope %.4f" % (n, perp[:n].mean()))
                # tail diagnostics
                print("   per-pattern slope: median %.4f, p99 %.4f, max %.4f, skew-ish mean-median %.4f"
                      % (np.median(perp), np.percentile(perp, 99), perp.max(), perp.mean() - np.median(perp)))

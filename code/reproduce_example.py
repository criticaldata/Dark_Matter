"""Reproduce every number in the worked example (main Section 4.4, SI S5).

The Dark Matter of Clinical Data. Requires Python 3 with numpy and scipy.
Run:  python reproduce_example.py
"""
from fractions import Fraction as F
from math import sqrt

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

# Model M1 (Assumption 4 holds): case mix, symmetric channel, shared mechanism
q = {0: F(3, 10), 1: F(7, 10)}          # P(Z*=1 | E=e)
a1, b1 = F(9, 10), F(1, 10)             # P(Z=1 | Z*=1), P(Z=1 | Z*=0)
mu0, mu1 = F(1, 5), F(4, 5)             # P(Y=1 | Z*=0), P(Y=1 | Z*=1)

rate = {}
print("Exact quantities (Steps 4-6, Table S3, joint tables):")
for e in (0, 1):
    pz1 = a1 * q[e] + b1 * (1 - q[e])
    pi1 = a1 * q[e] / pz1                    # P(Z*=1 | Z=1, e)
    pi0 = (1 - a1) * q[e] / (1 - pz1)        # P(Z*=1 | Z=0, e)
    y1 = mu0 + (mu1 - mu0) * pi1             # P(Y=1 | Z=1, e)
    y0 = mu0 + (mu1 - mu0) * pi0             # P(Y=1 | Z=0, e)
    rate[(1, e)], rate[(0, e)] = y1, y0
    joint = [y1 * pz1, (1 - y1) * pz1, y0 * (1 - pz1), (1 - y0) * (1 - pz1)]
    print(f"  site {e}: P(Z=1)={pz1}, posteriors {pi1}, {pi0}; "
          f"P(Y=1|Z=1)={y1}, P(Y=1|Z=0)={y0}; P(Y,Z) at (1,1),(0,1),(1,0),(0,0) = "
          + ", ".join(str(x) for x in joint))

gap1, gap0 = rate[(1, 1)] - rate[(1, 0)], rate[(0, 1)] - rate[(0, 0)]
print(f"Gap in Z=1 cell: {gap1} = {float(gap1):.4f}; Z=0 cell: {gap0}")

# Power: two-sided two-proportion test, pooled variance under the null
p1, p2 = float(rate[(1, 0)]), float(rate[(1, 1)])
pb = (p1 + p2) / 2
za, zb = norm.ppf(0.975), norm.ppf(0.80)
n = (za * sqrt(2 * pb * (1 - pb)) + zb * sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2 / (p1 - p2) ** 2
print(f"Patients per site in Z=1 stratum for 80% power: {n:.1f}")

# Sensitivity analysis: symmetric site-invariant error rate eps
pz = {0: 0.34, 1: 0.66}
obs = {k: float(v) for k, v in rate.items()}


def posteriors(eps):
    out = {}
    for e in (0, 1):
        qe = (pz[e] - eps) / (1 - 2 * eps)
        out[(1, e)] = (1 - eps) * qe / pz[e]
        out[(0, e)] = eps * qe / (1 - pz[e])
    return out


def one_cell_gap(eps):
    p = posteriors(eps)
    return p[(1, 1)] - p[(1, 0)]


print("One-cell bound |pi_1(1) - pi_0(1)| (Table S5):")
for eps in (0.01, 0.02, 0.05, 0.063, 0.10, 0.15):
    print(f"  eps={eps}: {one_cell_gap(eps):.4f}")
thr = brentq(lambda x: one_cell_gap(x) - float(gap1), 1e-3, 0.2)
print(f"Threshold eps at which the one-cell bound equals the gap: {thr:.4f}")

print("Two-cell least-squares fit of (mu0, mu1) under Assumption 4:")
for eps in (0.05, 0.063, 0.08, 0.10, 0.12):
    p = posteriors(eps)
    keys = [(1, 0), (0, 0), (1, 1), (0, 1)]
    A = np.array([[1 - p[k], p[k]] for k in keys])
    y = np.array([obs[k] for k in keys])
    mu, *_ = np.linalg.lstsq(A, y, rcond=None)
    print(f"  eps={eps}: mu0={mu[0]:.3f}, mu1={mu[1]:.3f}, "
          f"max |residual|={np.abs(A @ mu - y).max():.4f}")


# ---------------------------------------------------------------------------
# Theorem 1(i): two completions of the untested analyte law (SI S3 instance)
# ---------------------------------------------------------------------------
def H(p):
    p = np.asarray(p, float).ravel(); p = p[p > 0]
    return float(-(p * np.log2(p)).sum())


def MI(pxy):
    pxy = np.asarray(pxy, float)
    return H(pxy.sum(1)) + H(pxy.sum(0)) - H(pxy)


def I_min(P):
    """Williams-Beer redundancy of sources (axis 1, axis 2) about Y (axis 0)."""
    py, tot = P.sum((1, 2)), 0.0
    for y in range(P.shape[0]):
        spec = []
        for ax in (2, 1):
            pya = P.sum(ax); pa = pya.sum(0); m = pya[y] > 0
            spec.append(float((pya[y][m] / py[y] * np.log2(pya[y][m] / pa[m] / py[y])).sum()))
        tot += py[y] * min(spec)
    return tot


obs_s1 = np.array([[0.30, 0.10], [0.05, 0.55]])   # P(Y=y, V=v | S=1)
pY_s0 = np.array([0.7, 0.3])                      # P(Y=y | S=0)
print("Theorem 1(i): completions agreeing on all observed quantities")
for name, q in (("A (V independent of Y)", [[0.5, 0.5], [0.5, 0.5]]),
                ("B (V equals Y with prob. 0.95)", [[0.95, 0.05], [0.05, 0.95]])):
    P = np.zeros((2, 2, 2))                       # P[y, v, s]
    P[:, :, 1] = 0.5 * obs_s1
    for y in range(2):
        P[y, :, 0] = 0.5 * pY_s0[y] * np.array(q[y])
    total, iyv, iys = MI(P.reshape(2, 4)), MI(P.sum(2)), MI(P.sum(1))
    red = I_min(P)
    print(f"  {name}: I(Y;V,S)={total:.3f}  I(Y;V)={iyv:.3f}  I(Y;S)={iys:.3f}  "
          f"Red={red:.3f}  Unq_V={iyv - red:.3f}  Syn={total - iyv - iys + red:.3f}")

"""Numerical check of Proposition 5 (two tests and a case-mix difference).

For random models with binary physiology, two sites with different case mix,
two conditionally independent binary tests with site-invariant error rates,
and a site-specific outcome mechanism, fit the model to its own observed
distribution from many random starts and confirm that every exact fit
recovers the true parameters (after fixing labels by sens + spec > 1).
Requires numpy and scipy.  Run:  python check_two_tests.py
"""
import numpy as np
from scipy.optimize import least_squares

rng = np.random.default_rng(3)


def observed(theta):
    q0, q1, s1, c1, s2, c2, m00, m10, m01, m11 = theta   # m_ke = P(Y=1 | Z*=k, E=e)
    out = []
    for e, q in ((0, q0), (1, q1)):
        mu = {0: (m00, m01)[e], 1: (m10, m11)[e]}
        P = np.zeros((2, 2, 2))                          # P[v1, v2, y | e]
        for z, pz in ((1, q), (0, 1 - q)):
            p1 = s1 if z else 1 - c1
            p2 = s2 if z else 1 - c2
            for v1 in (0, 1):
                for v2 in (0, 1):
                    for y in (0, 1):
                        P[v1, v2, y] += (pz * (p1 if v1 else 1 - p1) * (p2 if v2 else 1 - p2)
                                         * (mu[z] if y else 1 - mu[z]))
        out.append(P)
    return out


def exact_fits(target, starts=40):
    fits = []
    resid = lambda x: np.concatenate([(a - b).ravel() for a, b in zip(observed(x), target)])
    for _ in range(starts):
        r = least_squares(resid, rng.uniform(0.05, 0.95, 10), bounds=(1e-6, 1 - 1e-6))
        if r.cost < 1e-20:
            x = r.x
            if x[2] + x[3] < 1:                          # relabel Z*
                x = np.array([1 - x[0], 1 - x[1], 1 - x[3], 1 - x[2], 1 - x[5], 1 - x[4],
                              x[7], x[6], x[9], x[8]])
            fits.append(x)
    return fits


worst = 0.0
for t in range(15):
    q0, q1 = rng.uniform(0.15, 0.85, 2)
    while abs(q0 - q1) < 0.15:
        q1 = rng.uniform(0.15, 0.85)
    s1, c1, s2, c2 = rng.uniform(0.65, 0.97, 4)
    m = rng.uniform(0.05, 0.95, 4)
    if t % 2 == 0:
        m[2], m[3] = m[0], m[1]                          # Assumption 4 holds
    theta = np.array([q0, q1, s1, c1, s2, c2, *m])
    fits = exact_fits(observed(theta))
    err = max(np.abs(x - theta).max() for x in fits)
    worst = max(worst, err)
    print(f"model {t:2d} ({'A4 holds' if t % 2 == 0 else 'A4 fails'}): "
          f"{len(fits)} exact fits, max deviation from truth {err:.1e}")
print(f"worst deviation over all models and fits: {worst:.1e}")

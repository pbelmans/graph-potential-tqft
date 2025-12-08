"""
To illustrate that how the TQFT approach to computing period sequences of graph potentials,
we implement Theorem C over a truncated power series ring.
"""

import functools
import itertools

N = 20

R = PowerSeriesRing(QQ, "t", default_prec=N)
t = R.gen()


def laplace(f):
    """Laplace transform of power series"""
    # t is actually t^2 in the original, hence we need factorial(2*e)
    d = f.degree()
    return R([factorial(2 * i) * f[i] for i in range(d + 1)]) + O(t ** (d + 1))


@functools.lru_cache(None)
def cached_binomial(n, k):
    """Cached version of binomial coefficients"""
    return binomial(n, k)


def entry(i, j):
    """Entry of the matrix A: series expansion from proof of Lemma 4.8"""
    # reindex around center of matrix
    i, j = i - (2 * N + 1), j - (2 * N + 1)

    if (i - j) % 2 != 0:
        return 0

    coefficients = [0] * (4 * N + 1)

    for a, b in itertools.product(range(2 * N + 1), repeat=2):
        c, d = a - i, b - j

        if c < 0 or c > 2 * N or d < 0 or d > 2 * N:
            continue
        if (a + b) % 2 != 0 or (c + d) % 2 != 0:
            continue

        n, m = (a + b) / 2, (c + d) / 2

        coefficients[n + m] += (
            cached_binomial(2 * n, n)
            * cached_binomial(2 * m, m)
            / (factorial(a) * factorial(b) * factorial(c) * factorial(d))
        )

    return R(coefficients) + O(t ** (N + 1))


"""The two linear operators from Theorem C"""
# we can precompute these, as they are independent of the genus
A = matrix(R, 4 * N + 1, 4 * N + 1, entry)
S = matrix(QQ, 4 * N + 1, 4 * N + 1, lambda i, j: 1 if (i + j) == (4 * N + 2) else 0)


def period(g, d=0):
    """Theorem C"""
    return R((A ** (g - 1) * S ** ((g + d + 1) % 2)).trace())


def regularized_period(g, d=0):
    """Regularized classical period from Theorem C"""
    return laplace(period(g, d))

"""
To illustrate that the TQFT approach in `tqft.sage` is much more efficient,
we try to compute 11 entries in the period sequence of the graph potential
for the ladder graph in genus `g=3`.

Already the 10th coefficient (equal to 435456000) takes almost a minute on a
MacBook Pro M4, using Sage 10.5

The memory usage also grows quickly, showing that this is not a viable approach,
already in the first non-trivial case of interest; because `g=2` is
the intersection of 2 quadrics in projective 5-space, for which other methods exist.
"""

R = LaurentPolynomialRing(QQbar, "a,b,c,d,e,f")
R.inject_variables()

# graph potential from ladder graph
# doubled edges labelled a,b,d,e, so c,f is perfect matching
W = (
    1 / (a * b * c)
    + (b * c) / a
    + (a * c) / b
    + (a * b) / c
    + c * d * e
    + c / (d * e)
    + d / (c * e)
    + e / (c * d)
    + a * b * f
    + a / (b * f)
    + b / (a * f)
    + f / (a * b)
    + d * e * f
    + d / (e * f)
    + e / (d * f)
    + f / (d * e)
)

S = PowerSeriesRing(ZZ, "t")
t = S.gen()
print(S([(W**i).constant_coefficient() for i in range(0, 11, 2)]) + O(t ^ 6))

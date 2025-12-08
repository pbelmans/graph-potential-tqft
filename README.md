# Code accompanying "Graph potentials and topological quantum field theories"

This code accompanies the paper
[_Graph potentials and topological quantum field theories_](https://arxiv.org/abs/2205.07244),
by Pieter Belmans, Sergey Galkin, and Swarnava Mukhopadhyay.

The main goal is to give a proof-of-concept implementation of
the approach to compute classical periods of graph potentials
using topological quantum field theories,
see Theorem C for details.
This is implemented in `tqft.sage`

For comparison we also give
the naive approach of taking constant coefficients of powers of graph potentials,
which is much slower.
Reported timings are taken from a MacBook Pro M4, using Sage 10.5.

In order to cite this code,
please use

```bibitex
@software{tqft-graph-potential,
  author = {Belmans, Pieter and Galkin, Sergey and Mukhopadhyay, Swarnava},
  title  = {Code accompanying ``Graph potentials and topological quantum field theories''},
  url    = {https://github.com/pbelmans/tqft-graph-potential},
  year   = {2025},
}
```

## The TQFT approach

The following methods constitute the public interface of `tqft.sage`:

* `period(g, d=0)`: computes the (unregularized) period in genus `g` and parity `d`
* `regularized_period(g, d=0)`: computes the (regularized) period in genus `g` and parity `d`

Note that we have performed a substitution, replacing `t^2` by `t`,
as the period sequence only has even terms.

The code computes the matrices:

* `A` (as an approximation over a truncated power series ring)
* `S` (a binary matrix, with coefficients in a truncated power series ring)

once, so that the partition function of the TQFT can be computed up to a certain precision
in a truncated power series ring using Theorem C.

```sage
sage: %time load("tqft.sage")
CPU times: user 5.06 s, sys: 6.12 ms, total: 5.07 s
Wall time: 5.07 s
sage: %time regularized_period(3, 1)
CPU times: user 2.48 s, sys: 5.25 ms, total: 2.48 s
Wall time: 2.48 s
1 + 384*t^2 + 23040*t^3 + 3265920*t^4 + 435456000*t^5 + 68263641600*t^6 + 11300889600000*t^7 + 1984905402480000*t^8 + 363141494876160000*t^9 + 68740862681387280384*t^10 + 13374642698336861061120*t^11 + 2662767003536289520206336*t^12 + 540523656892326293054668800*t^13 + 111563632918634635298003558400*t^14 + 23360913428124338852025263063040*t^15 + 4953724620693464623365974246716800*t^16 + 1062191426353323530320788867571138560*t^17 + 230019469144455151388113309541761843200*t^18 + 50253274294240010831049317252037540249600*t^19 + 11066677786514528833678926891485928296785920*t^20 + O(t^21)
```

It is expected that using a language more oriented towards high-performance
and applying various optimizations can significantly improve the computation.
This is however only a proof-of-concept.

The file `output.py` computes the period sequences as reported in `output.txt`,
and is part of a continuous integration on GitHub to showcase the code.

## Comparison to naive approach

To showcase the efficiency of the TQFT approach,
we also provide an implementation of the naive approach,
of taking the constant coefficient of powers of the graph potential.
This is implemented in `naive.tqft`.

```sage
sage: %time load("naive.sage")
Defining a, b, c, d, e, f
1 + 384*t^2 + 23040*t^3 + 3265920*t^4 + 435456000*t^5 + O(t^6)
CPU times: user 1min 11s, sys: 197 ms, total: 1min 11s
Wall time: 1min 11s
```

This approach quickly becomes infeasible, already taking almost a minute for the 10th coefficient,
for the first non-trivial interesting case, namely the (odd) graph potential in genus 3.

[![CI](https://github.com/pbelmans/graph-potential-tqft/actions/workflows/check.yml/badge.svg)](https://github.com/pbelmans/graph-potential-tqft/actions/workflows/check.yml)

# Code accompanying "Graph potentials and topological quantum field theories"

This code accompanies the paper
[_Graph potentials and topological quantum field theories_](https://arxiv.org/abs/2205.07244),
by Pieter Belmans, Sergey Galkin, and Swarnava Mukhopadhyay,
accepted for publication in Proceedings of the London Mathematical Society.

The main goal is to implement the approach to compute classical periods of graph potentials
using topological quantum field theories,
see Theorem C for details.
There are two implementations: one in Sage (`tqft.sage`) and one in Julia (`tqft.jl`).

For comparison we also give the naive approach of taking constant coefficients
of powers of graph potentials, which is much slower.
Reported timings are taken from a MacBook Pro M4, using Sage 10.5 and Julia 1.11.

In order to cite this code,
please use

```bibtex
@software{tqft-graph-potential,
  author = {Belmans, Pieter and Galkin, Sergey and Mukhopadhyay, Swarnava},
  title  = {Code accompanying ``Graph potentials and topological quantum field theories''},
  url    = {https://github.com/pbelmans/graph-potential-tqft},
  year   = {2026},
}
```

## The TQFT approach: Sage

The following methods constitute the public interface of `tqft.sage`:

* `period(g, d=0)`: computes the (unregularized) period in genus `g` and parity `d`
* `regularized_period(g, d=0)`: computes the (regularized) period in genus `g` and parity `d`

Note that we have performed a substitution, replacing `t^2` by `t`,
as the period sequence only has even terms.

The code computes the matrices:

* `A` (as an approximation over a truncated power series ring)
* `S` (an anti-diagonal 0/1 matrix, with entries in a truncated power series ring)

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

The file `output.sage` computes the period sequences as reported in `output.txt`,
and is part of a continuous integration on GitHub to check the code.
It computes the periods in Tables 1 and 2 of the paper,
at up to twice the precision reported in the paper.

## Comparison to naive approach: Sage

To showcase the efficiency of the TQFT approach,
we also provide an implementation of the naive approach in Sage,
of taking the constant coefficient of powers of the graph potential.
This is implemented in `naive.sage`.

```sage
sage: %time load("naive.sage")
Defining a, b, c, d, e, f
1 + 384*t^2 + 23040*t^3 + 3265920*t^4 + 435456000*t^5 + O(t^6)
CPU times: user 1min 11s, sys: 197 ms, total: 1min 11s
Wall time: 1min 11s
```

This approach quickly becomes infeasible for the first non-trivial case of interest,
namely the (odd) graph potential in genus 3: computing only 6 coefficients already takes over a minute.

## The TQFT approach: Julia

The Julia implementation in `tqft.jl` is a reimplementation of `tqft.sage` using
[Nemo.jl](https://nemocas.org/), which wraps the FLINT library for fast exact arithmetic.
The public interface mirrors the Sage version:

```julia
julia> include("tqft.jl")
julia> @time tqft = TQFT(20);
  0.008748 seconds (165.73 k allocations: 6.157 MiB)
julia> @time regularized_period(tqft, 3, 1)
  0.066633 seconds (870.50 k allocations: 37.147 MiB, 11.78% gc time)
1 + 384*t^2 + 23040*t^3 + 3265920*t^4 + 435456000*t^5 + 68263641600*t^6 + 11300889600000*t^7 + 1984905402480000*t^8 + 363141494876160000*t^9 + 68740862681387280384*t^10 + 13374642698336861061120*t^11 + 2662767003536289520206336*t^12 + 540523656892326293054668800*t^13 + 111563632918634635298003558400*t^14 + 23360913428124338852025263063040*t^15 + 4953724620693464623365974246716800*t^16 + 1062191426353323530320788867571138560*t^17 + 230019469144455151388113309541761843200*t^18 + 50253274294240010831049317252037540249600*t^19 + 11066677786514528833678926891485928296785920*t^20 + O(t^21)
```

Timings are reported after JIT compilation (i.e., after an initial warmup run).

The speedup over Sage is largely due to FLINT's optimized polynomial arithmetic and
a block-diagonal decomposition of the matrix `A`: since `A[i,j] = 0` whenever `i - j`
is odd, the matrix splits into two independent blocks `A_even` and `A_odd`, reducing
the cost of each matrix power from roughly $(4N+1)^3$ to $(2N+1)^3 + (2N)^3$,
a factor of about 4.

## Comparison to naive approach: Julia

We also provide a Julia reimplementation of the naive approach in `naive.jl`.
Instead of working in a multivariate polynomial ring, it represents the graph potential
as a dictionary mapping exponent tuples to integer coefficients, using Nemo's
arbitrary-precision integers. This avoids symbolic ring overhead and keeps the
representation sparse throughout.

```julia
julia> include("naive.jl")
julia> @time naive_periods(6)
  1.737 seconds (28185699 allocations: 513.766 MiB)
1 + 384*t^2 + 23040*t^3 + 3265920*t^4 + 435456000*t^5 + O(t^6)
```

This is already significantly faster than the Sage naive implementation (~71 seconds),
thanks to working directly with sparse dictionary representations
and Nemo's exact integer arithmetic.
Nevertheless, the approach scales very poorly with the number of terms,
confirming that the TQFT method is essential for computing longer period sequences.

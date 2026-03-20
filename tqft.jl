"""
TQFT approach to computing period sequences of graph potentials.
Translation of tqft.sage using Nemo.jl with block-diagonal decomposition.

Usage:
  tqft = TQFT(20);
  regularized_period(tqft, 3, 1)
"""

using Nemo: QQ, ZZ, QQRelPowerSeriesRing, QQRelPowerSeriesRingElem, MatElem, ZZRingElem, QQFieldElem, power_series_ring, matrix, binomial, precision, coeff, parent, tr

struct TQFT
  ring::QQRelPowerSeriesRing
  variable::QQRelPowerSeriesRingElem
  A_even::MatElem{QQRelPowerSeriesRingElem}
  A_odd::MatElem{QQRelPowerSeriesRingElem}
  precision::Int
end

"""
Build the TQFT block-decomposed matrices for precision `N`.

A has block structure A = A_even ⊕ A_odd because entry(i,j)=0 when i-j is odd.
This gives ~4× speedup on matrix multiply: (2N+1)³ + (2N)³ vs (4N+1)³.
"""
function TQFT(N::Int = 20)
  prec = N + 1
  ring, variable = power_series_ring(QQ, prec, "t")

  max_k = 4N
  factorial_table = Vector{ZZRingElem}(undef, max_k + 1)
  factorial_table[1] = ZZ(1)
  for k in 1:max_k
    factorial_table[k + 1] = factorial_table[k] * k
  end

  inverse_factorial = [QQ(1, factorial_table[k + 1]) for k in 0:max_k]
  central_binomial = [binomial(ZZ(2n), ZZ(n)) for n in 0:2N]

  # Weight table: weight[x+1, y+1] = C(x+y, (x+y)/2) / (x! * y!)
  weight = Matrix{QQFieldElem}(undef, 2N + 1, 2N + 1)
  for x in 0:2N, y in 0:2N
    weight[x + 1, y + 1] = iseven(x + y) ?
      QQ(central_binomial[(x + y) >> 1 + 1]) * inverse_factorial[x + 1] * inverse_factorial[y + 1] :
      QQ(0)
  end

  function make_entry(i_orig::Int, j_orig::Int)
    ri, rj = i_orig - (2N + 1), j_orig - (2N + 1)
    coefficients = [QQ(0) for _ in 1:prec]

    for a in 0:2N
      c = a - ri
      (0 <= c <= 2N) || continue
      for b in (a & 1):2:2N
        d = b - rj
        (0 <= d <= 2N) || continue
        idx = (a + b) >> 1 + (c + d) >> 1
        idx < prec || continue
        coefficients[idx + 1] += weight[a + 1, b + 1] * weight[c + 1, d + 1]
      end
    end

    ring(coefficients, prec, prec, 0)
  end

  # Even block: indices 0, 2, 4, ..., 4N → dim = 2N+1
  # Odd block:  indices 1, 3, 5, ..., 4N-1 → dim = 2N
  dim_even, dim_odd = 2N + 1, 2N
  A_even = matrix(ring, dim_even, dim_even,
    [make_entry(2(k - 1), 2(l - 1)) for k in 1:dim_even for l in 1:dim_even])
  A_odd = matrix(ring, dim_odd, dim_odd,
    [make_entry(2(k - 1) + 1, 2(l - 1) + 1) for k in 1:dim_odd for l in 1:dim_odd])

  TQFT(ring, variable, A_even, A_odd, N)
end

"""
Anti-diagonal trace: compute trace(M · S) where S is the anti-diagonal matrix
with S[i,j] = 1 iff i+j = 4N+2 (0-based). Avoids forming S or doing matrix multiply.
"""
function antidiag_trace(M_even, M_odd, N::Int)
  sum(M_even[p, 2N + 3 - p] for p in 2:(2N + 1)) + sum(M_odd[p, 2N + 2 - p] for p in 2:2N)
end

"""Compute the (unregularized) period via Theorem C."""
function period(tqft::TQFT, g::Int, d::Int = 0)
  M_even = tqft.A_even^(g - 1)
  M_odd  = tqft.A_odd^(g - 1)

  iseven(g + d + 1) && return tr(M_even) + tr(M_odd)
  antidiag_trace(M_even, M_odd, tqft.precision)
end

"""Laplace transform: multiply i-th coefficient by (2i)!."""
function laplace(f)
  prec = precision(f)
  parent(f)(QQFieldElem[factorial(ZZ(2i)) * coeff(f, i) for i in 0:prec-1], prec, prec, 0)
end

"""Regularized classical period from Theorem C."""
function regularized_period(tqft::TQFT, g::Int, d::Int = 0)
  return laplace(period(tqft, g, d))
end

nothing

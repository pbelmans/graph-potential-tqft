"""
Naive approach: compute period sequence by taking constant coefficients
of powers of the graph potential (ladder graph, genus 3).

Uses a dictionary-based Laurent polynomial over ZZ for fast integer arithmetic.
This is inherently slow for large powers — the TQFT approach is much more efficient.
"""

using Nemo

const LaurentPoly = Dict{NTuple{6,Int}, ZZRingElem}

function lpoly_mul(f::LaurentPoly, g::LaurentPoly)::LaurentPoly
    result = LaurentPoly()
    for (e1, c1) in f, (e2, c2) in g
        e = (e1[1]+e2[1], e1[2]+e2[2], e1[3]+e2[3], e1[4]+e2[4], e1[5]+e2[5], e1[6]+e2[6])
        result[e] = get!(result, e, ZZ(0)) + c1 * c2
    end
    return result
end

function const_coeff(f::LaurentPoly)::ZZRingElem
    return get(f, (0,0,0,0,0,0), ZZ(0))
end

"""
Compute the first `n_terms` of the period sequence for the ladder graph (genus 3, odd).
Returns a Nemo power series: Σ const_coeff(W^{2i}) * t^i for i = 0, ..., n_terms-1.
"""
function naive_periods(n_terms::Int = 6)
    one = ZZ(1)
    W = LaurentPoly()
    terms = [
        (-1, -1, -1,  0,  0,  0),  # 1/(abc)
        (-1,  1,  1,  0,  0,  0),  # bc/a
        ( 1, -1,  1,  0,  0,  0),  # ac/b
        ( 1,  1, -1,  0,  0,  0),  # ab/c
        ( 0,  0,  1,  1,  1,  0),  # cde
        ( 0,  0,  1, -1, -1,  0),  # c/(de)
        ( 0,  0, -1,  1, -1,  0),  # d/(ce)
        ( 0,  0, -1, -1,  1,  0),  # e/(cd)
        ( 1,  1,  0,  0,  0,  1),  # abf
        ( 1, -1,  0,  0,  0, -1),  # a/(bf)
        (-1,  1,  0,  0,  0, -1),  # b/(af)
        (-1, -1,  0,  0,  0,  1),  # f/(ab)
        ( 0,  0,  0,  1,  1,  1),  # def
        ( 0,  0,  0,  1, -1, -1),  # d/(ef)
        ( 0,  0,  0, -1,  1, -1),  # e/(df)
        ( 0,  0,  0, -1, -1,  1),  # f/(de)
    ]
    for t in terms
        W[t] = get(W, t, ZZ(0)) + one
    end

    R, t = power_series_ring(ZZ, n_terms, "t")
    coeffs = ZZRingElem[]

    # Iteratively compute W^0, W^2, W^4, ... by multiplying by W^2 each time
    W2 = lpoly_mul(W, W)
    current = LaurentPoly((0,0,0,0,0,0) => ZZ(1))

    for i in 0:n_terms-1
        push!(coeffs, const_coeff(current))
        if i < n_terms - 1
            current = lpoly_mul(current, W2)
        end
    end

    return R(coeffs, length(coeffs), n_terms, 0)
end

if abspath(PROGRAM_FILE) == @__FILE__
    println("Computing naive periods (6 terms)...")
    @time result = naive_periods(6)
    println(result)
end

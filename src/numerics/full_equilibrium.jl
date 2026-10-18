module FullEquilibrium

using ArbNumerics
using ..IntervalSafety
const IS = IntervalSafety

"""
Interval Gaussian elimination of the ENTIRE rectangular augmented system
`transpose(A) * lambda = objective_gradient`.

Returns `(:ok, multipliers)`, `(:stationarity, nothing)` for proved
incompatibility, or `(:rank, nothing)` when verification fails. An `:ok`
result encloses every possible solution; it does not prove existence.

All rows participate in elimination. Verified nonzero interval pivots are
required. Remaining rows are tested for incompatibility, then back substitution
encloses the multipliers. Failure to verify a pivot is inconclusive, not proof
of incompatibility. Inputs are not mutated.
"""
function enclose_multipliers(::Type{T}, A, objective_gradient) where {T<:ArbReal}
    unknowns, equations = size(A)
    length(objective_gradient) == equations || throw(DimensionMismatch())
    equations >= unknowns || return (:rank, nothing)
    all(IS.finite_interval, A) && all(IS.finite_interval, objective_gradient) ||
        return (:rank, nothing)
    work = Matrix{T}(undef, equations, unknowns + 1)
    for row in 1:equations
        for column in 1:unknowns
            work[row, column] = A[column, row]
        end
        work[row, unknowns + 1] = objective_gradient[row]
    end

    # Keep the augmented matrix and its native numeric objects alive throughout
    # elimination and back substitution, including temporary scalar C calls.
    GC.@preserve work begin
        for column in 1:unknowns
            pivot_row = 0
            best_strength = -Inf
            for row in column:equations
                value = work[row, column]
                IS.zero_excluded(value) || continue
                strength = Float64(IS.interval_min(abs(lowerbound(value)), abs(upperbound(value))))
                if strength > best_strength
                    pivot_row = row
                    best_strength = strength
                end
            end
            pivot_row == 0 && return (:rank, nothing)
            if pivot_row != column
                for entry in 1:(unknowns + 1)
                    work[column, entry], work[pivot_row, entry] =
                        work[pivot_row, entry], work[column, entry]
                end
            end
            for row in (column + 1):equations
                iszero(work[row, column]) && continue
                factor = work[row, column] / work[column, column]
                IS.finite_interval(factor) || return (:rank, nothing)
                for entry in (column + 1):(unknowns + 1)
                    work[row, entry] -= factor * work[column, entry]
                    IS.finite_interval(work[row, entry]) || return (:rank, nothing)
                end
                # This coefficient is algebraically zero for every point system:
                # a - (a/p)*p = 0. Recording that identity avoids a spurious radius.
                work[row, column] = T(0)
            end
        end

        # All unknown coefficients in these rows were eliminated. A right-hand
        # side excluding zero proves the full system incompatible throughout B.
        for row in (unknowns + 1):equations
            IS.zero_excluded(work[row, unknowns + 1]) &&
                return (:stationarity, nothing)
        end

        multipliers = Vector{T}(undef, unknowns)
        for row in unknowns:-1:1
            right = work[row, unknowns + 1]
            for column in (row + 1):unknowns
                right -= work[row, column] * multipliers[column]
            end
            multipliers[row] = right / work[row, row]
            IS.finite_interval(multipliers[row]) || return (:rank, nothing)
        end
        return (:ok, multipliers)
    end
end

end

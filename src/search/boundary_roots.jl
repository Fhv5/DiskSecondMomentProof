module BoundaryRoots

export BoundaryRootWitness,
       VerifiedBoundaryRoot,
       registered_witnesses,
       verify_witness,
       witness_state_in_box,
       rational_text,
       qsqrt3_text

using ArbNumerics
using ...DiskSecondMomentProof
using ..LeafFilters

const LE = LeafFilters
const Q = DiskSecondMomentProof.Qsqrt3
const RationalBig = Rational{BigInt}

"""
Data in Q(sqrt(3)) proposed for a critical root on a contact-stratum boundary.
`permutation[old]` is the corresponding label in the target contact graph.
"""
struct BoundaryRootWitness
    graph6::String
    target_graph6::String
    angle_pi_multiples::Vector{RationalBig}
    tree_steps::Vector{Tuple{Q,Q}}
    reduced_multipliers::Vector{Q}
    extra_contacts::Vector{Tuple{Int,Int}}
    permutation::Vector{Int}
    expected_moment::Q
    expected_determinant::Q
end

"""A boundary-root witness whose equations have been checked."""
struct VerifiedBoundaryRoot
    witness::BoundaryRootWitness
    determinant::Q
end

"""Format a rational number without floating-point conversion."""
rational_text(value::Rational) =
    denominator(value) == 1 ? string(numerator(value)) :
    "$(numerator(value))/$(denominator(value))"

"""Format a Q(sqrt(3)) number without floating-point conversion."""
function qsqrt3_text(value::Q)
    iszero(value.b) && return rational_text(value.a)
    iszero(value.a) && return "$(rational_text(value.b))*sqrt(3)"
    sign = value.b > 0 ? "+" : "-"
    return "$(rational_text(value.a)) $sign $(rational_text(abs(value.b)))*sqrt(3)"
end

"""Return the registered witnesses applicable to one graph."""
function registered_witnesses(graph6::AbstractString)
    graph6 = String(graph6)
    sqrt_three = Q(0 // 1, 1 // 1)
    if graph6 == "EPVW"
        return BoundaryRootWitness[BoundaryRootWitness(
            "EPVW",
            "EQlw",
            RationalBig[2 // 3, 2 // 3, 2 // 3, 1 // 3],
            [
                (Q(-1), sqrt_three),
                (Q(-1), sqrt_three),
                (Q(-1), sqrt_three),
                (Q(1), sqrt_three),
            ],
            Q[Q(1 // 2), Q(1), Q(1 // 2)],
            [(1, 4)],
            [4, 1, 2, 6, 3, 5],
            Q(11),
            Q(-4608),
        )]
    end
    return BoundaryRootWitness[]
end

function tree_step(angle_pi_multiple::RationalBig)
    third_turns = angle_pi_multiple * 3
    denominator(third_turns) == 1 ||
        error("unsupported boundary-witness angle: expected a multiple of pi/3")
    direction = mod(numerator(third_turns), 6)
    sqrt_three = Q(0 // 1, 1 // 1)
    direction == 0 && return (Q(2), Q(0))
    direction == 1 && return (Q(1), sqrt_three)
    direction == 2 && return (Q(-1), sqrt_three)
    direction == 3 && return (Q(-2), Q(0))
    direction == 4 && return (Q(-1), -sqrt_three)
    return (Q(1), -sqrt_three)
end

"""Compute a determinant over Q(sqrt(3))."""
function qsqrt3_determinant(matrix::Matrix{Q})
    rows, columns = size(matrix)
    rows == columns || error("determinant requires a square matrix")
    rows == 0 && return Q(1)
    work = copy(matrix)
    denominator = Q(1)
    sign = Q(1)
    for pivot_index in 1:(rows - 1)
        relative_row = findfirst(
            row -> !iszero(work[row, pivot_index]),
            pivot_index:rows,
        )
        relative_row === nothing && return Q(0)
        pivot_row = pivot_index - 1 + relative_row
        if pivot_row != pivot_index
            work[pivot_index, :], work[pivot_row, :] =
                copy(work[pivot_row, :]), copy(work[pivot_index, :])
            sign = -sign
        end
        pivot = work[pivot_index, pivot_index]
        for row in (pivot_index + 1):rows, column in (pivot_index + 1):columns
            work[row, column] =
                (work[row, column] * pivot -
                 work[row, pivot_index] * work[pivot_index, column]) / denominator
        end
        denominator = pivot
        for row in (pivot_index + 1):rows
            work[row, pivot_index] = Q(0)
        end
    end
    return sign * work[rows, columns]
end

"""Verify the equations, contact inequalities and determinant of a boundary-root witness."""
function verify_witness(graph, tree, witness::BoundaryRootWitness)
    graph.code == witness.graph6 || error("boundary witness graph mismatch")
    length(tree) == graph.n - 1 || error("invalid rooted tree")
    length(witness.angle_pi_multiples) == length(tree) - 1 ||
        error("boundary witness angle count mismatch")
    length(witness.tree_steps) == length(tree) - 1 ||
        error("boundary witness step count mismatch")
    for index in eachindex(witness.tree_steps)
        witness.tree_steps[index] == tree_step(witness.angle_pi_multiples[index]) ||
            error("boundary witness angle and tree step disagree at index $index")
    end
    k = length(witness.tree_steps)
    raw = fill(Q(0), 2, graph.n)
    first_derivative = fill(Q(0), 2, graph.n, k)
    second_diagonal = fill(Q(0), 2, graph.n, k)
    raw[1, tree[1][2]] = Q(2)
    for (angle_index, (parent, child)) in enumerate(tree[2:end])
        dx, dy = witness.tree_steps[angle_index]
        dx * dx + dy * dy == Q(4) || error("boundary witness step is not a contact")
        raw[1, child] = raw[1, parent] + dx
        raw[2, child] = raw[2, parent] + dy
        first_derivative[:, child, :] .= first_derivative[:, parent, :]
        second_diagonal[:, child, :] .= second_diagonal[:, parent, :]
        first_derivative[1, child, angle_index] += -dy
        first_derivative[2, child, angle_index] += dx
        second_diagonal[1, child, angle_index] += -dx
        second_diagonal[2, child, angle_index] += -dy
    end

    centers = copy(raw)
    center_first = copy(first_derivative)
    center_second = copy(second_diagonal)
    for coordinate in 1:2
        centers[coordinate, :] .-= sum(raw[coordinate, :]) / graph.n
        for angle_index in 1:k
            center_first[coordinate, :, angle_index] .-=
                sum(first_derivative[coordinate, :, angle_index]) / graph.n
            center_second[coordinate, :, angle_index] .-=
                sum(second_diagonal[coordinate, :, angle_index]) / graph.n
        end
    end

    tree_edges = Set((min(left, right), max(left, right)) for (left, right) in tree)
    extra_edges = [edge for edge in graph.edges if !(edge in tree_edges)]
    q = length(extra_edges)
    length(witness.reduced_multipliers) == q ||
        error("boundary witness multiplier count mismatch")
    constraints = fill(Q(0), q)
    constraint_jacobian = fill(Q(0), q, k)
    constraint_hessians = fill(Q(0), q, k, k)
    for (row, (left, right)) in enumerate(extra_edges)
        dx = raw[1, left] - raw[1, right]
        dy = raw[2, left] - raw[2, right]
        constraints[row] = (dx * dx + dy * dy - Q(4)) / 2
        for first_angle in 1:k
            dax = first_derivative[1, left, first_angle] -
                  first_derivative[1, right, first_angle]
            day = first_derivative[2, left, first_angle] -
                  first_derivative[2, right, first_angle]
            constraint_jacobian[row, first_angle] = dx * dax + dy * day
            for second_angle in 1:k
                dbx = first_derivative[1, left, second_angle] -
                      first_derivative[1, right, second_angle]
                dby = first_derivative[2, left, second_angle] -
                      first_derivative[2, right, second_angle]
                value = dax * dbx + day * dby
                if first_angle == second_angle
                    ddx = second_diagonal[1, left, first_angle] -
                          second_diagonal[1, right, first_angle]
                    ddy = second_diagonal[2, left, first_angle] -
                          second_diagonal[2, right, first_angle]
                    value += dx * ddx + dy * ddy
                end
                constraint_hessians[row, first_angle, second_angle] = value
            end
        end
    end

    gradient = fill(Q(0), k)
    hessian = fill(Q(0), k, k)
    for first_angle in 1:k
        gradient[first_angle] = sum(
            centers[coordinate, vertex] *
            center_first[coordinate, vertex, first_angle]
            for coordinate in 1:2, vertex in 1:graph.n
        )
        for second_angle in 1:k
            value = sum(
                center_first[coordinate, vertex, first_angle] *
                center_first[coordinate, vertex, second_angle]
                for coordinate in 1:2, vertex in 1:graph.n
            )
            if first_angle == second_angle
                value += sum(
                    centers[coordinate, vertex] *
                    center_second[coordinate, vertex, first_angle]
                    for coordinate in 1:2, vertex in 1:graph.n
                )
            end
            hessian[first_angle, second_angle] = value
        end
    end

    residual = vcat(
        constraints,
        gradient - transpose(constraint_jacobian) * witness.reduced_multipliers,
    )
    all(iszero, residual) || error("boundary witness has a nonzero KKT residual")
    lagrangian_hessian = copy(hessian)
    for edge_index in 1:q
        lagrangian_hessian .-=
            witness.reduced_multipliers[edge_index] .*
            constraint_hessians[edge_index, :, :]
    end
    jacobian = fill(Q(0), q + k, k + q)
    jacobian[1:q, 1:k] .= constraint_jacobian
    jacobian[(q + 1):(q + k), 1:k] .= lagrangian_hessian
    jacobian[(q + 1):(q + k), (k + 1):(k + q)] .= -transpose(constraint_jacobian)
    determinant = qsqrt3_determinant(jacobian)
    determinant == witness.expected_determinant ||
        error("boundary witness KKT determinant mismatch")
    !iszero(determinant) || error("boundary witness KKT root is singular")

    squared_distance(left, right) =
        (centers[1, left] - centers[1, right])^2 +
        (centers[2, left] - centers[2, right])^2
    all(squared_distance(left, right) == Q(4) for (left, right) in graph.edges) ||
        error("boundary witness misses a declared contact")
    all(squared_distance(left, right) == Q(4)
        for (left, right) in witness.extra_contacts) ||
        error("boundary witness misses an asserted additional contact")
    moment = sum(value * value for value in centers) / 2
    moment == witness.expected_moment || error("boundary witness second moment mismatch")

    target_graph = LE.graph_data(witness.target_graph6)
    augmented_edges = vcat(graph.edges, witness.extra_contacts)
    mapped_edges = Set(
        (min(witness.permutation[left], witness.permutation[right]),
         max(witness.permutation[left], witness.permutation[right]))
        for (left, right) in augmented_edges
    )
    mapped_edges == Set(target_graph.edges) ||
        error("boundary witness augmented graph mismatch")
    return VerifiedBoundaryRoot(witness, determinant)
end

"""Enclose one Q(sqrt(3)) number in an Arb interval."""
function arb_qsqrt3(::Type{T}, value::Q) where {T<:ArbReal}
    safety = DiskSecondMomentProof.IntervalSafety
    return safety.enclosing_rational(T, value.a) +
           safety.enclosing_rational(T, value.b) * sqrt(T(3))
end

"""Return true when a verified witness state lies in an interval state box."""
function witness_state_in_box(
    ::Type{T},
    verified::VerifiedBoundaryRoot,
    angle_box,
    multiplier_box,
) where {T<:ArbReal}
    witness = verified.witness
    length(witness.angle_pi_multiples) == length(angle_box) || return false
    length(witness.reduced_multipliers) == length(multiplier_box) || return false
    for index in eachindex(angle_box)
        angle_enclosure = DiskSecondMomentProof.IntervalSafety.enclosing_rational(
            T, witness.angle_pi_multiples[index]) * T(pi)
        lowerbound(angle_enclosure) >= lowerbound(angle_box[index]) || return false
        upperbound(angle_enclosure) <= upperbound(angle_box[index]) || return false
    end
    for index in eachindex(multiplier_box)
        multiplier_enclosure = arb_qsqrt3(T, witness.reduced_multipliers[index])
        lowerbound(multiplier_enclosure) >= lowerbound(multiplier_box[index]) || return false
        upperbound(multiplier_enclosure) <= upperbound(multiplier_box[index]) || return false
    end
    return true
end

end

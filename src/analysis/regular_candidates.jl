using LinearAlgebra
using ArbNumerics

"""
Result of the interval route for one regular isolated candidate.

All fields used for mathematical decisions retain `ArbReal` enclosures. Any
floating-point linear algebra used internally is only allowed to propose a
pivot set or a witness; the proposal is accepted only after an interval check.
"""
struct IntrinsicHessianVerification{T<:ArbReal}
    status::Symbol
    kernel_dimension::Int
    internal_dimension::Int
    contact_pivot_columns::Vector{Int}
    rigid_coordinate_rows::Vector{Int}
    H_kernel::Matrix{T}
    H_internal::Matrix{T}
    negative_direction::Vector{T}
    negative_value::Union{Nothing,T}
end

struct RegularCandidateResult{T<:ArbReal}
    candidate_id::String
    name::String
    graph6::String
    n::Int
    m::Int
    precision::Int
    isolation_status::Symbol
    graph_geometry_status::Symbol
    minimum_noncontact_gap::Union{Nothing,T}
    rank_status::Symbol
    rank_A::Union{Nothing,Int}
    gauge_status::Symbol
    leaf_status::Symbol
    leaf_vertex::Union{Nothing,Int}
    leaf_witness::Union{Nothing,T}
    multiplier_status::Symbol
    lambda_box::Vector{T}
    hessian::Union{Nothing,IntrinsicHessianVerification{T}}
    objective_interval::Union{Nothing,T}
    objective_status::Symbol
    conclusion::Symbol
    notes::Vector{String}
end

@inline _contains_zero(x::ArbReal) = !IntervalSafety.zero_excluded(x)
@inline _strictly_positive(x::ArbReal) =
    IntervalSafety.finite_interval(x) && lowerbound(x) > zero(x)
@inline _strictly_negative(x::ArbReal) =
    IntervalSafety.finite_interval(x) && upperbound(x) < zero(x)

function _midpoint_float_matrix(A::AbstractMatrix{<:ArbReal})
    return [Float64(midpoint(A[i, j])) for i in axes(A, 1), j in axes(A, 2)]
end

function _plain_transpose(A::AbstractMatrix{T}) where {T}
    result = zeros(T, size(A, 2), size(A, 1))
    for i in axes(A, 1), j in axes(A, 2)
        result[j, i] = A[i, j]
    end
    return result
end

function _numeric_pivot_columns(A::AbstractMatrix{<:ArbReal}, target_rank::Int)
    M = _midpoint_float_matrix(A)
    rows, cols = size(M)
    pivot_columns = Int[]
    pivot_row = 1
    for col in 1:cols
        pivot_row > rows && break
        candidates = pivot_row:rows
        local_index = argmax(abs.(M[candidates, col]))
        selected_row = first(candidates) + local_index - 1
        pivot = M[selected_row, col]
        pivot == 0.0 && continue
        if selected_row != pivot_row
            M[pivot_row, :], M[selected_row, :] = copy(M[selected_row, :]), copy(M[pivot_row, :])
        end
        for row in (pivot_row + 1):rows
            factor = M[row, col] / M[pivot_row, col]
            M[row, col:cols] .-= factor .* M[pivot_row, col:cols]
        end
        push!(pivot_columns, col)
        pivot_row += 1
        length(pivot_columns) == target_rank && break
    end
    return pivot_columns
end

function _verified_determinant(A::AbstractMatrix{T}, columns::Vector{Int}) where {T<:ArbReal}
    isempty(columns) && return true, one(T)
    minor = Matrix{T}(A[:, columns])
    size(minor, 1) == size(minor, 2) || return false, zero(T)
    verified = IntervalSafety.verified_inverse(minor)
    verified === nothing && return false, zero(T)
    return true, last(verified)
end

function _first_verified_column_minor(A::AbstractMatrix{T}, rank::Int) where {T<:ArbReal}
    rows, cols = size(A)
    all(IntervalSafety.finite_interval, A) || return nothing, nothing
    rank == rows || return nothing, nothing
    rank == 0 && return Int[], one(T)
    rank > cols && return nothing, nothing

    proposed = _numeric_pivot_columns(A, rank)
    if length(proposed) == rank
        valid, determinant_box = _verified_determinant(A, proposed)
        valid && return proposed, determinant_box
    end

    choice = Vector{Int}(undef, rank)
    found_columns = Ref{Union{Nothing,Vector{Int}}}(nothing)
    found_det = Ref{Union{Nothing,T}}(nothing)

    function search(start::Int, depth::Int)
        found_columns[] === nothing || return
        if depth > rank
            columns = copy(choice)
            valid, determinant_box = _verified_determinant(A, columns)
            if valid
                found_columns[] = columns
                found_det[] = determinant_box
            end
            return
        end
        last_start = cols - (rank - depth)
        for col in start:last_start
            choice[depth] = col
            search(col + 1, depth + 1)
            found_columns[] === nothing || return
        end
    end

    search(1, 1)
    return found_columns[], found_det[]
end

function _interval_kernel_basis(A::AbstractMatrix{T}, pivot_columns::Vector{Int}) where {T<:ArbReal}
    rows, cols = size(A)
    length(pivot_columns) == rows || throw(ArgumentError("pivot set does not prove full row rank"))
    pivot_set = Set(pivot_columns)
    free_columns = [col for col in 1:cols if !(col in pivot_set)]
    K = zeros(T, cols, length(free_columns))
    for (j, col) in enumerate(free_columns)
        K[col, j] = one(T)
    end
    if !isempty(free_columns)
        pivot_block = Matrix{T}(A[:, pivot_columns])
        free_block = Matrix{T}(A[:, free_columns])
        verified = IntervalSafety.verified_inverse(pivot_block)
        verified === nothing && throw(ArgumentError("interval kernel pivot is unverified"))
        K[pivot_columns, :] .= -(first(verified) * free_block)
    end
    return K, free_columns
end

function _prove_interval_positive_definite(H::AbstractMatrix{T}) where {T<:ArbReal}
    n = size(H, 1)
    size(H, 2) == n || return false
    n == 0 && return true
    S = (Matrix{T}(H) + _plain_transpose(H)) / T(2)
    L = zeros(T, n, n)
    for j in 1:n
        diagonal = S[j, j]
        for k in 1:(j - 1)
            diagonal -= L[j, k]^2
        end
        IntervalSafety.finite_interval(diagonal) &&
            lowerbound(diagonal) > 0 || return false
        L[j, j] = sqrt(diagonal)
        for i in (j + 1):n
            value = S[i, j]
            for k in 1:(j - 1)
                value -= L[i, k] * L[j, k]
            end
            L[i, j] = value / L[j, j]
        end
    end
    return true
end

function _find_interval_negative_witness(H::AbstractMatrix{T}) where {T<:ArbReal}
    n = size(H, 1)
    n == 0 && return T[], nothing
    all(IntervalSafety.finite_interval, H) || return T[], nothing
    midpoint_matrix = _midpoint_float_matrix((Matrix{T}(H) + _plain_transpose(H)) / T(2))
    all(isfinite, midpoint_matrix) || return T[], nothing
    spectral = eigen(Symmetric(midpoint_matrix))
    order = sortperm(spectral.values)
    for index in order
        proposal = [T(value) for value in spectral.vectors[:, index]]
        q = dot(proposal, H * proposal)
        if upperbound(q) < 0
            return proposal, q
        end
    end
    return T[], nothing
end

function _lagrangian_hessian_at_candidate(obj::AbstractObjective, c::Configuration{T},
                                 g::ContactGraph, lambda::Vector{T}) where {T<:ArbReal}
    return lagrangian_hessian(ambient_hessian(obj, c), lambda, c, g)
end

function _verify_intrinsic_hessian(obj::AbstractObjective, c::Configuration{T},
                              g::ContactGraph, lambda::Vector{T},
                              contact_pivots::Vector{Int}) where {T<:ArbReal}
    A = contact_matrix(c, g)
    K, free_columns = _interval_kernel_basis(A, contact_pivots)
    kernel_dim = size(K, 2)

    if kernel_dim < 3
        return IntrinsicHessianVerification{T}(
            :invalid_rigid_dimension, kernel_dim, 0, contact_pivots, Int[],
            zeros(T, kernel_dim, kernel_dim), zeros(T, 0, 0), T[], nothing
        )
    end

    rigid_motions = euclidean_motions(c)
    rigid_coordinates = Matrix{T}(rigid_motions[free_columns, :])
    rigid_rows, _ = _first_verified_column_minor(_plain_transpose(rigid_coordinates), 3)
    if rigid_rows === nothing
        return IntrinsicHessianVerification{T}(
            :rigid_modes_not_separated, kernel_dim, max(0, kernel_dim - 3),
            contact_pivots, Int[], zeros(T, kernel_dim, kernel_dim),
            zeros(T, 0, 0), T[], nothing
        )
    end

    rigid_row_set = Set(rigid_rows)
    internal_indices = [i for i in 1:kernel_dim if !(i in rigid_row_set)]
    internal_dim = length(internal_indices)

    H_lagrange = _lagrangian_hessian_at_candidate(obj, c, g, lambda)
    H_kernel = Matrix{T}(transpose(K) * H_lagrange * K)
    H_internal = Matrix{T}(H_kernel[internal_indices, internal_indices])

    if internal_dim == 0
        return IntrinsicHessianVerification{T}(
            :rolling_rigid, kernel_dim, 0, contact_pivots, rigid_rows,
            H_kernel, H_internal, T[], nothing
        )
    end

    if _prove_interval_positive_definite(H_internal)
        return IntrinsicHessianVerification{T}(
            :positive_definite, kernel_dim, internal_dim, contact_pivots,
            rigid_rows, H_kernel, H_internal, T[], nothing
        )
    end

    coefficient_witness, negative_value = _find_interval_negative_witness(H_internal)
    if negative_value !== nothing
        kernel_coefficients = zeros(T, kernel_dim)
        kernel_coefficients[internal_indices] .= coefficient_witness
        physical_direction = K * kernel_coefficients
        return IntrinsicHessianVerification{T}(
            :negative_direction, kernel_dim, internal_dim, contact_pivots,
            rigid_rows, H_kernel, H_internal, physical_direction, negative_value
        )
    end

    return IntrinsicHessianVerification{T}(
        :inconclusive_or_degenerate, kernel_dim, internal_dim, contact_pivots,
        rigid_rows, H_kernel, H_internal, T[], nothing
    )
end

function _multiplier_sign_status(lambda_box::Vector{<:ArbReal})
    any(_strictly_negative, lambda_box) && return :negative_present
    all(_strictly_positive, lambda_box) && return :strictly_positive
    all(value -> lowerbound(value) >= 0, lambda_box) && return :nonnegative_proved
    return :sign_undetermined
end

function _leaf_obstruction(c::Configuration{T}, g::ContactGraph) where {T<:ArbReal}
    for vertex in 1:g.n
        degree(g, vertex) == 1 || continue
        edge = only(filter(e -> e.src == vertex || e.dst == vertex, g.edges))
        neighbour = edge.src == vertex ? edge.dst : edge.src
        px = c.centers[1, vertex]
        py = c.centers[2, vertex]
        ax = c.centers[1, vertex] - c.centers[1, neighbour]
        ay = c.centers[2, vertex] - c.centers[2, neighbour]
        determinant_box = px * ay - py * ax
        dot_box = px * ax + py * ay
        if !_contains_zero(determinant_box)
            return :excluded_nonparallel, vertex, determinant_box
        elseif upperbound(dot_box) < 0
            return :excluded_wrong_orientation, vertex, dot_box
        end
    end
    return :not_excluded, nothing, nothing
end

function _verify_contact_geometry(c::Configuration{T}, g::ContactGraph) where {T<:ArbReal}
    edge_set = Set(_edge_tuple.(g.edges))
    minimum_gap = nothing
    unresolved = false
    overlap = false
    for i in 1:(g.n - 1), j in (i + 1):g.n
        (i, j) in edge_set && continue
        gap = squared_distance(c, i, j) - T(4)
        IntervalSafety.finite_interval(gap) ||
            return :nonedge_strictness_unresolved, nothing
        if minimum_gap === nothing || lowerbound(gap) < lowerbound(minimum_gap)
            minimum_gap = gap
        end
        if upperbound(gap) < 0
            overlap = true
        elseif lowerbound(gap) <= 0
            unresolved = true
        end
    end
    overlap && return :nonedge_overlap_proved, minimum_gap
    unresolved && return :nonedge_strictness_unresolved, minimum_gap
    return :exact_graph_proved, minimum_gap
end

function _verify_gauge_transversality(
    c_box::Configuration{T},
    root_left::Int,
    root_right::Int,
) where {T<:ArbReal}
    orientation = c_box.centers[1, root_right] - c_box.centers[1, root_left]
    return lowerbound(orientation) > 0 ?
        :original_balance_proved : :gauge_transversality_unresolved
end

function _target_value(::Type{T}, n::Int) where {T<:ArbReal}
    return IntervalSafety.enclosing_rational(T, second_moment_target(n))
end

function _objective_comparison(value::T, target::T) where {T<:ArbReal}
    lowerbound(value) > upperbound(target) && return :strictly_above_target
    upperbound(value) < lowerbound(target) && return :strictly_below_target
    return :overlaps_target
end

"""
Analyze an angular root already isolated and transferred by `Refinement`.
The Cartesian enclosure is the image of that same root, not the result of a
second Cartesian Krawczyk solve.
"""
function verify_refined_candidate(refined, candidate_id::String; precision::Int=512)
    T = ArbReal{precision}
    g = ContactGraph(refined.graph6)
    n, m = g.n, num_contacts(g)
    c_box = refined.center_box
    lambda_box = refined.lambda_box
    notes = String[
        "The unique angular root was mapped directly to this centered Cartesian enclosure.",
        "The full contact-multiplier enclosure was computed from that same root image.",
    ]
    graph_status, minimum_gap = _verify_contact_geometry(c_box, g)
    root_left, root_right = first(refined.chart.tree)
    gauge_status = _verify_gauge_transversality(c_box, root_left, root_right)
    A_box = contact_matrix(c_box, g)
    valid_minor, _ = _verified_determinant(A_box, refined.contact_pivots)
    rank_status = valid_minor ? :full_row_rank_proved : :full_row_rank_unresolved
    rank_A = valid_minor ? m : nothing
    leaf_status, leaf_vertex, leaf_witness = _leaf_obstruction(c_box, g)
    multiplier_status = valid_minor && gauge_status === :original_balance_proved ?
        _multiplier_sign_status(lambda_box) : :not_checked

    hessian_result = nothing
    if valid_minor && multiplier_status !== :negative_present &&
       multiplier_status !== :not_checked
        hessian_result = _verify_intrinsic_hessian(
            SecondMomentObjective(), c_box, g, lambda_box, refined.contact_pivots,
        )
    end
    objective_value = refined.moment_box
    objective_status = _objective_comparison(objective_value, _target_value(T, n))

    conclusion = :survives_general_checks
    if graph_status === :nonedge_overlap_proved
        conclusion = :input_error
    elseif graph_status !== :exact_graph_proved ||
           gauge_status !== :original_balance_proved || !valid_minor
        conclusion = :special_treatment_required
    elseif leaf_status in (:excluded_nonparallel, :excluded_wrong_orientation)
        conclusion = :excluded_by_leaf_condition
    elseif multiplier_status === :negative_present
        conclusion = :excluded_by_first_order
    elseif hessian_result !== nothing && hessian_result.status === :negative_direction
        conclusion = :excluded_by_negative_intrinsic_direction
    elseif objective_status === :strictly_above_target
        conclusion = :excluded_by_objective_bound
    elseif objective_status === :strictly_below_target
        conclusion = :counterexample_to_announced_bound
    elseif objective_status === :overlaps_target
        conclusion = :target_equality_requires_exact_argument
    end

    return RegularCandidateResult{T}(
        candidate_id, "Search candidate $candidate_id", refined.graph6, n, m,
        precision, :angle_unique_transferred,
        graph_status, minimum_gap, rank_status, rank_A, gauge_status,
        leaf_status, leaf_vertex, leaf_witness,
        multiplier_status, lambda_box, hessian_result, objective_value,
        objective_status, conclusion, notes,
    )
end

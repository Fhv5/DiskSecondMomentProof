module IntervalSearch

export SearchOptions,
       SearchStats,
       AngularChart,
       IsolatedCandidate,
       GraphSearchResult,
       write_search_reports

using ArbNumerics
using LinearAlgebra
using ...DiskSecondMomentProof
using ..LeafFilters
using ..BoundaryRoots

const LE = LeafFilters
const BR = BoundaryRoots
const IS = DiskSecondMomentProof.IntervalSafety
const RationalBig = Rational{BigInt}

"""Precision and operational limits for the interval search."""
struct SearchOptions
    precision::Int
    max_boxes::Int
    max_depth::Int
    minimum_width::RationalBig
    krawczyk_depth::Int
end

SearchOptions(; precision=512, max_boxes=3_000_000, max_depth=140,
              minimum_width_power=60, krawczyk_depth=12) = begin
    precision > 0 || throw(ArgumentError("precision must be positive"))
    max_boxes > 0 || throw(ArgumentError("max_boxes must be positive"))
    max_depth >= 0 || throw(ArgumentError("max_depth cannot be negative"))
    minimum_width_power >= 0 || throw(ArgumentError("minimum_width_power cannot be negative"))
    krawczyk_depth >= 0 || throw(ArgumentError("krawczyk_depth cannot be negative"))
    SearchOptions(precision, max_boxes, max_depth,
                  RationalBig(1, BigInt(2)^minimum_width_power), krawczyk_depth)
end

"""Counters recording all terminal decisions in an interval search."""
mutable struct SearchStats
    visited::Int
    split::Int
    contact::Int
    overlap::Int
    rank::Int
    multiplier::Int
    stationarity::Int
    moment_balance::Int
    krawczyk::Int
    boundary_root::Int
    leaf_equilibrium::Int
    bridge_equilibrium::Int
    bridge_rotation::Int
    symmetry::Int
    bridge_krawczyk::Int
    isolated::Int
    unresolved::Int
end

SearchStats() = SearchStats(ntuple(_ -> 0, fieldcount(SearchStats))...)

"""Tree coordinates and symmetry restrictions used throughout one graph search."""
struct AngularChart
    selected_leaf::Union{Nothing,Int}
    tree::Vector{Tuple{Int,Int}}
    reflection_index::Int
    symmetry_constraints::Vector{Tuple{Int,Int}}
end

function angular_chart(graph; leaf=isempty(graph.leaves) ? nothing : first(graph.leaves))
    if leaf === nothing
        return AngularChart(nothing, rooted_tree(graph), 1, Tuple{Int,Int}[])
    end
    leaf in graph.leaves || throw(ArgumentError("selected anchor is not a graph leaf"))
    tree = LE.anchored_tree(graph, leaf)
    reflection, constraints = LE.angular_symmetry_data(graph, tree)
    return AngularChart(leaf, tree, reflection, constraints)
end

initial_angle_box(graph, chart::AngularChart) = LE.initial_box(graph.n, chart.reflection_index)

"""A unique interval-enclosed critical realization found in one angle box."""
struct IsolatedCandidate
    isolation_box::Vector{ArbReal}
    root_state_box::Vector{ArbReal}
    angle_box::Vector{ArbReal}
    center_box::Matrix{ArbReal}
    moment::ArbReal
    kind::Symbol
    chart::AngularChart
end

"""Complete or unresolved interval-search result for one contact graph."""
struct GraphSearchResult
    graph6::String
    n::Int
    contacts::Int
    expected_dimension::Int
    status::Symbol
    stats::SearchStats
    candidates::Vector{IsolatedCandidate}
    unresolved_boxes::Vector{LE.AngleBox}
    elapsed_seconds::Float64
    chart::AngularChart
end

function rooted_tree(graph)
    root_left, root_right = first(graph.edges)
    tree = Tuple{Int,Int}[(root_left, root_right)]
    seen = Set((root_left, root_right))
    queue = [root_left, root_right]
    cursor = 1
    while cursor <= length(queue)
        parent = queue[cursor]
        cursor += 1
        for child in sort(graph.adjacency[parent])
            child in seen && continue
            push!(seen, child)
            push!(queue, child)
            push!(tree, (parent, child))
        end
    end
    length(tree) == graph.n - 1 || error("graph $(graph.code) is disconnected")
    return tree
end

function initial_angle_box(n::Int)
    pi_upper = RationalBig(22, 7)
    lower = fill(-pi_upper, n - 2)
    upper = fill(pi_upper, n - 2)
    !isempty(lower) && (lower[1] = 0)
    return LE.AngleBox(lower, upper, 0)
end

function interval_angles(::Type{T}, box::LE.AngleBox) where {T<:ArbReal}
    return [IS.enclosing_rational_interval(T, box.lower[i], box.upper[i])
            for i in eachindex(box.lower)]
end

function tree_kinematics(::Type{T}, graph, tree, angles) where {T<:ArbReal}
    k = length(angles)
    raw = fill(T(0), 2, graph.n)
    first_derivative = fill(T(0), 2, graph.n, k)
    second_diagonal = fill(T(0), 2, graph.n, k)
    raw[1, tree[1][2]] = T(2)

    for (angle_index, (parent, child)) in enumerate(tree[2:end])
        angle = angles[angle_index]
        raw[1, child] = raw[1, parent] + T(2) * cos(angle)
        raw[2, child] = raw[2, parent] + T(2) * sin(angle)
        first_derivative[:, child, :] .= first_derivative[:, parent, :]
        second_diagonal[:, child, :] .= second_diagonal[:, parent, :]
        first_derivative[1, child, angle_index] += -T(2) * sin(angle)
        first_derivative[2, child, angle_index] +=  T(2) * cos(angle)
        second_diagonal[1, child, angle_index] += -T(2) * cos(angle)
        second_diagonal[2, child, angle_index] += -T(2) * sin(angle)
    end

    centers = copy(raw)
    centered_first = copy(first_derivative)
    centered_second = copy(second_diagonal)
    for coordinate in 1:2
        centers[coordinate, :] .-= sum(raw[coordinate, :]) / T(graph.n)
        for angle_index in 1:k
            centered_first[coordinate, :, angle_index] .-=
                sum(first_derivative[coordinate, :, angle_index]) / T(graph.n)
            centered_second[coordinate, :, angle_index] .-=
                sum(second_diagonal[coordinate, :, angle_index]) / T(graph.n)
        end
    end
    return raw, centers, first_derivative, centered_first, second_diagonal, centered_second
end

function reduced_terms(::Type{T}, graph, tree, angles) where {T<:ArbReal}
    raw, centers, raw_first, center_first, raw_second, center_second =
        tree_kinematics(T, graph, tree, angles)
    k = length(angles)
    tree_edges = Set((min(a, b), max(a, b)) for (a, b) in tree)
    extra_edges = [edge for edge in graph.edges if !(edge in tree_edges)]
    q = length(extra_edges)

    constraints = fill(T(0), q)
    constraint_jacobian = fill(T(0), q, k)
    constraint_hessians = fill(T(0), q, k, k)
    for (row, (left, right)) in enumerate(extra_edges)
        dx = raw[1, left] - raw[1, right]
        dy = raw[2, left] - raw[2, right]
        constraints[row] = (dx * dx + dy * dy - T(4)) / T(2)
        for a in 1:k
            dax = raw_first[1, left, a] - raw_first[1, right, a]
            day = raw_first[2, left, a] - raw_first[2, right, a]
            constraint_jacobian[row, a] = dx * dax + dy * day
            for b in 1:k
                dbx = raw_first[1, left, b] - raw_first[1, right, b]
                dby = raw_first[2, left, b] - raw_first[2, right, b]
                value = dax * dbx + day * dby
                if a == b
                    ddx = raw_second[1, left, a] - raw_second[1, right, a]
                    ddy = raw_second[2, left, a] - raw_second[2, right, a]
                    value += dx * ddx + dy * ddy
                end
                constraint_hessians[row, a, b] = value
            end
        end
    end

    moment = sum(value * value for value in centers) / T(2)
    gradient = fill(T(0), k)
    hessian = fill(T(0), k, k)
    for a in 1:k
        gradient[a] = sum(
            centers[coordinate, vertex] * center_first[coordinate, vertex, a]
            for coordinate in 1:2, vertex in 1:graph.n
        )
        for b in 1:k
            value = sum(
                center_first[coordinate, vertex, a] * center_first[coordinate, vertex, b]
                for coordinate in 1:2, vertex in 1:graph.n
            )
            if a == b
                value += sum(
                    centers[coordinate, vertex] * center_second[coordinate, vertex, a]
                    for coordinate in 1:2, vertex in 1:graph.n
                )
            end
            hessian[a, b] = value
        end
    end

    return (; raw, centers, moment, extra_edges, constraints,
            constraint_jacobian, constraint_hessians, gradient, hessian)
end

function reduced_kkt(::Type{T}, graph, tree, angles, multipliers) where {T<:ArbReal}
    terms = reduced_terms(T, graph, tree, angles)
    q, k = size(terms.constraint_jacobian)
    residual = vcat(
        terms.constraints,
        terms.gradient - transpose(terms.constraint_jacobian) * multipliers,
    )
    jacobian = fill(T(0), q + k, k + q)
    jacobian[1:q, 1:k] .= terms.constraint_jacobian
    lagrangian_hessian = copy(terms.hessian)
    for edge_index in 1:q
        lagrangian_hessian .-= multipliers[edge_index] .* terms.constraint_hessians[edge_index, :, :]
    end
    jacobian[(q + 1):(q + k), 1:k] .= lagrangian_hessian
    jacobian[(q + 1):(q + k), (k + 1):(k + q)] .= -transpose(terms.constraint_jacobian)
    return residual, jacobian
end

function contact_system(::Type{T}, graph, tree, angles) where {T<:ArbReal}
    terms = reduced_terms(T, graph, tree, angles)
    return terms.constraints, terms.constraint_jacobian
end

contains_zero(value::ArbReal) = !IS.zero_excluded(value)
strictly_outside_zero(value::ArbReal) = !contains_zero(value)

function interval_krawczyk(::Type{T}, evaluator, state_box) where {T<:ArbReal}
    residual_box, jacobian_box = evaluator(state_box)
    all(IS.finite_interval, residual_box) &&
        all(IS.finite_interval, jacobian_box) || return (:unavailable, state_box)
    any(strictly_outside_zero, residual_box) && return (:excluded, state_box)
    numeric_jacobian = LE.midpoint_float_matrix(jacobian_box)
    all(isfinite, numeric_jacobian) || return (:unavailable, state_box)
    size(numeric_jacobian, 1) == size(numeric_jacobian, 2) || return (:unavailable, state_box)
    inverse = try
        inv(numeric_jacobian)
    catch
        nothing
    end
    inverse === nothing && return (:unavailable, state_box)
    all(isfinite, inverse) || return (:unavailable, state_box)
    preconditioner = [T(inverse[row, column]) for row in axes(inverse, 1), column in axes(inverse, 2)]
    midpoint_state = [midpoint(value) for value in state_box]
    residual_midpoint, _ = evaluator(midpoint_state)
    all(IS.finite_interval, residual_midpoint) || return (:unavailable, state_box)
    delta = state_box - midpoint_state
    identity_matrix = Matrix{T}(I, length(state_box), length(state_box))
    image = midpoint_state - preconditioner * residual_midpoint +
            (identity_matrix - preconditioner * jacobian_box) * delta
    all(IS.finite_interval, image) || return (:unavailable, state_box)
    for index in eachindex(state_box)
        if upperbound(image[index]) < lowerbound(state_box[index]) ||
           lowerbound(image[index]) > upperbound(state_box[index])
            return (:excluded, image)
        end
    end
    interior = all(
        lowerbound(image[index]) > lowerbound(state_box[index]) &&
        upperbound(image[index]) < upperbound(state_box[index])
        for index in eachindex(state_box)
    )
    return (interior ? :unique : :inconclusive, image)
end

function candidate_geometry_clear(::Type{T}, graph, tree, angle_box) where {T<:ArbReal}
    terms = reduced_terms(T, graph, tree, angle_box)
    for (left, right) in graph.nonedges
        gap = LE.squared_gap(terms.raw, left, right, T(4))
        lowerbound(gap) > 0 || return false
    end
    return true
end

function increment!(stats::SearchStats, reason::Symbol)
    reason == :contact && (stats.contact += 1)
    reason == :overlap && (stats.overlap += 1)
    reason == :rank && (stats.rank += 1)
    reason == :multiplier && (stats.multiplier += 1)
    reason == :stationarity && (stats.stationarity += 1)
    reason == :moment_balance && (stats.moment_balance += 1)
    reason == :krawczyk && (stats.krawczyk += 1)
    reason == :boundary_root && (stats.boundary_root += 1)
    reason == :leaf_equilibrium && (stats.leaf_equilibrium += 1)
    reason == :bridge_equilibrium && (stats.bridge_equilibrium += 1)
    reason == :bridge_rotation && (stats.bridge_rotation += 1)
    reason == :symmetry && (stats.symmetry += 1)
    reason == :bridge_krawczyk && (stats.bridge_krawczyk += 1)
end

function interval_text(value::ArbReal)
    return "[$(lowerbound(value)), $(upperbound(value))]"
end

const SUMMARY_HEADER = "graph6\tedges\texpected_dimension\tstatus\tvisited\tsplit\tsymmetry\tleaf_equilibrium\tbridge_equilibrium\tbridge_rotation\tbridge_krawczyk\tcontact\toverlap\tmultiplier\tstationarity\tmoment_balance\tkrawczyk\tboundary_root\trank_terminal\tisolated\tunresolved\telapsed_seconds"
const CANDIDATE_HEADER = "graph6\tcandidate\tkind\tmoment_interval\tcenter_intervals"
const BOUNDARY_ROOT_HEADER = "graph6\ttarget_graph6\toccurrences\tangle_pi_multiples\treduced_multipliers\textra_contacts\tpermutation\texact_moment\texact_kkt_determinant"

"""Write complete interval-search reports from one finished run."""
function write_search_reports(results, n::Int, output_dir::String)
    mkpath(output_dir)
    summary_path = joinpath(output_dir, "n$n-summary.tsv")
    candidates_path = joinpath(output_dir, "n$n-candidates.tsv")
    boundary_roots_path = joinpath(output_dir, "n$n-boundary-roots.tsv")
    open(summary_path, "w") do io
        println(io, SUMMARY_HEADER)
        for result in results
            println(io, join((result.graph6, result.contacts, result.expected_dimension,
                              result.status, result.stats.visited, result.stats.split,
                              result.stats.symmetry, result.stats.leaf_equilibrium,
                              result.stats.bridge_equilibrium, result.stats.bridge_rotation,
                              result.stats.bridge_krawczyk, result.stats.contact,
                              result.stats.overlap,
                              result.stats.multiplier, result.stats.stationarity,
                              result.stats.moment_balance, result.stats.krawczyk,
                              result.stats.boundary_root, result.stats.rank,
                              length(result.candidates),
                              length(result.unresolved_boxes), result.elapsed_seconds), '\t'))
        end
    end
    open(candidates_path, "w") do io
        println(io, CANDIDATE_HEADER)
        for result in results
            for (index, candidate) in enumerate(result.candidates)
                centers = join([
                    "$(vertex):($(interval_text(candidate.center_box[1, vertex])),$(interval_text(candidate.center_box[2, vertex])))"
                    for vertex in 1:result.n
                ], ";")
                println(io, join((result.graph6, index, candidate.kind,
                                  interval_text(candidate.moment), centers), '\t'))
            end
        end
    end
    open(boundary_roots_path, "w") do io
        println(io, BOUNDARY_ROOT_HEADER)
        for result in results
            result.stats.boundary_root == 0 && continue
            graph = LE.graph_data(result.graph6)
            tree = result.chart.tree
            for witness in BR.registered_witnesses(result.graph6)
                verified = BR.verify_witness(graph, tree, witness)
                angles = join(BR.rational_text.(witness.angle_pi_multiples), ",")
                multipliers = join(BR.qsqrt3_text.(witness.reduced_multipliers), ",")
                extra_contacts = join(["$left-$right" for (left, right) in witness.extra_contacts], ",")
                permutation = join(witness.permutation, ",")
                println(io, join((
                    result.graph6,
                    witness.target_graph6,
                    result.stats.boundary_root,
                    angles,
                    multipliers,
                    extra_contacts,
                    permutation,
                    BR.qsqrt3_text(witness.expected_moment),
                    BR.qsqrt3_text(verified.determinant),
                ), '\t'))
            end
        end
    end
    return summary_path, candidates_path, boundary_roots_path
end

end

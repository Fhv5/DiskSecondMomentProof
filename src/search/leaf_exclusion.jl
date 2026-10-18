module LeafExclusion

using ArbNumerics
using LinearAlgebra
using ...DiskSecondMomentProof: ContactGraph
import ...DiskSecondMomentProof: IntervalSafety

export LeafAngleBox,
       LeafExclusionOptions,
       LeafExclusionStatistics,
       LeafExclusionRun,
       LeafGraphExclusion,
       LeafCensusExclusion,
       write_leaf_exclusion_reports

"""
Interval branch-and-prune for the contact-leaf obstruction.

For every connected graph with a leaf, the procedure covers all of its labelled
realizations (up to translations, rotations, and reflection) by boxes of
spanning-tree angles.  A box is removed only after an ArbReal interval
calculation proves that at least one necessary condition for an ambient local
minimum is impossible throughout that box.

A completed run with status `proved_excluded` is a proof.  A run that reaches
an operational limit is reported as `unresolved`; it is never silently treated
as an exclusion.
This file contains library code only.  Argument parsing, terminal output, and
filesystem policy belong to the command-line scripts.
"""

const RationalBig = Rational{BigInt}
const PACKAGE_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct LeafGraphData
    code::String
    n::Int
    edges::Vector{Tuple{Int,Int}}
    adjacency::Vector{Vector{Int}}
    nonedges::Vector{Tuple{Int,Int}}
    leaves::Vector{Int}
    bridge_cuts::Vector{Tuple{Int,Vector{Int},Int,Int}}
end

"""Rational angular box used by the exhaustive branch-and-prune cover."""
struct LeafAngleBox
    lower::Vector{RationalBig}
    upper::Vector{RationalBig}
    depth::Int
end

"""Counters recording how boxes were handled during leaf exclusion."""
mutable struct LeafExclusionStatistics
    visited::Int
    split::Int
    excluded_contact::Int
    excluded_overlap::Int
    excluded_leaf_equilibrium::Int
    excluded_bridge_equilibrium::Int
    excluded_bridge_rotation::Int
    excluded_symmetry::Int
    excluded_krawczyk::Int
    excluded_multiplier::Int
    excluded_stationarity::Int
    excluded_moment_balance::Int
    unresolved_rank::Int
    unresolved_limit::Int
end

LeafExclusionStatistics() = LeafExclusionStatistics(
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

"""Result of applying the obstruction through one selected leaf."""
struct LeafExclusionRun
    leaf::Int
    status::Symbol
    stats::LeafExclusionStatistics
    elapsed_seconds::Float64
    unresolved::Vector{LeafAngleBox}
end

"""Result for one graph in the contact-leaf exclusion."""
struct LeafGraphExclusion
    graph6::String
    n::Int
    leaves::Vector{Int}
    status::Symbol
    runs::Vector{LeafExclusionRun}
end

"""Result for a complete supplied census of connected contact graphs."""
struct LeafCensusExclusion
    n::Int
    total_graphs::Int
    leaf_graphs::Int
    proved_graphs::Int
    status::Symbol
    graph_results::Vector{LeafGraphExclusion}
    elapsed_seconds::Float64
end

"""Precision and operational limits for interval leaf exclusion."""
struct LeafExclusionOptions
    precision::Int
    max_boxes::Int
    max_depth::Int
    minimum_width::RationalBig
    all_leaves::Bool
end

LeafExclusionOptions(; precision::Int=512,
                     max_boxes::Int=250_000,
                     max_depth::Int=80,
                     minimum_width_power::Int=24,
                     all_leaves::Bool=false) = begin
    precision > 0 || throw(ArgumentError("precision must be positive"))
    max_boxes > 0 || throw(ArgumentError("max_boxes must be positive"))
    max_depth >= 0 || throw(ArgumentError("max_depth cannot be negative"))
    minimum_width_power >= 0 ||
        throw(ArgumentError("minimum_width_power cannot be negative"))
    LeafExclusionOptions(
        precision,
        max_boxes,
        max_depth,
        RationalBig(1, BigInt(2)^minimum_width_power),
        all_leaves,
    )
end

function graph_data(code::AbstractString)
    code = String(code)
    contact_graph = ContactGraph(code)
    n = contact_graph.n
    edges = sort([(edge.src, edge.dst) for edge in contact_graph.edges],
                 by=edge -> (edge[2], edge[1]))
    neighbors = [Int[] for _ in 1:n]
    edge_set = Set(edges)
    for (left, right) in edges
        push!(neighbors[left], right)
        push!(neighbors[right], left)
    end
    nonedges = [(i, j) for i in 1:(n - 1) for j in (i + 1):n if !((i, j) in edge_set)]
    leaves = [i for i in 1:n if length(neighbors[i]) == 1]
    bridge_cuts = Tuple{Int,Vector{Int},Int,Int}[]
    for (edge_index, (left, right)) in enumerate(edges)
        reached = Set([left])
        queue = [left]
        cursor = 1
        while cursor <= length(queue)
            vertex = queue[cursor]
            cursor += 1
            for next in neighbors[vertex]
                ((vertex == left && next == right) || (vertex == right && next == left)) && continue
                next in reached && continue
                push!(reached, next)
                push!(queue, next)
            end
        end
        right in reached && continue
        side = sort!(collect(reached))
        inside, outside = left, right
        if length(side) > n - length(side)
            side = [vertex for vertex in 1:n if !(vertex in reached)]
            inside, outside = right, left
        end
        push!(bridge_cuts, (edge_index, side, inside, outside))
    end
    return LeafGraphData(code, n, edges, neighbors, nonedges, leaves, bridge_cuts)
end

function anchored_tree(graph::LeafGraphData, leaf::Int)
    neighbor = only(graph.adjacency[leaf])
    tree = Tuple{Int,Int}[(leaf, neighbor)]
    seen = Set((leaf, neighbor))
    queue = [leaf, neighbor]
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
    length(tree) == graph.n - 1 || error("graph $(graph.code) is not connected")
    return tree
end

function initial_box(n::Int, reflection_index::Int)
    pi_upper = RationalBig(22, 7)
    dimensions = n - 2
    lower = fill(-pi_upper, dimensions)
    upper = fill(pi_upper, dimensions)
    dimensions > 0 && (lower[reflection_index] = 0)
    return LeafAngleBox(lower, upper, 0)
end

function swapped_vertex(vertex::Int, first::Int, second::Int)
    vertex == first && return second
    vertex == second && return first
    return vertex
end

function vertex_swap_is_automorphism(graph::LeafGraphData, first::Int, second::Int)
    original = Set(graph.edges)
    mapped = Set{Tuple{Int,Int}}()
    for (left, right) in graph.edges
        new_left = swapped_vertex(left, first, second)
        new_right = swapped_vertex(right, first, second)
        push!(mapped, (min(new_left, new_right), max(new_left, new_right)))
    end
    return mapped == original
end

const VERIFIED_SYMMETRY_PATTERNS = Dict{String,Symbol}(
    "D@{" => :P0, "DB{" => :P0, "D\x60[" => :P0,
    "E?Fg" => :P0, "E?Nw" => :P0, "E?dg" => :P0,
    "E?lw" => :P0, "E@ow" => :P0, "E@po" => :P0,
    "E@pw" => :P0, "EGEw" => :P0, "E_Lw" => :P0,
    "E?No" => :P1, "E@NW" => :P1, "ECXo" => :P1,
    "E\x60LW" => :P1,
    "DBw" => :P2, "DJk" => :P2, "E?lo" => :P2,
    "E@lo" => :P2,
    "D?{" => :P3, "E?Bw" => :P3, "E?Fw" => :P3,
    "EGFw" => :P3,
)

function tree_descendant(vertex::Int, ancestor::Int, parent_of::Dict{Int,Int})
    current = vertex
    while true
        current == ancestor && return true
        haskey(parent_of, current) || return false
        current = parent_of[current]
    end
end

function verified_symmetry_pattern(graph::LeafGraphData, tree, groups,
                                   reflection_index::Int, parent, child)
    pattern = get(VERIFIED_SYMMETRY_PATTERNS, graph.code, :unsupported)
    pattern == :unsupported && return false
    grouped = Set(index for group in groups for index in group)
    parent_of = Dict(last(edge) => first(edge) for edge in tree)
    for group in groups, first_index in group, second_index in group
        parent[first_index] == parent[second_index] || return false
        vertex_swap_is_automorphism(
            graph, child[first_index], child[second_index],
        ) || return false
    end
    is_tree_leaf(index) = !(child[index] in parent)

    if pattern == :P0
        return !(reflection_index in grouped) &&
               all(is_tree_leaf, grouped)
    elseif pattern == :P1
        reflection_index in grouped && return false
        for group in groups, index in group
            tree_descendant(child[reflection_index], child[index], parent_of) &&
                return false
            for other_group in groups
                other_group === group && continue
                any(other -> tree_descendant(
                    child[other], child[index], parent_of,
                ), other_group) && return false
            end
        end
        return true
    elseif pattern == :P2
        length(groups) == 1 || return false
        length(only(groups)) == 2 || return false
        reflection_index in grouped && return false
        first_index, second_index = only(groups)
        first_child, second_child = child[first_index], child[second_index]
        origin = parent[first_index]
        reflection_parent = parent[reflection_index]
        tip = child[reflection_index]
        return reflection_parent == first_child &&
               length(Set((origin, first_child, second_child, tip))) == 4 &&
               tip in graph.adjacency[first_child] &&
               tip in graph.adjacency[second_child]
    elseif pattern == :P3
        return length(grouped) == length(child) &&
               reflection_index in grouped &&
               all(is_tree_leaf, grouped)
    end
    return false
end

"""Return only symmetry reductions whose tree pattern is proved for this census."""
function angular_symmetry_data(graph::LeafGraphData, tree)
    dimension = graph.n - 2
    parent = [tree[index + 1][1] for index in 1:dimension]
    child = [tree[index + 1][2] for index in 1:dimension]
    interchangeable = falses(dimension, dimension)
    for first in 1:(dimension - 1), second in (first + 1):dimension
        if parent[first] == parent[second] &&
           vertex_swap_is_automorphism(graph, child[first], child[second])
            interchangeable[first, second] = true
            interchangeable[second, first] = true
        end
    end
    groups = Vector{Vector{Int}}()
    seen = falses(dimension)
    for start in 1:dimension
        seen[start] && continue
        group = [index for index in 1:dimension if index == start || interchangeable[start, index]]
        for index in group
            seen[index] = true
        end
        length(group) > 1 && push!(groups, group)
    end

    symmetric_indices = Set{Int}()
    for group in groups
        union!(symmetric_indices, group)
    end
    reflection_index = something(findfirst(index -> !(index in symmetric_indices), 1:dimension), 1)
    constraints = Tuple{Int,Int}[]
    for group in groups
        ordered = filter(!=(reflection_index), group)
        if length(ordered) >= 2
            for index in 1:(length(ordered) - 1)
                push!(constraints, (ordered[index], ordered[index + 1]))
            end
        end
    end
    if !isempty(constraints) &&
       !verified_symmetry_pattern(graph, tree, groups, reflection_index,
                                  parent, child)
        return reflection_index, Tuple{Int,Int}[]
    end
    return reflection_index, constraints
end

function symmetry_possible(box::LeafAngleBox, constraints)
    for (first, second) in constraints
        box.lower[first] > box.upper[second] && return false
    end
    return true
end

function interval_from_rationals(::Type{T}, lo::RationalBig, hi::RationalBig) where {T<:ArbReal}
    return IntervalSafety.enclosing_rational_interval(T, lo, hi)
end

function interval_angles(::Type{T}, box::LeafAngleBox) where {T<:ArbReal}
    return [interval_from_rationals(T, box.lower[i], box.upper[i]) for i in eachindex(box.lower)]
end

function centers_from_angles(::Type{T}, graph::LeafGraphData, tree, angles) where {T<:ArbReal}
    raw = fill(T(0), 2, graph.n)
    placed = falses(graph.n)
    leaf, neighbor = tree[1]
    placed[leaf] = true
    placed[neighbor] = true
    raw[1, neighbor] = T(2)
    for (angle, (parent, child)) in zip(angles, tree[2:end])
        placed[parent] || error("invalid rooted tree order")
        placed[child] && error("tree child was already placed")
        raw[1, child] = raw[1, parent] + T(2) * cos(angle)
        raw[2, child] = raw[2, parent] + T(2) * sin(angle)
        placed[child] = true
    end
    centers = copy(raw)
    for coordinate in 1:2
        mean_coordinate = sum(raw[coordinate, :]) / T(graph.n)
        centers[coordinate, :] .-= mean_coordinate
    end
    return raw, centers
end

function centers_with_angle_derivatives(::Type{T}, graph::LeafGraphData, tree, angles) where {T<:ArbReal}
    dimension = length(angles)
    raw = fill(T(0), 2, graph.n)
    derivative = fill(T(0), 2, graph.n, dimension)
    leaf, neighbor = tree[1]
    raw[1, neighbor] = T(2)
    for (angle_index, (parent, child)) in enumerate(tree[2:end])
        raw[1, child] = raw[1, parent] + T(2) * cos(angles[angle_index])
        raw[2, child] = raw[2, parent] + T(2) * sin(angles[angle_index])
        derivative[:, child, :] .= derivative[:, parent, :]
        derivative[1, child, angle_index] += -T(2) * sin(angles[angle_index])
        derivative[2, child, angle_index] += T(2) * cos(angles[angle_index])
    end
    return raw, derivative
end

function contact_matrix_interval(::Type{T}, centers, graph::LeafGraphData) where {T<:ArbReal}
    matrix = fill(T(0), length(graph.edges), 2 * graph.n)
    for (row, (left, right)) in enumerate(graph.edges)
        dx = centers[1, left] - centers[1, right]
        dy = centers[2, left] - centers[2, right]
        matrix[row, 2left - 1] = dx
        matrix[row, 2left] = dy
        matrix[row, 2right - 1] = -dx
        matrix[row, 2right] = -dy
    end
    return matrix
end

contains_zero(value::ArbReal) =
    !IntervalSafety.zero_excluded(value)
strictly_outside_zero(value::ArbReal) = !contains_zero(value)

function squared_gap(centers, left::Int, right::Int, four)
    dx = centers[1, left] - centers[1, right]
    dy = centers[2, left] - centers[2, right]
    return dx * dx + dy * dy - four
end

function leaf_equilibrium_possible(centers, raw_centers, graph::LeafGraphData, leaf_limit)
    for leaf in graph.leaves
        neighbor = only(graph.adjacency[leaf])
        ux = raw_centers[1, leaf] - raw_centers[1, neighbor]
        uy = raw_centers[2, leaf] - raw_centers[2, neighbor]
        px = centers[1, leaf]
        py = centers[2, leaf]
        cross = px * uy - py * ux
        contains_zero(cross) || return false
        coefficient = (px * ux + py * uy) / 4
        upperbound(coefficient) < 0 && return false
        lowerbound(coefficient) > upperbound(leaf_limit) && return false
    end
    return true
end

function bridge_rotation_possible(centers, graph::LeafGraphData)
    all_vertices = 1:graph.n
    for (_, side, _, outside) in graph.bridge_cuts
        complement = [vertex for vertex in all_vertices if !(vertex in side)]
        side_x = sum(centers[1, vertex] for vertex in side) / length(side)
        side_y = sum(centers[2, vertex] for vertex in side) / length(side)
        other_x = sum(centers[1, vertex] for vertex in complement) / length(complement)
        other_y = sum(centers[2, vertex] for vertex in complement) / length(complement)
        vx = side_x - centers[1, outside]
        vy = side_y - centers[2, outside]
        wx = centers[1, outside] - other_x
        wy = centers[2, outside] - other_y
        first_derivative = -wx * vy + wy * vx
        contains_zero(first_derivative) || return false
        second_derivative_without_positive_factor = -(wx * vx + wy * vy)
        upperbound(second_derivative_without_positive_factor) < 0 && return false
    end
    return true
end

function bridge_rotation_system(centers, derivative, graph::LeafGraphData)
    dimension = size(derivative, 3)
    residuals = Vector{eltype(centers)}(undef, length(graph.bridge_cuts))
    jacobian = Matrix{eltype(centers)}(undef, length(graph.bridge_cuts), dimension)
    all_vertices = 1:graph.n
    for (row, (_, side, _, outside)) in enumerate(graph.bridge_cuts)
        complement = [vertex for vertex in all_vertices if !(vertex in side)]
        sx = sum(centers[1, vertex] for vertex in side) / length(side)
        sy = sum(centers[2, vertex] for vertex in side) / length(side)
        tx = sum(centers[1, vertex] for vertex in complement) / length(complement)
        ty = sum(centers[2, vertex] for vertex in complement) / length(complement)
        vx, vy = sx - centers[1, outside], sy - centers[2, outside]
        wx, wy = centers[1, outside] - tx, centers[2, outside] - ty
        residuals[row] = -wx * vy + wy * vx
        for angle_index in 1:dimension
            dsx = sum(derivative[1, vertex, angle_index] for vertex in side) / length(side)
            dsy = sum(derivative[2, vertex, angle_index] for vertex in side) / length(side)
            dtx = sum(derivative[1, vertex, angle_index] for vertex in complement) / length(complement)
            dty = sum(derivative[2, vertex, angle_index] for vertex in complement) / length(complement)
            dvx = dsx - derivative[1, outside, angle_index]
            dvy = dsy - derivative[2, outside, angle_index]
            dwx = derivative[1, outside, angle_index] - dtx
            dwy = derivative[2, outside, angle_index] - dty
            jacobian[row, angle_index] =
                -dwx * vy - wx * dvy + dwy * vx + wy * dvx
        end
    end
    return residuals, jacobian
end

function krawczyk_excludes_bridge_root(::Type{T}, graph::LeafGraphData, tree,
                                       angle_box) where {T<:ArbReal}
    dimension = length(angle_box)
    length(graph.bridge_cuts) >= dimension || return false
    raw_box, derivative_box = centers_with_angle_derivatives(T, graph, tree, angle_box)
    _, jacobian_box = bridge_rotation_system(raw_box, derivative_box, graph)
    all(IntervalSafety.finite_interval, jacobian_box) || return false
    numeric_jacobian = midpoint_float_matrix(jacobian_box)
    all(isfinite, numeric_jacobian) || return false
    row_factorization = qr(transpose(numeric_jacobian), ColumnNorm())
    rows = Vector(row_factorization.p)[1:dimension]
    square_numeric = numeric_jacobian[rows, :]
    R_numeric = try
        inv(square_numeric)
    catch
        nothing
    end
    R_numeric === nothing && return false
    all(isfinite, R_numeric) || return false
    R = [T(R_numeric[row, column]) for row in 1:dimension, column in 1:dimension]

    center_angles = [midpoint(value) for value in angle_box]
    raw_mid, derivative_mid = centers_with_angle_derivatives(T, graph, tree, center_angles)
    residual_mid, _ = bridge_rotation_system(raw_mid, derivative_mid, graph)
    all(IntervalSafety.finite_interval, residual_mid) || return false
    selected_residual = residual_mid[rows]
    selected_jacobian = jacobian_box[rows, :]
    delta = angle_box - center_angles
    identity_matrix = Matrix{T}(I, dimension, dimension)
    image = center_angles - R * selected_residual +
            (identity_matrix - R * selected_jacobian) * delta
    all(IntervalSafety.finite_interval, image) || return false
    for index in 1:dimension
        if upperbound(image[index]) < lowerbound(angle_box[index]) ||
           lowerbound(image[index]) > upperbound(angle_box[index])
            return true
        end
    end
    return false
end

function midpoint_float_matrix(matrix)
    return [Float64(midpoint(matrix[row, column]))
            for row in axes(matrix, 1), column in axes(matrix, 2)]
end

function plain_transpose(matrix::AbstractMatrix{T}) where {T}
    result = Matrix{T}(undef, size(matrix, 2), size(matrix, 1))
    for row in axes(matrix, 1), column in axes(matrix, 2)
        result[column, row] = matrix[row, column]
    end
    return result
end

function widest_dimension(box::LeafAngleBox)
    widths = box.upper .- box.lower
    _, index = findmax(widths)
    return index, widths[index]
end

function split_box(box::LeafAngleBox, dimension::Int)
    middle = (box.lower[dimension] + box.upper[dimension]) / 2
    left_upper = copy(box.upper)
    left_upper[dimension] = middle
    right_lower = copy(box.lower)
    right_lower[dimension] = middle
    return (
        LeafAngleBox(copy(box.lower), left_upper, box.depth + 1),
        LeafAngleBox(right_lower, copy(box.upper), box.depth + 1),
    )
end

function increment_exclusion!(stats::LeafExclusionStatistics, reason::Symbol)
    reason == :contact && (stats.excluded_contact += 1)
    reason == :overlap && (stats.excluded_overlap += 1)
    reason == :leaf_equilibrium && (stats.excluded_leaf_equilibrium += 1)
    reason == :bridge_equilibrium && (stats.excluded_bridge_equilibrium += 1)
    reason == :bridge_rotation && (stats.excluded_bridge_rotation += 1)
    reason == :symmetry && (stats.excluded_symmetry += 1)
    reason == :krawczyk && (stats.excluded_krawczyk += 1)
    reason == :multiplier && (stats.excluded_multiplier += 1)
    reason == :stationarity && (stats.excluded_stationarity += 1)
    reason == :moment_balance && (stats.excluded_moment_balance += 1)
end

function box_text(box::LeafAngleBox)
    intervals = ["$(box.lower[i]):$(box.upper[i])" for i in eachindex(box.lower)]
    return join(intervals, ";")
end

"""Write leaf-filter decisions and unresolved interval boxes."""
function write_leaf_exclusion_reports(result::LeafCensusExclusion,
                                      output_directory::AbstractString)
    mkpath(output_directory)
    summary_path = joinpath(output_directory,
                            "n$(result.n)-leaf-interval-summary.tsv")
    unresolved_path = joinpath(output_directory,
                               "n$(result.n)-leaf-interval-unresolved.tsv")
    open(summary_path, "w") do io
        println(io, "graph6\tleaf\tstatus\tvisited\tsplit\tcontact\toverlap\tleaf_equilibrium\tbridge_equilibrium\tbridge_rotation\tsymmetry\tkrawczyk\tmultiplier\tstationarity\tmoment_balance\tunresolved_rank\tunresolved_limit\tseconds")
        for graph_result in result.graph_results, run in graph_result.runs
            s = run.stats
            println(io, join((graph_result.graph6, run.leaf, run.status, s.visited, s.split,
                              s.excluded_contact, s.excluded_overlap,
                              s.excluded_leaf_equilibrium,
                              s.excluded_bridge_equilibrium,
                              s.excluded_bridge_rotation,
                              s.excluded_symmetry,
                              s.excluded_krawczyk,
                              s.excluded_multiplier, s.excluded_stationarity,
                              s.excluded_moment_balance, s.unresolved_rank,
                              s.unresolved_limit,
                              round(run.elapsed_seconds; digits=6)), '\t'))
        end
    end
    open(unresolved_path, "w") do io
        println(io, "graph6\tleaf\tdepth\tangle_intervals_exact")
        for graph_result in result.graph_results, run in graph_result.runs, box in run.unresolved
            println(io, "$(graph_result.graph6)\t$(run.leaf)\t$(box.depth)\t$(box_text(box))")
        end
    end
    return summary_path, unresolved_path
end

end

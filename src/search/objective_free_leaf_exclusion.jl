module ObjectiveFreeLeafExclusion

using ArbNumerics
using ..DiskSecondMomentProof

const LE = DiskSecondMomentProof.CandidateSearch.LeafExclusion

export prove_leaf_census

"""Test necessary equilibrium and leaf-force bounds across contact bridges."""
function bridge_equilibrium_possible(centers, raw_centers, graph, leaf_limit)
    for (edge_index, side, inside, outside) in graph.bridge_cuts
        sx = sum(centers[1, vertex] for vertex in side)
        sy = sum(centers[2, vertex] for vertex in side)
        ux = raw_centers[1, inside] - raw_centers[1, outside]
        uy = raw_centers[2, inside] - raw_centers[2, outside]
        LE.contains_zero(sx * uy - sy * ux) || return false
        coefficient = (sx * ux + sy * uy) / 4
        upperbound(coefficient) < 0 && return false
        left, right = graph.edges[edge_index]
        if left in graph.leaves || right in graph.leaves
            lowerbound(coefficient) > upperbound(leaf_limit) && return false
        end
    end
    return true
end

"""Intersect multiplier enclosures with nonnegative and leaf-specific bounds."""
function nonnegative_multipliers(::Type{T}, enclosure, graph, leaf_limit) where {T<:ArbReal}
    all(DiskSecondMomentProof.IntervalSafety.finite_interval, enclosure) ||
        return :unresolved
    bounded = Vector{T}(undef, length(enclosure))
    for (index, value) in enumerate(enclosure)
        left, right = graph.edges[index]
        lower = DiskSecondMomentProof.IntervalSafety.interval_max(lowerbound(value), zero(T))
        upper = upperbound(value)
        if left in graph.leaves || right in graph.leaves
            upper = DiskSecondMomentProof.IntervalSafety.interval_min(upper, upperbound(leaf_limit))
        end
        lower <= upper || return nothing
        bounded[index] = setinterval(lower, upper)
    end
    return bounded
end

"""Apply necessary local-minimum conditions without an objective target bound."""
function classify_box(
    ::Type{T},
    graph,
    tree,
    symmetry_constraints,
    box,
    leaf_limit,
) where {T<:ArbReal}
    LE.symmetry_possible(box, symmetry_constraints) || return :symmetry
    angles = LE.interval_angles(T, box)
    raw_centers, centers = LE.centers_from_angles(T, graph, tree, angles)
    tree_set = Set((min(a, b), max(a, b)) for (a, b) in tree)

    for (left, right) in graph.edges
        (left, right) in tree_set && continue
        LE.strictly_outside_zero(LE.squared_gap(raw_centers, left, right, T(4))) &&
            return :contact
    end
    for (left, right) in graph.nonedges
        upperbound(LE.squared_gap(raw_centers, left, right, T(4))) <= 0 &&
            return :overlap
    end

    LE.leaf_equilibrium_possible(centers, raw_centers, graph, leaf_limit) ||
        return :leaf_equilibrium
    bridge_equilibrium_possible(centers, raw_centers, graph, leaf_limit) ||
        return :bridge_equilibrium
    LE.bridge_rotation_possible(centers, graph) || return :bridge_rotation
    if box.depth >= 18 && LE.krawczyk_excludes_bridge_root(T, graph, tree, angles)
        return :krawczyk
    end

    matrix = LE.contact_matrix_interval(T, raw_centers, graph)
    objective_gradient = vec(centers)
    solve_status, enclosure = DiskSecondMomentProof.FullEquilibrium.enclose_multipliers(
        T, matrix, objective_gradient)
    solve_status == :stationarity && return :stationarity
    solve_status == :rank && return :subdivide_rank
    multipliers = nonnegative_multipliers(T, enclosure, graph, leaf_limit)
    multipliers === :unresolved && return :subdivide_rank
    multipliers === nothing && return :multiplier
    residual = LE.plain_transpose(matrix) * multipliers - objective_gradient
    any(LE.strictly_outside_zero, residual) && return :stationarity

    moment = sum(value * value for value in centers) / T(2)
    LE.strictly_outside_zero(T(2) * sum(multipliers) - moment) &&
        return :moment_balance
    return :subdivide_possible
end

"""Cover the anchored angle domain for one selected leaf."""
function prove_for_leaf(::Type{T}, graph, leaf, options; trace=nothing) where {T<:ArbReal}
    tree = LE.anchored_tree(graph, leaf)
    reflection_index, symmetry_constraints = LE.angular_symmetry_data(graph, tree)
    leaf_limit = DiskSecondMomentProof.IntervalSafety.enclosing_rational(
        T, (graph.n - 1) // graph.n)
    stack = LE.LeafAngleBox[LE.initial_box(graph.n, reflection_index)]
    unresolved = LE.LeafAngleBox[]
    stats = LE.LeafExclusionStatistics()
    started = time()

    while !isempty(stack) && stats.visited < options.max_boxes
        box = pop!(stack)
        stats.visited += 1
        reason = classify_box(
            T,
            graph,
            tree,
            symmetry_constraints,
            box,
            leaf_limit,
        )
        if reason in (
            :contact,
            :overlap,
            :leaf_equilibrium,
            :bridge_equilibrium,
            :bridge_rotation,
            :symmetry,
            :krawczyk,
            :multiplier,
            :stationarity,
            :moment_balance,
        )
            trace !== nothing && trace(box, reason)
            LE.increment_exclusion!(stats, reason)
            continue
        end

        dimension, width = LE.widest_dimension(box)
        if box.depth >= options.max_depth || width <= options.minimum_width
            trace !== nothing && trace(box, :unresolved)
            reason == :subdivide_rank && (stats.unresolved_rank += 1)
            push!(unresolved, box)
            continue
        end
        first_child, second_child = LE.split_box(box, dimension)
        trace !== nothing && trace(box, :split)
        push!(stack, second_child)
        push!(stack, first_child)
        stats.split += 1
    end
    if !isempty(stack)
        stats.unresolved_limit += length(stack)
        append!(unresolved, stack)
    end
    status = isempty(unresolved) ? :proved_excluded : :unresolved
    return LE.LeafExclusionRun(leaf, status, stats, time() - started, unresolved)
end

"""Attempt a complete exclusion for one graph containing a contact leaf."""
function prove_leaf_graph(graph, options; trace_directory=nothing)
    isempty(graph.leaves) && return LE.LeafGraphExclusion(
        graph.code,
        graph.n,
        Int[],
        :not_applicable,
        LE.LeafExclusionRun[],
    )
    T = Core.apply_type(ArbReal, options.precision)
    runs = LE.LeafExclusionRun[]
    leaves = options.all_leaves ? graph.leaves : graph.leaves[1:1]
    for leaf in leaves
        run = if trace_directory === nothing
            prove_for_leaf(T, graph, leaf, options)
        else
            path = DiskSecondMomentProof.BoxRecords.box_record_path(
                trace_directory, "leaf", graph.code; leaf=leaf)
            mkpath(dirname(path))
            open(path, "w") do io
                prove_for_leaf(T, graph, leaf, options;
                    trace=(box, reason) ->
                        DiskSecondMomentProof.BoxRecords.record_box(io, box, reason))
            end
        end
        push!(runs, run)
        run.status == :proved_excluded && break
    end
    status = any(run -> run.status == :proved_excluded, runs) ?
        :proved_excluded : :unresolved
    return LE.LeafGraphExclusion(
        graph.code,
        graph.n,
        copy(graph.leaves),
        status,
        runs,
    )
end

"""Exclude leaf graphs, reserving only the two-disk target; retain unresolved boxes."""
function prove_leaf_census(codes; options=LE.LeafExclusionOptions(),
                           verbose=false, trace_directory=nothing, reserved_exact_targets=String[])
    graphs = LE.graph_data.(String.(codes))
    isempty(graphs) && throw(ArgumentError("the census cannot be empty"))
    n = first(graphs).n
    n in (2, 3, 4, 5, 6) || throw(ArgumentError("the proof supports two, three, four, five or six disks"))
    all(graph -> graph.n == n, graphs) || throw(ArgumentError("mixed disk counts"))
    reserved = Set(String.(reserved_exact_targets))
    if n == 2
        length(graphs) == 1 && only(graphs).code == "A_" && reserved == Set(["A_"]) ||
            throw(ArgumentError("the two-disk contact graph must be reserved for algebraic target analysis"))
    else
        isempty(reserved) || throw(ArgumentError("only the two-disk target may bypass leaf exclusion"))
    end
    leaf_graphs = filter(graph -> !isempty(graph.leaves) && !(graph.code in reserved), graphs)
    started = time()
    results = LE.LeafGraphExclusion[]
    for (index, graph) in enumerate(leaf_graphs)
        verbose && println("[$index/$(length(leaf_graphs))] $(graph.code)")
        push!(results, prove_leaf_graph(graph, options;
            trace_directory=trace_directory))
    end
    proved = count(result -> result.status == :proved_excluded, results)
    status = proved == length(results) ? :proved_excluded : :unresolved
    return LE.LeafCensusExclusion(
        first(graphs).n,
        length(graphs),
        length(leaf_graphs),
        proved,
        status,
        results,
        time() - started,
    )
end

end

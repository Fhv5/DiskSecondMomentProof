module MinimumCandidateSearch

using ArbNumerics
using ..DiskSecondMomentProof
const LS = DiskSecondMomentProof.CandidateSearch.LeaflessSearch
const LE = DiskSecondMomentProof.CandidateSearch.LeafExclusion
const BR = DiskSecondMomentProof.CandidateSearch.BoundaryRoots

"""Enclose contact multipliers and reject incompatible stationarity or signs."""
function nonnegative_full_multipliers(::Type{T}, graph, terms) where {T<:ArbReal}
    matrix = LE.contact_matrix_interval(T, terms.raw, graph)
    objective_gradient = vec(terms.centers)
    solve_status, enclosure = DiskSecondMomentProof.FullEquilibrium.enclose_multipliers(
        T, matrix, objective_gradient)
    solve_status == :stationarity && return (:stationarity, nothing)
    solve_status == :rank && return (:rank, nothing)

    multipliers = T[]
    for value in enclosure
        DiskSecondMomentProof.IntervalSafety.finite_interval(value) ||
            return (:rank, nothing)
        upperbound(value) < 0 && return (:multiplier, nothing)
        lower = DiskSecondMomentProof.IntervalSafety.interval_max(lowerbound(value), zero(T))
        upper = upperbound(value)
        lower <= upper || return (:multiplier, nothing)
        push!(multipliers, setinterval(lower, upper))
    end

    residual = LE.plain_transpose(matrix) * multipliers - objective_gradient
    any(LS.strictly_outside_zero, residual) && return (:stationarity, nothing)
    return (:possible, multipliers)
end

"""Contract an already isolated root enclosure using interval Krawczyk iterations."""
function refine_unique_image(::Type{T}, evaluator, image; iterations=8) where {T<:ArbReal}
    refined = image
    for _ in 1:iterations
        status, next_image = LS.interval_krawczyk(T, evaluator, refined)
        status == :unique || break
        refined = next_image
    end
    return refined
end

"""Apply objective-free geometric and first-order filters to an angle box."""
function classify_candidate_box(
    ::Type{T},
    graph,
    tree,
    boundary_roots,
    box,
    options,
) where {T<:ArbReal}
    angles = LS.interval_angles(T, box)
    terms = LS.reduced_terms(T, graph, tree, angles)
    any(LS.strictly_outside_zero, terms.constraints) && return (:contact, nothing)
    for (left, right) in graph.nonedges
        upperbound(LE.squared_gap(terms.raw, left, right, T(4))) <= 0 &&
            return (:overlap, nothing)
    end

    multiplier_status, full_multipliers =
        nonnegative_full_multipliers(T, graph, terms)
    multiplier_status == :rank && return (:rank, nothing)
    multiplier_status == :stationarity && return (:stationarity, nothing)
    multiplier_status == :multiplier && return (:multiplier, nothing)
    LS.strictly_outside_zero(T(2) * sum(full_multipliers) - terms.moment) &&
        return (:moment_balance, nothing)

    box.depth < options.krawczyk_depth && return (:subdivide, nothing)
    q = length(terms.extra_edges)
    k = length(angles)
    if q == k
        evaluator = state -> LS.contact_system(T, graph, tree, state)
        status, image = LS.interval_krawczyk(T, evaluator, angles)
        status == :excluded && return (:krawczyk, nothing)
        if status == :unique
            image = refine_unique_image(T, evaluator, image)
            any(root -> BR.witness_state_in_box(T, root, angles, ArbReal[]),
                boundary_roots) && return (:boundary_root, nothing)
            LS.candidate_geometry_clear(T, graph, tree, image) ||
                return (:subdivide, nothing)
            candidate_terms = LS.reduced_terms(T, graph, tree, image)
            return (:isolated, LS.IsolatedCandidate(
                angles,
                image,
                image,
                candidate_terms.centers,
                candidate_terms.moment,
                :rigid_contact_root,
            ))
        end
    elseif q < k
        extra_indices = [findfirst(==(edge), graph.edges) for edge in terms.extra_edges]
        reduced_multiplier_box = [full_multipliers[index] for index in extra_indices]
        state_box = vcat(angles, reduced_multiplier_box)
        evaluator = state -> begin
            state_angles = state[1:k]
            state_multipliers = state[(k + 1):(k + q)]
            LS.reduced_kkt(T, graph, tree, state_angles, state_multipliers)
        end
        status, image = LS.interval_krawczyk(T, evaluator, state_box)
        status == :excluded && return (:krawczyk, nothing)
        if status == :unique
            image = refine_unique_image(T, evaluator, image)
            any(root -> BR.witness_state_in_box(
                    T, root, angles, reduced_multiplier_box), boundary_roots) &&
                return (:boundary_root, nothing)
            candidate_angles = image[1:k]
            LS.candidate_geometry_clear(T, graph, tree, candidate_angles) ||
                return (:subdivide, nothing)
            candidate_terms = LS.reduced_terms(T, graph, tree, candidate_angles)
            return (:isolated, LS.IsolatedCandidate(
                state_box,
                image,
                candidate_angles,
                candidate_terms.centers,
                candidate_terms.moment,
                :flexible_critical_root,
            ))
        end
    else
        return (:rank, nothing)
    end
    return (:subdivide, nothing)
end

"""Cover one graph's angle domain and retain isolated roots and unresolved boxes."""
function search_candidate_graph(graph6::AbstractString, options::LS.SearchOptions;
                                trace=nothing)
    graph6 = String(graph6)
    graph = LE.graph_data(graph6)
    graph.n in (3, 4, 5, 6) || throw(ArgumentError("the proof supports three, four, five or six disks"))
    isempty(graph.leaves) || throw(ArgumentError("candidate isolation requires a leafless graph"))
    tree = LS.rooted_tree(graph)
    boundary_roots = [
        BR.verify_witness(graph, tree, witness)
        for witness in BR.registered_witnesses(graph6)
    ]
    stack = LE.LeafAngleBox[LS.initial_angle_box(graph.n)]
    unresolved = LE.LeafAngleBox[]
    candidates = LS.IsolatedCandidate[]
    stats = LS.SearchStats()
    started = time()
    T = ArbReal{options.precision}

    while !isempty(stack) && stats.visited < options.max_boxes
        box = pop!(stack)
        stats.visited += 1
        reason, candidate = classify_candidate_box(
            T,
            graph,
            tree,
            boundary_roots,
            box,
            options,
        )
        if reason == :isolated
            trace !== nothing && trace(box, :isolated)
            push!(candidates, candidate)
            stats.isolated += 1
            continue
        elseif reason != :subdivide && reason != :rank
            trace !== nothing && trace(box, reason)
            LS.increment!(stats, reason)
            continue
        end

        dimension, width = LE.widest_dimension(box)
        if box.depth >= options.max_depth || width <= options.minimum_width
            trace !== nothing && trace(box, :unresolved)
            reason == :rank && (stats.rank += 1)
            push!(unresolved, box)
            stats.unresolved += 1
            continue
        end
        first_child, second_child = LE.split_box(box, dimension)
        trace !== nothing && trace(box, :split)
        push!(stack, second_child)
        push!(stack, first_child)
        stats.split += 1
    end

    if !isempty(stack)
        append!(unresolved, stack)
        stats.unresolved += length(stack)
    end
    status = isempty(unresolved) ? :complete : :unresolved
    return LS.GraphSearchResult(
        graph6,
        graph.n,
        length(graph.edges),
        2 * graph.n - 3 - length(graph.edges),
        status,
        stats,
        candidates,
        unresolved,
        time() - started,
    )
end


end

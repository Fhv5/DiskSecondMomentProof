module Refinement

using ArbNumerics
using LinearAlgebra
using JSON3
using ..DiskSecondMomentProof: Configuration, ContactGraph, contact_matrix,
    num_contacts, gradient, SecondMomentObjective, evaluate
import ..DiskSecondMomentProof

const LS = DiskSecondMomentProof.CandidateSearch.LeaflessSearch
const LE = DiskSecondMomentProof.CandidateSearch.LeafExclusion

export RefinedCandidate, refine_candidate, write_refinement_report

"""One isolated angular root transferred to centered Cartesian interval boxes."""
struct RefinedCandidate{T<:ArbReal}
    graph6::String
    kind::Symbol
    isolation_box::Vector{T}
    root_state_box::Vector{T}
    angle_box::Vector{T}
    center_box::Configuration{T}
    lambda_box::Vector{T}
    contact_pivots::Vector{Int}
    moment_box::T
end

function _interval_record(value::ArbReal)
    return DiskSecondMomentProof.IntervalSafety.interval_record(value)
end

function _box_contained(inner, outer)
    length(inner) == length(outer) || return false
    return all(
        lowerbound(outer[i]) <= lowerbound(inner[i]) &&
        upperbound(inner[i]) <= upperbound(outer[i])
        for i in eachindex(inner)
    )
end

"""Write interval endpoints for each transferred root without decimal seeds."""
function write_refinement_report(records, path::AbstractString)
    mkpath(dirname(path))
    open(path, "w") do io
        for (identifier, candidate) in records
            JSON3.write(io, Dict(
                "candidate_id" => identifier,
                "graph6" => candidate.graph6,
                "kind" => string(candidate.kind),
                "angular_state" => _interval_record.(candidate.root_state_box),
                "centers" => [
                    [_interval_record(candidate.center_box.centers[axis, vertex]) for axis in 1:2]
                    for vertex in 1:size(candidate.center_box.centers, 2)
                ],
                "multipliers" => _interval_record.(candidate.lambda_box),
                "moment" => _interval_record(candidate.moment_box),
                "contact_pivots" => candidate.contact_pivots,
            ))
            println(io)
        end
    end
    return path
end

function _root_evaluator(::Type{T}, graph, tree, kind::Symbol, k::Int, q::Int) where {T<:ArbReal}
    if kind === :rigid_contact_root && q == k
        return state -> LS.contact_system(T, graph, tree, state)
    elseif kind === :flexible_critical_root && q < k
        return state -> LS.reduced_kkt(T, graph, tree, state[1:k], state[(k + 1):(k + q)])
    end
    throw(ArgumentError("candidate kind does not match the reduced system"))
end

"""
Refine the same unique angular root found by the search and map it to Cartesian
centers and complete contact multipliers. No independent Cartesian root search
or decimal-root matching is performed.
"""
function refine_candidate(graph6::AbstractString, candidate::LS.IsolatedCandidate;
                          precision::Int=512, iterations::Int=16)
    T = ArbReal{precision}
    graph = LE.graph_data(graph6)
    tree = LS.rooted_tree(graph)
    k = graph.n - 2
    q = length(graph.edges) - (graph.n - 1)
    evaluator = _root_evaluator(T, graph, tree, candidate.kind, k, q)
    source_box = T.(candidate.isolation_box)
    state = T.(candidate.root_state_box)
    expected_length = candidate.kind === :rigid_contact_root ? k : k + q
    length(source_box) == expected_length && length(state) == expected_length ||
        throw(ArgumentError("candidate state has the wrong dimension"))

    status, _ = LS.interval_krawczyk(T, evaluator, source_box)
    status === :unique || throw(ArgumentError("angular root isolation could not be reverified"))
    _box_contained(state, source_box) ||
        throw(ArgumentError("refined angular root is outside its isolation box"))

    for _ in 1:iterations
        next_status, next_state = LS.interval_krawczyk(T, evaluator, state)
        next_status === :unique || break
        state = next_state
    end

    angles = Vector{T}(state[1:k])
    terms = LS.reduced_terms(T, graph, tree, angles)
    center_box = Configuration(Matrix{T}(terms.centers))
    contact_graph = ContactGraph(String(graph6))
    A = contact_matrix(center_box, contact_graph)
    m = num_contacts(contact_graph)
    pivots, _ = DiskSecondMomentProof._first_verified_column_minor(A, m)
    pivots === nothing && throw(ArgumentError("contact rank remains unresolved after refinement"))

    gradient_box = gradient(SecondMomentObjective(), center_box)
    square = DiskSecondMomentProof._plain_transpose(A[:, pivots])
    verified = DiskSecondMomentProof.IntervalSafety.verified_inverse(square)
    verified === nothing && throw(ArgumentError("multiplier system has no verified interval inverse"))
    lambda_box = Vector{T}(first(verified) * gradient_box[pivots])
    residual = gradient_box - DiskSecondMomentProof._plain_transpose(A) * lambda_box
    all(DiskSecondMomentProof._contains_zero, residual) ||
        throw(ArgumentError("transferred multiplier enclosure is incompatible with stationarity"))

    return RefinedCandidate(
        String(graph6), candidate.kind, source_box, state, angles,
        center_box, lambda_box, pivots, evaluate(SecondMomentObjective(), center_box),
    )
end

end

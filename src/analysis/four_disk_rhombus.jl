"""Evidence linking a five-contact endpoint to the whole rhombus family."""
struct FourDiskEndpointResult
    candidate_id::String
    graph6::String
    n::Int
    family_id::String
    contact_rank::Int
    minimum_noncontact_gap::Qsqrt3
    multipliers::Vector{Qsqrt3}
    objective::Qsqrt3
    target::Qsqrt3
    boundary_to_endpoint_permutations::Dict{String,Vector{Int}}
    conclusion::Symbol
end

const _RHOMBUS_EDGES = Set(((1, 2), (1, 3), (2, 4), (3, 4)))

function _four_disk_rhombus_terms(entry::Dict{String,Any})
    get(entry, "representation", "") == "exact_family" ||
        throw(ArgumentError("rhombus input must use the `exact_family` representation"))
    graph = ContactGraph(string(entry["graph6"]))
    graph.n == 4 && Set(_edge_tuple.(graph.edges)) == _RHOMBUS_EDGES ||
        throw(ArgumentError("rhombus coverage requires the declared four-cycle labeling"))
    ContactGraph(4, _entry_edges(entry)).edges == graph.edges ||
        throw(ArgumentError("rhombus contact list mismatch"))
    family = entry["family"]::Dict{String,Any}
    get(family, "model", "") == "four_disk_rhombus" ||
        throw(ArgumentError("rhombus coverage requires the four_disk_rhombus model"))
    domain = family["domain"]::Dict{String,Any}
    get(domain, "relation", "") == "b^2=4-a^2" &&
        get(domain, "positive_axes", false) === true ||
        throw(ArgumentError("rhombus coverage requires positive a,b and b^2=4-a^2"))
    interval = domain["a_squared_interval"]::Dict{String,Any}
    _parse_rational(string(interval["lower"])) == 1//1 &&
        _parse_rational(string(interval["upper"])) == 3//1 &&
        interval["lower_closed"] === false && interval["upper_closed"] === false ||
        throw(ArgumentError("rhombus interior must be exactly 1<a^2<3"))
    raw = family["centers"]::Dict{String,Any}
    Set(keys(raw)) == Set(string.(1:4)) ||
        throw(ArgumentError("rhombus center labels must be exactly 1:4"))
    centers = _parse_arc_centers(raw, 4)
    a = _arc(_RationalBig[0, 1], _RationalBig[0])
    b = _arc(_RationalBig[0], _RationalBig[1])
    z = _arc(0)
    expected = [(-a, z), (z, b), (z, -b), (a, z)]
    all(centers[row, i] == expected[i][row] for i in 1:4, row in 1:2) ||
        throw(ArgumentError("rhombus centers do not equal the exhaustive family"))
    return graph, family, centers, a, b
end

"""
Bind the input to the circle-intersection argument in
docs/proof.md, including both endpoint graphs.
"""
function verify_four_disk_rhombus_coverage(entry::Dict{String,Any})
    _four_disk_rhombus_terms(entry)
    return :complete_penny_family
end

function _verify_four_disk_rhombus_family(entry, family, graph)
    graph, family, centers, a, b = _four_disk_rhombus_terms(entry)
    all(iszero, vec(sum(centers, dims=2))) || error("rhombus centroid identity failed")
    for edge in graph.edges
        distance = sum((centers[row, edge.src] - centers[row, edge.dst])^2 for row in 1:2)
        distance == _arc(4) || error("rhombus contact identity failed")
    end
    # Both gaps are strictly positive exactly on the declared squared interval.
    gaps = Dict((1, 4) => _arc(4) * a^2 - _arc(4),
                (2, 3) => _arc(12) - _arc(4) * a^2)
    for (i, j) in _all_nonedges(graph)
        gap = sum((centers[row, i] - centers[row, j])^2 for row in 1:2) - _arc(4)
        gap == gaps[(i, j)] || error("rhombus noncontact identity failed")
    end
    multipliers = _parse_rational_multiplier_map(family["multipliers"], graph)
    all(==(1//2), multipliers) || throw(ArgumentError("rhombus multipliers must equal 1/2"))
    matrix = _arc_contact_matrix(centers, graph)
    columns = Int.(family["rank_minor_columns"])
    columns == [1, 2, 3, 5] || throw(ArgumentError("rhombus rank minor columns changed"))
    _ring_determinant(Matrix(matrix[:, columns])) == _arc(-2) * a^3 * b ||
        error("rhombus rank determinant identity failed")
    # a,b>0 on the whole domain, so this determinant never vanishes.
    reaction = zeros(_ArcScalar, 8)
    for row in axes(matrix, 1), column in axes(matrix, 2)
        reaction[column] += matrix[row, column] * multipliers[row]
    end
    reaction == vec(centers) || error("rhombus stationarity identity failed")
    direction = _parse_arc_centers(family["flat_direction"], 4)
    all(iszero, matrix * vec(direction)) &&
        all(iszero, vec(sum(direction, dims=2))) ||
        throw(ArgumentError("rhombus flat direction is not centered and tangent"))
    rotation_pairing = sum(-centers[2, i] * direction[1, i] +
                            centers[1, i] * direction[2, i] for i in 1:4)
    iszero(rotation_pairing) && sum(value^2 for value in direction) == _arc(8) ||
        throw(ArgumentError("rhombus flat direction must be nonzero and orthogonal to rigid motions"))
    iszero(_arc_second_variation(direction, graph, multipliers)) ||
        throw(ArgumentError("rhombus intrinsic second variation is not identically zero"))
    _arc_objective(centers) == _arc(4) && second_moment_target(4) == 4//1 ||
        error("rhombus objective identity failed")
    return FamilyResult(string(entry["critical_id"]), string(entry["graph6"]), 4,
        :four_disk_rhombus, "1<a^2<3; a>0; b>0; b^2=4-a^2", 4, 2,
        multipliers, 4//1, second_moment_target(4), true, true, true,
        :exact_squared_parameter_positivity, true, :full_row_rank_proved, columns,
        :exact_flat_intrinsic_direction, nothing, :target_family_exactly_verified)
end

function _rhombus_evaluate(value::_ArcScalar, a::Qsqrt3, b::Qsqrt3)
    function polynomial(coefficients)
        result = zero(Qsqrt3)
        for coefficient in reverse(coefficients)
            result = result * a + Qsqrt3(coefficient)
        end
        return result
    end
    return polynomial(value.p) + b * polynomial(value.q)
end

function _four_disk_boundary_centers(family_entry, boundary::Int)
    _, _, symbolic, _, _ = _four_disk_rhombus_terms(family_entry)
    a, b = boundary == 1 ? (Qsqrt3(1), Qsqrt3(0//1, 1//1)) :
           boundary == 3 ? (Qsqrt3(0//1, 1//1), Qsqrt3(1)) :
           throw(ArgumentError("rhombus boundary must be 1 or 3"))
    return [_rhombus_evaluate(symbolic[row, i], a, b) for row in 1:2, i in 1:4]
end

_four_disk_distance(centers, i, j) = sum((centers[row, i] - centers[row, j])^2 for row in 1:2)

"""Check a rigid contact graph as a boundary of the verified equality family."""
function verify_four_disk_endpoint(entry::Dict{String,Any}, family_entry::Dict{String,Any})
    verify_family(family_entry)
    get(entry, "representation", "") == "exact_qsqrt3" ||
        throw(ArgumentError("rhombus endpoint must use the Q(sqrt(3)) representation"))
    graph = ContactGraph(string(entry["graph6"]))
    graph.n == 4 && length(graph.edges) == 5 ||
        throw(ArgumentError("rhombus endpoint requires four vertices and five contacts"))
    ContactGraph(4, _entry_edges(entry)).edges == graph.edges ||
        throw(ArgumentError("rhombus endpoint contact list mismatch"))
    _validate_qsqrt3_centers(entry["centros"], 4, string(entry["critical_id"]))
    centers = _parse_qsqrt3_centers(entry, 4)
    all(iszero, vec(sum(centers, dims=2))) || throw(ArgumentError("rhombus endpoint is not centered"))
    edges = Set(_edge_tuple.(graph.edges))
    minimum_gap = nothing
    for i in 1:3, j in (i + 1):4
        distance = _four_disk_distance(centers, i, j)
        if (i, j) in edges
            distance == Qsqrt3(4) || throw(ArgumentError("rhombus endpoint contact equation failed"))
        else
            minimum_gap = distance - Qsqrt3(4)
            _q3_sign(minimum_gap) > 0 || throw(ArgumentError("rhombus endpoint noncontact failed"))
        end
    end
    boundary_permutations = Dict{String,Vector{Int}}()
    for boundary in (1, 3)
        representative = _four_disk_boundary_centers(family_entry, boundary)
        for permutation in CandidateSearch.TargetRealizations.all_permutations(4)
            if all(_four_disk_distance(representative, i, j) ==
                   _four_disk_distance(centers, permutation[i], permutation[j])
                   for i in 1:3 for j in (i + 1):4)
                boundary_permutations[string(boundary)] = permutation
                break
            end
        end
        haskey(boundary_permutations, string(boundary)) ||
            throw(ArgumentError("target does not represent both rhombus family boundaries"))
    end
    matched_permutation = boundary_permutations["1"]
    base_edges = Set(minmax(matched_permutation[i], matched_permutation[j]) for (i, j) in _RHOMBUS_EDGES)
    matrix = _qsqrt3_contact_matrix(centers, graph)
    rank = _qsqrt3_rank(matrix)
    rank == 5 || throw(ArgumentError("rhombus endpoint contact rank is not five"))
    multipliers = _parse_edge_multiplier_map(entry["multipliers"], graph)
    expected = Qsqrt3[(edge.src, edge.dst) in base_edges ? Qsqrt3(1//2) : Qsqrt3(0) for edge in graph.edges]
    multipliers == expected && transpose(matrix) * multipliers == vec(centers) ||
        throw(ArgumentError("rhombus endpoint requires four half multipliers and one zero"))
    objective = _qsqrt3_second_moment(centers)
    objective == Qsqrt3(second_moment_target(4)) || error("rhombus endpoint objective failed")
    return FourDiskEndpointResult(string(entry["critical_id"]), string(entry["graph6"]), 4,
        string(family_entry["critical_id"]), rank, minimum_gap, multipliers,
        objective, Qsqrt3(second_moment_target(4)),
        boundary_permutations, :target_endpoint_exactly_verified)
end

"""Require both parts of the four-disk equality set and their linkage."""
function assert_four_disk_closure(analysis)
    endpoints = filter(result -> result.route == :exact_target_endpoint, analysis.results)
    families = filter(result -> result.route == :exact_family &&
        result.evidence isa FamilyResult && result.evidence.model == :four_disk_rhombus,
        analysis.results)
    length(endpoints) == length(families) == 1 ||
        error("Four-disk closure requires exactly one endpoint and one rhombus family")
    endpoint, family = only(endpoints), only(families)
    endpoint.evidence isa FourDiskEndpointResult &&
        endpoint.conclusion == endpoint.evidence.conclusion == :target_endpoint_exactly_verified &&
        family.conclusion == family.evidence.conclusion == :target_family_exactly_verified &&
        endpoint.evidence.family_id == family.candidate_id == family.evidence.candidate_id &&
        endpoint.candidate_id == endpoint.evidence.candidate_id &&
        endpoint.evidence.objective == endpoint.evidence.target == Qsqrt3(4) &&
        family.evidence.objective == family.evidence.target == 4//1 &&
        all((family.evidence.centroid_verified, family.evidence.contacts_verified,
             family.evidence.noncontacts_verified, family.evidence.stationarity_verified)) &&
        family.evidence.rank_status == :full_row_rank_proved &&
        family.evidence.second_order_status == :exact_flat_intrinsic_direction &&
        Set(keys(endpoint.evidence.boundary_to_endpoint_permutations)) == Set(["1", "3"]) &&
        all(sort(permutation) == collect(1:4) for permutation in values(endpoint.evidence.boundary_to_endpoint_permutations)) &&
        endpoint.evidence.contact_rank == 5 &&
        endpoint.evidence.minimum_noncontact_gap == Qsqrt3(8) ||
        error("Four-disk equality evidence is incomplete or inconsistent")
    return Set([endpoint.candidate_id, family.candidate_id])
end

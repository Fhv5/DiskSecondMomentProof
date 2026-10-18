const _RationalBig = Rational{BigInt}

"""
An element p(a) + s*q(a) of Q[a,s]/(s^2 - (4-a^2)).

Polynomial coefficients are stored in ascending powers of a. This polynomial
arithmetic is sufficient for the explicit one-parameter family in the
input data and avoids introducing a computer-algebra dependency.
"""
struct _ArcScalar
    p::Vector{_RationalBig}
    q::Vector{_RationalBig}
end

function _trim_polynomial(coefficients::Vector{_RationalBig})
    result = copy(coefficients)
    while length(result) > 1 && iszero(result[end])
        pop!(result)
    end
    isempty(result) && push!(result, BigInt(0) // BigInt(1))
    return result
end

_arc(p::Vector{_RationalBig}, q::Vector{_RationalBig}) =
    _ArcScalar(_trim_polynomial(p), _trim_polynomial(q))
_arc(value::_RationalBig) = _arc(_RationalBig[value], _RationalBig[0//1])
_arc(value::Integer) = _arc(BigInt(value) // BigInt(1))

function _poly_add(left::Vector{_RationalBig}, right::Vector{_RationalBig})
    result = zeros(_RationalBig, max(length(left), length(right)))
    result[1:length(left)] .+= left
    result[1:length(right)] .+= right
    return _trim_polynomial(result)
end

function _poly_negate(poly::Vector{_RationalBig})
    return _trim_polynomial(-poly)
end

function _poly_multiply(left::Vector{_RationalBig}, right::Vector{_RationalBig})
    result = zeros(_RationalBig, length(left) + length(right) - 1)
    for (i, x) in enumerate(left), (j, y) in enumerate(right)
        result[i + j - 1] += x * y
    end
    return _trim_polynomial(result)
end

function _multiply_by_four_minus_a_squared(poly::Vector{_RationalBig})
    return _poly_multiply(poly, _RationalBig[4//1, 0//1, -1//1])
end

Base.zero(::Type{_ArcScalar}) = _arc(0)
Base.zero(::_ArcScalar) = _arc(0)
Base.one(::Type{_ArcScalar}) = _arc(1)
Base.one(::_ArcScalar) = _arc(1)
Base.iszero(value::_ArcScalar) = all(iszero, value.p) && all(iszero, value.q)
Base.:(==)(left::_ArcScalar, right::_ArcScalar) = left.p == right.p && left.q == right.q
Base.:+(left::_ArcScalar, right::_ArcScalar) =
    _arc(_poly_add(left.p, right.p), _poly_add(left.q, right.q))
Base.:-(value::_ArcScalar) = _arc(_poly_negate(value.p), _poly_negate(value.q))
Base.:-(left::_ArcScalar, right::_ArcScalar) = left + (-right)
function Base.:*(left::_ArcScalar, right::_ArcScalar)
    rational_part = _poly_add(_poly_multiply(left.p, right.p),
                              _multiply_by_four_minus_a_squared(_poly_multiply(left.q, right.q)))
    radical_part = _poly_add(_poly_multiply(left.p, right.q),
                             _poly_multiply(left.q, right.p))
    return _arc(rational_part, radical_part)
end
Base.:*(coefficient::_RationalBig, value::_ArcScalar) = _arc(
    _RationalBig[coefficient * term for term in value.p],
    _RationalBig[coefficient * term for term in value.q])
Base.:*(value::_ArcScalar, coefficient::_RationalBig) = coefficient * value
Base.:/(value::_ArcScalar, divisor::Integer) = (BigInt(1) // BigInt(divisor)) * value
Base.:^(value::_ArcScalar, exponent::Integer) =
    exponent < 0 ? throw(ArgumentError("negative powers are unsupported")) :
    Base.power_by_squaring(value, exponent)

function _parse_polynomial(raw::Vector)
    isempty(raw) && return _RationalBig[0//1]
    return _trim_polynomial(_RationalBig[_parse_rational(string(value)) for value in raw])
end

function _parse_arc_scalar(raw::Dict{String,Any})
    Set(keys(raw)) == Set(("p", "q")) ||
        throw(ArgumentError("an arc scalar must contain exactly p and q"))
    raw["p"] isa Vector && raw["q"] isa Vector ||
        throw(ArgumentError("arc coefficients p and q must be arrays"))
    return _arc(_parse_polynomial(raw["p"]), _parse_polynomial(raw["q"]))
end

function _parse_arc_centers(raw_centers::Dict{String,Any}, n::Int)
    centers = Matrix{_ArcScalar}(undef, 2, n)
    for i in 1:n, coordinate in 1:2
        centers[coordinate, i] = _parse_arc_scalar(raw_centers[string(i)][coordinate])
    end
    return centers
end

function _arc_contact_matrix(centers::Matrix{_ArcScalar}, graph::ContactGraph)
    matrix = zeros(_ArcScalar, num_contacts(graph), 2 * graph.n)
    for (row, edge) in enumerate(graph.edges), coordinate in 1:2
        value = centers[coordinate, edge.src] - centers[coordinate, edge.dst]
        matrix[row, 2 * edge.src - 2 + coordinate] = value
        matrix[row, 2 * edge.dst - 2 + coordinate] = -value
    end
    return matrix
end

function _arc_objective(centers::Matrix{_ArcScalar})
    total = zero(_ArcScalar)
    for value in centers
        total += value^2
    end
    return total / 2
end

function _polynomial_value(poly::Vector{_RationalBig}, x::_RationalBig)
    result = BigInt(0) // BigInt(1)
    for coefficient in reverse(poly)
        result = result * x + coefficient
    end
    return result
end

"""Prove p(a)>0 throughout an open rational interval for degree at most two."""
function _positive_polynomial_on_open_interval(poly::Vector{_RationalBig},
                                               lower::_RationalBig,
                                               upper::_RationalBig)
    lower < upper || return false
    p = _trim_polynomial(poly)
    degree = length(p) - 1
    degree <= 2 || throw(ArgumentError("polynomial positivity currently supports degree at most two"))
    all(iszero, p) && return false

    left_value = _polynomial_value(p, lower)
    right_value = _polynomial_value(p, upper)
    left_value >= 0 && right_value >= 0 || return false

    if degree == 2 && p[3] > 0
        vertex = -p[2] / (2 * p[3])
        if lower < vertex < upper
            _polynomial_value(p, vertex) > 0 || return false
        end
    end
    return true
end

function _parse_rational_multiplier_map(raw::Dict{String,Any}, graph::ContactGraph)
    expected_keys = Set("$(edge.src)-$(edge.dst)" for edge in graph.edges)
    Set(keys(raw)) == expected_keys ||
        throw(ArgumentError("family multiplier keys must match the contact edges"))
    return _RationalBig[_parse_rational(string(raw["$(edge.src)-$(edge.dst)"]))
                          for edge in graph.edges]
end

"""Verification data over the domain of a continuous critical family."""
struct ExactFamilyResult
    candidate_id::String
    graph6::String
    n::Int
    model::Symbol
    parameter_domain::String
    contact_count::Int
    noncontact_count::Int
    multipliers::Vector{_RationalBig}
    objective::_RationalBig
    target::_RationalBig
    centroid_verified::Bool
    contacts_verified::Bool
    noncontacts_verified::Bool
    noncontact_method::Symbol
    stationarity_verified::Bool
    rank_status::Symbol
    rank_witness_columns::Vector{Int}
    second_order_status::Symbol
    negative_second_variation::Union{Nothing,_RationalBig}
    conclusion::Symbol
end

# Evidence files use this concrete type name.
const FamilyResult = ExactFamilyResult

function _all_nonedges(graph::ContactGraph)
    edges = Set(_edge_tuple.(graph.edges))
    return [(i, j) for i in 1:(graph.n - 1) for j in (i + 1):graph.n if !((i, j) in edges)]
end

function _verify_balanced_star_family(entry::Dict{String,Any}, family::Dict{String,Any},
                                      graph::ContactGraph)
    center = Int(family["center_vertex"])
    leaves = Int.(family["leaf_vertices"])
    length(unique(leaves)) == length(leaves) || throw(ArgumentError("leaf vertices must be distinct"))
    Set(vcat(leaves, [center])) == Set(1:graph.n) ||
        throw(ArgumentError("the center and leaves must partition all vertices"))

    domain = family["domain"]::Dict{String,Any}
    _parse_rational(string(domain["leaf_squared_norm"])) == 4//1 ||
        throw(ArgumentError("balanced-star leaves must have squared norm 4"))
    raw_vector_sum = domain["leaf_vector_sum"]
    raw_vector_sum isa Vector && length(raw_vector_sum) == 2 &&
        all(value -> _parse_rational(string(value)) == 0//1, raw_vector_sum) ||
        throw(ArgumentError("balanced-star vectors must sum exactly to zero"))
    domain["strict_graph_nonedges"] === true ||
        throw(ArgumentError("the family domain must impose every graph noncontact strictly"))

    extra_edges = ContactEdge[]
    for raw_edge in get(family, "extra_contacts", Any[])
        push!(extra_edges, ContactEdge(Int(raw_edge[1]), Int(raw_edge[2])))
    end
    all(edge -> edge.src != center && edge.dst != center, extra_edges) ||
        throw(ArgumentError("balanced-star extra contacts must join leaves"))
    expected_graph = ContactGraph(graph.n,
        vcat([ContactEdge(leaf, center) for leaf in leaves], extra_edges))
    expected_graph.edges == graph.edges ||
        throw(ArgumentError("balanced-star description does not generate the declared graph"))

    multipliers = _parse_rational_multiplier_map(family["multipliers"], graph)
    positions = [zeros(_RationalBig, length(leaves)) for _ in 1:graph.n]
    for (index, leaf) in enumerate(leaves)
        positions[leaf][index] = 1//1
    end
    reactions = [zeros(_RationalBig, length(leaves)) for _ in 1:graph.n]
    for (edge, multiplier) in zip(graph.edges, multipliers)
        difference = positions[edge.src] - positions[edge.dst]
        reactions[edge.src] += multiplier * difference
        reactions[edge.dst] -= multiplier * difference
    end
    zero_mod_balance(coefficients) = all(value -> value == coefficients[1], coefficients)
    centroid_verified = zero_mod_balance(reduce(+, positions))
    stationarity_verified = all(vertex -> zero_mod_balance(reactions[vertex] - positions[vertex]),
                                1:graph.n)
    centroid_verified || throw(ArgumentError("balanced-star centroid identity failed"))
    stationarity_verified || throw(ArgumentError("balanced-star stationarity identity failed"))

    objective = BigInt(2 * length(leaves)) // BigInt(1)
    target = second_moment_target(graph.n)
    objective > target || throw(ArgumentError("balanced-star family is not above the target"))
    return FamilyResult(string(entry["critical_id"]), string(entry["graph6"]), graph.n,
        :balanced_star, "semialgebraic domain declared in the input",
        num_contacts(graph), length(_all_nonedges(graph)), multipliers, objective, target,
        centroid_verified, true, true, :encoded_in_exact_domain, stationarity_verified,
        :not_checked, Int[],
        :not_checked, nothing,
        :excluded_by_exact_whole_family_objective)
end

function _ring_determinant(matrix::Matrix{_ArcScalar})
    rows, columns = size(matrix)
    rows == columns || throw(DimensionMismatch("determinant requires a square matrix"))
    rows == 0 && return one(_ArcScalar)
    rows == 1 && return matrix[1, 1]
    value = zero(_ArcScalar)
    for column in 1:columns
        remaining_columns = [index for index in 1:columns if index != column]
        minor = Matrix{_ArcScalar}(matrix[2:end, remaining_columns])
        term = matrix[1, column] * _ring_determinant(minor)
        value += isodd(column) ? term : -term
    end
    return value
end

function _verify_arc_rank_witness(matrix::Matrix{_ArcScalar}, family::Dict{String,Any},
                                  lower::_RationalBig, upper::_RationalBig)
    raw = family["rank_minor"]::Dict{String,Any}
    columns = Int.(raw["columns"])
    length(columns) == size(matrix, 1) && length(unique(columns)) == length(columns) ||
        throw(ArgumentError("arc-family rank minor has the wrong size"))
    all(column -> 1 <= column <= size(matrix, 2), columns) ||
        throw(ArgumentError("arc-family rank minor column is out of range"))
    determinant = _ring_determinant(Matrix{_ArcScalar}(matrix[:, columns]))

    constant = _parse_rational(string(raw["constant"]))
    iszero(constant) && throw(ArgumentError("arc-family rank-minor constant must be nonzero"))
    expected = _arc(constant)
    for raw_factor in raw["positive_factors"]
        factor = _parse_polynomial(raw_factor["coefficients"])
        _positive_polynomial_on_open_interval(factor, lower, upper) ||
            throw(ArgumentError("arc-family rank factor is not positive on the domain"))
        multiplicity = Int(raw_factor["multiplicity"])
        multiplicity > 0 || throw(ArgumentError("rank-factor multiplicity must be positive"))
        expected *= _arc(factor, _RationalBig[0//1])^multiplicity
    end
    determinant == expected ||
        throw(ArgumentError("arc-family rank-minor determinant identity failed"))
    return columns
end

function _arc_second_variation(direction::Matrix{_ArcScalar}, graph::ContactGraph,
                               multipliers::Vector{_RationalBig})
    means = [sum(direction[coordinate, :]) / graph.n for coordinate in 1:2]
    value = zero(_ArcScalar)
    for vertex in 1:graph.n, coordinate in 1:2
        value += (direction[coordinate, vertex] - means[coordinate])^2
    end
    for (edge, multiplier) in zip(graph.edges, multipliers), coordinate in 1:2
        difference = direction[coordinate, edge.src] - direction[coordinate, edge.dst]
        value -= multiplier * difference^2
    end
    return value
end

function _verify_arc_family(entry::Dict{String,Any}, family::Dict{String,Any},
                            graph::ContactGraph)
    domain = family["domain"]::Dict{String,Any}
    string(domain["relation"]) == "s^2=4-a^2" ||
        throw(ArgumentError("unsupported arc relation"))
    interval = domain["a_interval"]::Dict{String,Any}
    interval["lower_closed"] === false && interval["upper_closed"] === false ||
        throw(ArgumentError("the arc parameter interval must be open"))
    lower = _parse_rational(string(interval["lower"]))
    upper = _parse_rational(string(interval["upper"]))
    (-2//1 <= lower < upper <= 2//1) ||
        throw(ArgumentError("the declared interval does not keep s^2=4-a^2 positive"))

    raw_centers = family["centers"]
    raw_centers isa Dict{String,Any} && Set(keys(raw_centers)) == Set(string.(1:graph.n)) ||
        throw(ArgumentError("arc-family center labels must be exactly 1:$(graph.n)"))
    centers = _parse_arc_centers(raw_centers, graph.n)
    centroid_verified = all(iszero, vec(sum(centers, dims=2)))
    centroid_verified || throw(ArgumentError("arc-family centroid identity failed"))

    edge_set = Set(_edge_tuple.(graph.edges))
    contacts_verified = true
    noncontacts_verified = true
    for i in 1:(graph.n - 1), j in (i + 1):graph.n
        distance_squared = sum((centers[row, i] - centers[row, j])^2 for row in 1:2)
        gap = distance_squared - _arc(4)
        if (i, j) in edge_set
            contacts_verified &= iszero(gap)
        else
            iszero(gap.q) || throw(ArgumentError("noncontact gap contains an unresolved radical term"))
            noncontacts_verified &= _positive_polynomial_on_open_interval(gap.p, lower, upper)
        end
    end
    contacts_verified || throw(ArgumentError("arc-family contact identity failed"))
    noncontacts_verified || throw(ArgumentError("arc-family noncontact inequality failed"))

    multipliers = _parse_rational_multiplier_map(family["multipliers"], graph)
    matrix = _arc_contact_matrix(centers, graph)
    rank_witness_columns = _verify_arc_rank_witness(matrix, family, lower, upper)
    gradient_vector = vec(centers)
    reaction = zeros(_ArcScalar, 2 * graph.n)
    for row in axes(matrix, 1), column in axes(matrix, 2)
        reaction[column] += matrix[row, column] * multipliers[row]
    end
    stationarity_verified = reaction == gradient_vector
    stationarity_verified || throw(ArgumentError("arc-family stationarity identity failed"))

    haskey(family, "negative_direction") ||
        throw(ArgumentError("arc family requires a negative-direction witness"))
    direction = _parse_arc_centers(family["negative_direction"], graph.n)
    all(iszero, matrix * vec(direction)) ||
        throw(ArgumentError("arc-family negative direction is not tangent to every contact"))
    second_variation = _arc_second_variation(direction, graph, multipliers)
    iszero(second_variation.q) && length(second_variation.p) == 1 &&
        second_variation.p[1] < 0 ||
        throw(ArgumentError("arc-family second variation is not a negative constant"))
    negative_second_variation = second_variation.p[1]

    objective_expression = _arc_objective(centers)
    iszero(objective_expression.q) && length(objective_expression.p) == 1 ||
        throw(ArgumentError("arc-family objective is not constant"))
    objective = objective_expression.p[1]
    target = second_moment_target(graph.n)
    objective > target || throw(ArgumentError("arc family is not above the target"))
    domain_text = "$(lower) < a < $(upper), s^2 = 4-a^2"
    return FamilyResult(string(entry["critical_id"]), string(entry["graph6"]), graph.n,
        :quadratic_arc, domain_text, num_contacts(graph), length(_all_nonedges(graph)),
        multipliers, objective, target, centroid_verified, contacts_verified,
        noncontacts_verified, :exact_polynomial_positivity, stationarity_verified,
        :full_row_rank_proved, rank_witness_columns,
        :exact_negative_direction, negative_second_variation,
        :excluded_by_exact_negative_intrinsic_direction)
end

"""
    verify_family(entry)

Verify every member of a critical family described by the dataset.  No
floating-point values, sampling, identifier allow-list, or declared objective
value participates in the conclusion.
"""
function verify_family(entry::Dict{String,Any})
    candidate_id = string(entry["critical_id"])
    string(entry["representation"]) == "exact_family" ||
        throw(ArgumentError("$candidate_id must use the `exact_family` representation"))
    family = entry["family"]::Dict{String,Any}
    graph = ContactGraph(string(entry["graph6"]))
    listed_graph = ContactGraph(graph.n, _entry_edges(entry))
    listed_graph.edges == graph.edges || throw(ArgumentError("family contact list mismatch"))

    model = string(family["model"])
    model == "balanced_star" && return _verify_balanced_star_family(entry, family, graph)
    model == "quadratic_arc" && return _verify_arc_family(entry, family, graph)
    model == "four_disk_rhombus" && return _verify_four_disk_rhombus_family(entry, family, graph)
    throw(ArgumentError("unsupported family model: $model"))
end

# Public API alias.
const verify_exact_family = verify_family

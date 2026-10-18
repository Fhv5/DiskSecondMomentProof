_parse_rational(text::AbstractString) = begin
    pieces = split(text, '/')
    length(pieces) == 1 && return parse(BigInt, pieces[1]) // BigInt(1)
    length(pieces) == 2 || throw(ArgumentError("invalid rational: $text"))
    parse(BigInt, pieces[1]) // parse(BigInt, pieces[2])
end

_parse_q3(raw::Dict{String,Any}) = Qsqrt3(_parse_rational(raw["a"]), _parse_rational(raw["b"]))

function _qsqrt3_rank(A::AbstractMatrix{Qsqrt3})
    M = copy(A)
    rows, cols = size(M)
    pivot_row = 1
    for col in 1:cols
        pivot = findfirst(row -> !iszero(M[row, col]), pivot_row:rows)
        pivot === nothing && continue
        selected = pivot_row + pivot - 1
        if selected != pivot_row
            M[pivot_row, :], M[selected, :] = copy(M[selected, :]), copy(M[pivot_row, :])
        end
        scale = M[pivot_row, col]
        M[pivot_row, :] ./= Ref(scale)
        for row in 1:rows
            row == pivot_row && continue
            factor = M[row, col]
            iszero(factor) || (M[row, :] .-= factor .* M[pivot_row, :])
        end
        pivot_row += 1
        pivot_row > rows && break
    end
    return pivot_row - 1
end

function _qsqrt3_linear_solution(A::AbstractMatrix{Qsqrt3}, b::Vector{Qsqrt3})
    rows, cols = size(A)
    length(b) == rows || throw(DimensionMismatch())
    augmented = hcat(A, b)
    pivot_row = 1
    pivot_columns = Int[]
    for col in 1:cols
        pivot = findfirst(row -> !iszero(augmented[row, col]), pivot_row:rows)
        pivot === nothing && continue
        selected = pivot_row + pivot - 1
        if selected != pivot_row
            augmented[pivot_row, :], augmented[selected, :] =
                copy(augmented[selected, :]), copy(augmented[pivot_row, :])
        end
        scale = augmented[pivot_row, col]
        augmented[pivot_row, :] ./= Ref(scale)
        for row in 1:rows
            row == pivot_row && continue
            factor = augmented[row, col]
            iszero(factor) || (augmented[row, :] .-= factor .* augmented[pivot_row, :])
        end
        push!(pivot_columns, col)
        pivot_row += 1
        pivot_row > rows && break
    end
    for row in 1:rows
        all(iszero, augmented[row, 1:cols]) && !iszero(augmented[row, end]) && return nothing
    end
    x = zeros(Qsqrt3, cols)
    for (row, col) in enumerate(pivot_columns)
        x[col] = augmented[row, end]
    end
    return x
end

"""Verification data for an announced rigid minimizer in Q(sqrt(3))."""
struct ExactMinimizerResult
    candidate_id::String
    graph6::String
    n::Int
    contact_rank::Int
    rolling_dimension_modulo_rigid_motions::Int
    # A complete contact graph has no noncontact pairs.
    minimum_noncontact_gap::Union{Nothing,Qsqrt3}
    multipliers::Vector{Qsqrt3}
    objective::Qsqrt3
    target::Qsqrt3
    conclusion::Symbol
end

# Evidence files use this concrete type name.
const TargetResult = ExactMinimizerResult

function _qsqrt3_contact_matrix(centers::Matrix{Qsqrt3}, g::ContactGraph)
    A = zeros(Qsqrt3, num_contacts(g), 2 * g.n)
    for (row, edge) in enumerate(g.edges)
        for coordinate in 1:2
            value = centers[coordinate, edge.src] - centers[coordinate, edge.dst]
            A[row, 2 * edge.src - 2 + coordinate] = value
            A[row, 2 * edge.dst - 2 + coordinate] = -value
        end
    end
    return A
end

function _qsqrt3_second_moment(centers::Matrix{Qsqrt3})
    n = size(centers, 2)
    means = [sum(centers[row, :]) / n for row in 1:2]
    total = zero(Qsqrt3)
    for i in 1:n, row in 1:2
        total += (centers[row, i] - means[row])^2
    end
    return total / 2
end

function _parse_qsqrt3_centers(entry::Dict{String,Any}, n::Int)
    centers = zeros(Qsqrt3, 2, n)
    raw_centers = entry["centros"]::Dict{String,Any}
    for i in 1:n, coordinate in 1:2
        centers[coordinate, i] = _parse_q3(raw_centers[string(i)][coordinate])
    end
    return centers
end

function _parse_edge_multiplier_map(raw::Dict{String,Any}, g::ContactGraph)
    expected_keys = Set("$(edge.src)-$(edge.dst)" for edge in g.edges)
    Set(keys(raw)) == expected_keys ||
        throw(ArgumentError("multiplier keys must match the contact edges"))
    return Qsqrt3[_parse_rational(raw["$(edge.src)-$(edge.dst)"]) for edge in g.edges]
end

"""Verify an announced rigid minimizer in `Q(sqrt(3))`."""
function verify_target(entry::Dict{String,Any})
    candidate_id = string(entry["critical_id"])
    graph6 = string(entry["graph6"])
    string(entry["representation"]) == "exact_qsqrt3" ||
        throw(ArgumentError("$candidate_id is not represented in Q(sqrt(3))"))
    g = ContactGraph(graph6)
    n = g.n
    listed = ContactGraph(n, _entry_edges(entry))
    g.edges == listed.edges || throw(ArgumentError("witness contact list mismatch"))

    centers = _parse_qsqrt3_centers(entry, n)
    all(iszero, vec(sum(centers, dims=2))) || throw(ArgumentError("witness is not centered"))

    edge_set = Set(_edge_tuple.(g.edges))
    minimum_gap = nothing
    for i in 1:(n-1), j in (i+1):n
        distance_squared = sum((centers[row, i] - centers[row, j])^2 for row in 1:2)
        if (i, j) in edge_set
            distance_squared == Qsqrt3(4) || throw(ArgumentError("contact equation failed at ($i,$j)"))
        else
            gap = distance_squared - Qsqrt3(4)
            _q3_sign(gap) > 0 || throw(ArgumentError("noncontact inequality failed at ($i,$j)"))
            minimum_gap = minimum_gap === nothing || _q3_sign(gap - minimum_gap) < 0 ? gap : minimum_gap
        end
    end

    A = _qsqrt3_contact_matrix(centers, g)
    rank_A = _qsqrt3_rank(A)
    gradient_vector = vec(centers)
    lambda = if haskey(entry, "multipliers")
        _parse_edge_multiplier_map(entry["multipliers"], g)
    else
        solution = _qsqrt3_linear_solution(transpose(A), gradient_vector)
        solution === nothing && throw(ArgumentError("stationarity system is inconsistent"))
        solution
    end
    transpose(A) * lambda == gradient_vector || throw(ArgumentError("multiplier balance failed"))
    all(value -> _q3_sign(value) > 0, lambda) || throw(ArgumentError("multipliers are not strictly positive"))

    objective = _qsqrt3_second_moment(centers)
    target = Qsqrt3(second_moment_target(n))
    objective == target || throw(ArgumentError("second moment does not equal target"))
    rolling_dimension = 2n - rank_A - 3
    rolling_dimension == 0 || throw(ArgumentError("announced minimizer is not rolling rigid"))

    return TargetResult(candidate_id, graph6, n, rank_A, rolling_dimension,
                             minimum_gap, lambda, objective, target,
                             :target_minimizer_exactly_verified)
end

# Public API alias.
const verify_exact_minimizer = verify_target

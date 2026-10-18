@inline _edge_tuple(e::ContactEdge) = (e.src, e.dst)

const _SUPPORTED_REPRESENTATIONS = Set(("decimal", "exact_qsqrt3", "exact_family"))

function _candidate_representation(data::Dict{String,Any}, entry::Dict{String,Any})
    representation = string(get(entry, "representation", get(data, "default_representation", "decimal")))
    representation in _SUPPORTED_REPRESENTATIONS ||
        throw(ArgumentError("unsupported candidate representation: $representation"))
    return representation
end

_is_rational_string(value) = value isa AbstractString &&
    occursin(r"^[+-]?[0-9]+(?:/[1-9][0-9]*)?$", value)

function _validate_center_labels(centers::Dict{String,Any}, n::Int, candidate_id::String)
    Set(keys(centers)) == Set(string.(1:n)) ||
        throw(ArgumentError("$candidate_id: center labels must be exactly 1:$n"))
end

function _validate_decimal_centers(centers, n::Int, candidate_id::String; field="centros")
    centers isa Dict{String,Any} ||
        throw(ArgumentError("$candidate_id: `$field` must be an object"))
    _validate_center_labels(centers, n, candidate_id)
    for i in 1:n
        xy = centers[string(i)]
        xy isa Vector && length(xy) == 2 ||
            throw(ArgumentError("$candidate_id: center $i must have two coordinates"))
        all(value -> value isa AbstractString, xy) ||
            throw(ArgumentError("$candidate_id: approximate coordinates must be decimal strings, never JSON numbers"))
    end
end


function _validate_qsqrt3_centers(centers, n::Int, candidate_id::String)
    centers isa Dict{String,Any} ||
        throw(ArgumentError("$candidate_id: `centros` must be an object"))
    _validate_center_labels(centers, n, candidate_id)
    for i in 1:n
        xy = centers[string(i)]
        xy isa Vector && length(xy) == 2 ||
            throw(ArgumentError("$candidate_id: center $i must have two coordinates"))
        for coordinate in xy
            coordinate isa Dict{String,Any} && Set(keys(coordinate)) == Set(("a", "b")) ||
                throw(ArgumentError("$candidate_id: every Q(sqrt(3)) coordinate must contain rational a and b"))
            _is_rational_string(coordinate["a"]) && _is_rational_string(coordinate["b"]) ||
                throw(ArgumentError("$candidate_id: invalid rational coefficient in Q(sqrt(3)) coordinate"))
        end
    end
end

function _entry_edges(entry::Dict{String,Any})
    haskey(entry, "contactos") || throw(ArgumentError("candidate is missing `contactos`"))
    result = ContactEdge[]
    for raw_edge in entry["contactos"]
        length(raw_edge) == 2 || throw(ArgumentError("every contact must contain two vertices"))
        raw_edge[1] isa Integer && raw_edge[2] isa Integer ||
            throw(ArgumentError("contact vertices must be integers"))
        push!(result, ContactEdge(Int(raw_edge[1]), Int(raw_edge[2])))
    end
    return result
end

"""
    validate_candidate_dataset(data)

Checks the structural and combinatorial integrity of a candidate dataset. It
does not infer contacts from approximate distances and makes no geometric
conclusion.
"""
function validate_candidate_dataset(data::Dict{String,Any})
    for key in ("discos", "graphs", "total_entries")
        haskey(data, key) || throw(ArgumentError("dataset is missing `$key`"))
    end

    data["discos"] isa Integer || throw(ArgumentError("`discos` must be an integer"))
    n = Int(data["discos"])
    entries = data["graphs"]
    entries isa Vector || throw(ArgumentError("`graphs` must be an array"))
    Int(data["total_entries"]) == length(entries) ||
        throw(ArgumentError("`total_entries` does not match the number of entries"))

    seen_ids = Set{String}()
    for (entry_index, entry_any) in enumerate(entries)
        entry_any isa Dict{String,Any} ||
            throw(ArgumentError("entry $entry_index is not a JSON object"))
        entry = entry_any::Dict{String,Any}

        for key in ("critical_id", "graph6", "contactos")
            haskey(entry, key) || throw(ArgumentError("entry $entry_index is missing `$key`"))
        end

        candidate_id = string(entry["critical_id"])
        isempty(candidate_id) && throw(ArgumentError("entry $entry_index has an empty critical_id"))
        candidate_id in seen_ids && throw(ArgumentError("duplicate critical_id: $candidate_id"))
        push!(seen_ids, candidate_id)

        decoded = ContactGraph(string(entry["graph6"]))
        decoded.n == n ||
            throw(ArgumentError("$candidate_id: graph6 has $(decoded.n) vertices, expected $n"))

        listed_edges = _entry_edges(entry)
        listed_graph = ContactGraph(n, listed_edges)
        length(listed_edges) == length(listed_graph.edges) ||
            throw(ArgumentError("$candidate_id: repeated contact edge"))
        Set(_edge_tuple.(decoded.edges)) == Set(_edge_tuple.(listed_graph.edges)) ||
            throw(ArgumentError("$candidate_id: `contactos` does not match graph6"))

        representation = _candidate_representation(data, entry)
        if representation == "decimal"
            haskey(entry, "centros") || throw(ArgumentError("$candidate_id: decimal candidate is missing `centros`"))
            _validate_decimal_centers(entry["centros"], n, candidate_id)
        elseif representation == "exact_qsqrt3"
            haskey(entry, "centros") || throw(ArgumentError("$candidate_id: Q(sqrt(3)) candidate is missing `centros`"))
            _validate_qsqrt3_centers(entry["centros"], n, candidate_id)
        elseif representation == "exact_family"
            haskey(entry, "family") && entry["family"] isa Dict{String,Any} ||
                throw(ArgumentError("$candidate_id: family is missing its `family` object"))
            haskey(entry, "centros") &&
                throw(ArgumentError("$candidate_id: families must not expose a sample as proof coordinates"))
            if haskey(entry, "sample_centers")
                _validate_decimal_centers(entry["sample_centers"], n, candidate_id; field="sample_centers")
            end
        end
    end
    return true
end

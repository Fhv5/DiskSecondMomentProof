module Census

using ..DiskSecondMomentProof: ContactGraph, ContactEdge, connected_graph_path
import ..DiskSecondMomentProof

export is_outerplanar, filter_planar_graphs, compare_with_stored_census, verify_four_disk_planar_census

function _validate_graph6(code::AbstractString, n::Int)
    bytes = codeunits(code)
    isempty(bytes) && throw(ArgumentError("graph6 code is empty"))
    bit_count = n * (n - 1) ÷ 2
    length(bytes) == 1 + cld(bit_count, 6) ||
        throw(ArgumentError("invalid graph6 length for n=$n: $code"))
    all(byte -> 63 <= byte <= 126, bytes) ||
        throw(ArgumentError("non-graph6 character: $code"))
    Int(bytes[1]) - 63 == n ||
        throw(ArgumentError("graph6 order does not equal n=$n: $code"))
    unused = 6 * cld(bit_count, 6) - bit_count
    unused == 0 || (Int(bytes[end]) - 63) & ((1 << unused) - 1) == 0 ||
        throw(ArgumentError("nonzero graph6 padding: $code"))
    return ContactGraph(code)
end

function _connected(graph::ContactGraph)
    reached = falses(graph.n)
    reached[1] = true
    changed = true
    while changed
        changed = false
        for edge in graph.edges
            if reached[edge.src] != reached[edge.dst]
                reached[edge.src] = reached[edge.dst] = true
                changed = true
            end
        end
    end
    return all(reached)
end

function _orders(n::Int, visit; fixed_first::Bool=true)
    order = collect(1:n)
    function permute(position::Int)
        if position > n
            return visit(order)
        end
        for index in position:n
            order[position], order[index] = order[index], order[position]
            permute(position + 1) && return true
            order[position], order[index] = order[index], order[position]
        end
        return false
    end
    return permute(fixed_first ? 2 : 1)
end

function _chords_do_not_cross(graph::ContactGraph, order::Vector{Int})
    position = zeros(Int, graph.n)
    for (index, vertex) in enumerate(order)
        position[vertex] = index
    end
    for i in 1:(length(graph.edges) - 1), j in (i + 1):length(graph.edges)
        first = graph.edges[i]
        second = graph.edges[j]
        length(Set((first.src, first.dst, second.src, second.dst))) == 4 || continue
        a, b = sort((position[first.src], position[first.dst]))
        c, d = sort((position[second.src], position[second.dst]))
        (a < c < b < d || c < a < d < b) && return false
    end
    return true
end

"""Decide outerplanarity by testing every circular vertex order, for n=2, n=3, n=4, n=5 or n=6."""
function is_outerplanar(graph::ContactGraph)
    graph.n in (2, 3, 4, 5, 6) || throw(ArgumentError("outerplanarity census supports n=2,3,4,5,6"))
    return _orders(graph.n, order -> _chords_do_not_cross(graph, order))
end

function _canonical_signature(graph::ContactGraph)
    edge_set = Set((edge.src, edge.dst) for edge in graph.edges)
    best = Ref(typemax(UInt64))
    _orders(graph.n, order -> begin
        value = UInt64(0)
        for j in 2:graph.n, i in 1:(j - 1)
            left, right = minmax(order[i], order[j])
            value = (value << 1) | UInt64((left, right) in edge_set)
        end
        best[] = min(best[], value)
        return false
    end; fixed_first=false)
    return best[]
end

"""Filter the graph6 output of plantri/shortg using the Julia graph tests."""
function filter_planar_graphs(raw_codes::AbstractVector{<:AbstractString}, n::Int)
    n in (2, 3, 4, 5, 6) || throw(ArgumentError("the census supports n=2,3,4,5,6"))
    codes = String.(raw_codes)
    graphs = [_validate_graph6(code, n) for code in codes]
    all(_connected, graphs) || throw(ArgumentError("input contains a disconnected graph"))
    signatures = _canonical_signature.(graphs)
    length(unique(signatures)) == length(signatures) ||
        throw(ArgumentError("input contains isomorphic duplicate graphs"))
    return [code for (code, graph) in zip(codes, graphs) if is_outerplanar(graph)]
end

"""Compare filtered graph classes with the stored proof input, ignoring labels."""
function compare_with_stored_census(codes::AbstractVector{<:AbstractString}, n::Int)
    stored = filter(!isempty, strip.(readlines(connected_graph_path(n))))
    length(codes) == length(stored) || return false
    generated_keys = Set(_canonical_signature(_validate_graph6(code, n)) for code in codes)
    stored_keys = Set(_canonical_signature(_validate_graph6(code, n)) for code in stored)
    return length(generated_keys) == length(codes) && generated_keys == stored_keys
end

"""Enumerate every connected four-vertex class and verify the planar K4 obstruction."""
function verify_four_disk_planar_census(raw_codes)
    codes = String.(raw_codes)
    graphs = [_validate_graph6(code, 4) for code in codes]
    all(_connected, graphs) || throw(ArgumentError("four-disk raw census contains a disconnected graph"))
    pairs = [ContactEdge(i, j) for i in 1:3 for j in (i + 1):4]
    complete = Set{UInt64}()
    for mask in 0:63
        edges = [edge for (index, edge) in enumerate(pairs) if mask & (1 << (index - 1)) != 0]
        graph = ContactGraph(4, edges)
        _connected(graph) && push!(complete, _canonical_signature(graph))
    end
    signatures = _canonical_signature.(graphs)
    length(signatures) == length(unique(signatures)) == length(complete) == 6 &&
        Set(signatures) == complete ||
        throw(ArgumentError("four-disk raw census does not contain exactly all six connected classes"))
    k4 = only([code for (code, graph) in zip(codes, graphs) if length(graph.edges) == 6])
    squared_contact = BigInt(4)//BigInt(1)
    off_diagonal = (squared_contact + squared_contact - squared_contact)/2
    gram = [DiskSecondMomentProof.Qsqrt3(i == j ? squared_contact : off_diagonal) for i in 1:3, j in 1:3]
    rank = DiskSecondMomentProof._qsqrt3_rank(gram)
    rank == 3 || error("K4 Gram obstruction failed")
    return (status=:complete_connected_census, connected_classes=6,
            k4_graph6=k4, gram_rank=rank, ambient_dimension=2)
end

end

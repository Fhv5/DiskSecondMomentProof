"""
    ContactEdge

Represents an undirected active contact between two disks (src, dst) with src < dst.
"""
struct ContactEdge
    src::Int
    dst::Int

    function ContactEdge(u::Int, v::Int)
        @assert u != v "Self-loops are not valid contact edges."
        return u < v ? new(u, v) : new(v, u)
    end
end

"""
    ContactGraph

Represents the contact graph of an n-disk cluster with a list of unique ContactEdges.
"""
struct ContactGraph
    n::Int
    edges::Vector{ContactEdge}

    function ContactGraph(n::Int, edges::Vector{ContactEdge})
        unique_edges = sort(unique(edges), by = e -> (e.src, e.dst))
        for e in unique_edges
            @assert 1 <= e.src <= n && 1 <= e.dst <= n "Edge vertices must be within 1..n."
        end
        return new(n, unique_edges)
    end
end

"""
    decode_graph6(graph6_str::AbstractString)

Decodes a graph6 format string into the number of vertices `n` and a list of `ContactEdge`s.
Parse the combinatorial contact list without distance thresholds.
"""
function decode_graph6(graph6_str::AbstractString)
    clean_s = strip(graph6_str)
    @assert !isempty(clean_s) "graph6 string cannot be empty."
    bytes = [Int(c) - 63 for c in clean_s]
    n = bytes[1]
    @assert 1 <= n <= 62 "Unsupported vertex count in simple graph6 header: n = $n."

    bits = Int[]
    for b in bytes[2:end]
        for shift in 5:-1:0
            push!(bits, (b >> shift) & 1)
        end
    end

    edges = ContactEdge[]
    idx = 1
    for j in 1:n
        for i in 1:(j - 1)
            if idx <= length(bits) && bits[idx] == 1
                push!(edges, ContactEdge(i, j))
            end
            idx += 1
        end
    end
    return n, edges
end

"""
    ContactGraph(graph6_str::AbstractString)

Constructs a ContactGraph directly from a graph6 format string.
"""
function ContactGraph(graph6_str::AbstractString)
    n, edges = decode_graph6(graph6_str)
    return ContactGraph(n, edges)
end

"""
    num_contacts(g::ContactGraph)

Returns the total number of active contact edges m in the graph.
"""
function num_contacts(g::ContactGraph)
    return length(g.edges)
end

"""
    degree(g::ContactGraph, v::Int)

Returns the degree (number of incident contacts) of disk v in graph g.
"""
function degree(g::ContactGraph, v::Int)
    return count(e -> e.src == v || e.dst == v, g.edges)
end

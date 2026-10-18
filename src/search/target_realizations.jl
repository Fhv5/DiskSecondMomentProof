module TargetRealizations

using ...DiskSecondMomentProof: Qsqrt3,
                                    ContactGraph,
                                    ContactEdge,
                                    second_moment_target,
                                    read_json_file,
                                    _q3_sign

export Point,
       triangle_construction,
       enumerate_embeddings,
       classify_embedding,
       geometric_classes,
       analyze_target,
       second_moment,
       ExactPoint,
       exact_objective

const Q3 = Qsqrt3
const SQRT3 = Q3(0//1, 1//1)

"""Planar point with coordinates in `Q(sqrt(3))`."""
struct ExactPoint
    x::Q3
    y::Q3
end

# Evidence files use this concrete type name.
const Point = ExactPoint

Base.:+(p::Point, q::Point) = Point(p.x + q.x, p.y + q.y)
Base.:-(p::Point, q::Point) = Point(p.x - q.x, p.y - q.y)
Base.:/(p::Point, k::Integer) = Point(p.x / k, p.y / k)

square(x) = x * x
distance_squared(p::Point, q::Point) = square(p.x - q.x) + square(p.y - q.y)

function qsqrt3_sign(x::Q3)
    return _q3_sign(x)
end

edge_tuple(edge::ContactEdge) = (edge.src, edge.dst)
edge_set(graph::ContactGraph) = Set(edge_tuple(edge) for edge in graph.edges)
canonical_edge(i::Int, j::Int) = i < j ? (i, j) : (j, i)

function has_edge(edges::Set{Tuple{Int,Int}}, i::Int, j::Int)
    return canonical_edge(i, j) in edges
end

function graph_triangles(graph::ContactGraph)
    edges = edge_set(graph)
    triangles = NTuple{3,Int}[]
    for i in 1:(graph.n - 2), j in (i + 1):(graph.n - 1), k in (j + 1):graph.n
        has_edge(edges, i, j) && has_edge(edges, i, k) && has_edge(edges, j, k) &&
            push!(triangles, (i, j, k))
    end
    return triangles
end

"""
Find a triangle seed and an order in which every remaining vertex is attached
to both endpoints of an already constructed contact edge.  The chosen edges
form a spanning 2-tree subgraph; other graph edges are checked afterwards.
"""
function triangle_construction(graph::ContactGraph)
    edges = edge_set(graph)

    function extend(placed::Vector{Int}, remaining::Vector{Int}, steps)
        isempty(remaining) && return copy(steps)
        for vertex in remaining
            for a_index in 1:(length(placed) - 1), b_index in (a_index + 1):length(placed)
                u, v = placed[a_index], placed[b_index]
                if has_edge(edges, u, v) && has_edge(edges, vertex, u) && has_edge(edges, vertex, v)
                    next_placed = [placed; vertex]
                    next_remaining = filter(!=(vertex), remaining)
                    push!(steps, (vertex, u, v))
                    result = extend(next_placed, next_remaining, steps)
                    result === nothing || return result
                    pop!(steps)
                end
            end
        end
        return nothing
    end

    for seed in graph_triangles(graph)
        placed = collect(seed)
        remaining = [vertex for vertex in 1:graph.n if !(vertex in placed)]
        steps = extend(placed, remaining, Tuple{Int,Int,Int}[])
        steps === nothing || return seed, steps
    end
    error("graph $(graph.n)-vertex graph has no spanning triangle construction")
end

function equilateral_third_points(p::Point, q::Point)
    midpoint = (p + q) / 2
    dx = q.x - p.x
    dy = q.y - p.y
    offset = Point(-(SQRT3 * dy) / 2, (SQRT3 * dx) / 2)
    return midpoint + offset, midpoint - offset
end

"""Enumerate the branches of a contact edge or a fixed triangle construction."""
function enumerate_embeddings(graph::ContactGraph, seed, steps)
    if graph.n == 2
        seed == (1, 2) && isempty(steps) && graph.edges == [ContactEdge(1, 2)] ||
            throw(ArgumentError("the two-disk construction must be the single contact edge"))
        return [[Point(Q3(0), Q3(0)), Point(Q3(2), Q3(0))]]
    end
    coordinates = Dict{Int,Point}(
        seed[1] => Point(Q3(0), Q3(0)),
        seed[2] => Point(Q3(2), Q3(0)),
        seed[3] => Point(Q3(1), SQRT3),
    )
    embeddings = Vector{Vector{Point}}()

    function visit(step_index::Int)
        if step_index > length(steps)
            push!(embeddings, [coordinates[i] for i in 1:graph.n])
            return
        end
        vertex, u, v = steps[step_index]
        for point in equilateral_third_points(coordinates[u], coordinates[v])
            coordinates[vertex] = point
            visit(step_index + 1)
        end
        delete!(coordinates, vertex)
    end

    visit(1)
    return embeddings
end

"""Classify a bar-framework embedding against the declared penny graph."""
function classify_embedding(graph::ContactGraph, points::Vector{Point})
    edges = edge_set(graph)
    for (i, j) in edges
        distance_squared(points[i], points[j]) == Q3(4) || return :failed_required_contact
    end

    has_overlap = false
    has_extra_contact = false
    for i in 1:(graph.n - 1), j in (i + 1):graph.n
        has_edge(edges, i, j) && continue
        gap = distance_squared(points[i], points[j]) - Q3(4)
        sign = qsqrt3_sign(gap)
        has_overlap |= sign < 0
        has_extra_contact |= sign == 0
    end
    has_overlap && return :overlap
    has_extra_contact && return :extra_contact
    return :exact_penny_realization
end

function all_permutations(n::Int)
    result = Vector{Vector{Int}}()
    current = Int[]
    used = falses(n)
    function visit()
        if length(current) == n
            push!(result, copy(current))
            return
        end
        for value in 1:n
            used[value] && continue
            used[value] = true
            push!(current, value)
            visit()
            pop!(current)
            used[value] = false
        end
    end
    visit()
    return result
end

function graph_automorphisms(graph::ContactGraph)
    edges = edge_set(graph)
    automorphisms = Vector{Vector{Int}}()
    for permutation in all_permutations(graph.n)
        preserves_edges = true
        for i in 1:(graph.n - 1), j in (i + 1):graph.n
            if has_edge(edges, i, j) != has_edge(edges, permutation[i], permutation[j])
                preserves_edges = false
                break
            end
        end
        preserves_edges && push!(automorphisms, permutation)
    end
    return automorphisms
end

function geometrically_equivalent(
    first::Vector{Point},
    second::Vector{Point},
    automorphisms::Vector{Vector{Int}},
)
    for permutation in automorphisms
        same_distances = true
        for i in 1:(length(first) - 1), j in (i + 1):length(first)
            if distance_squared(first[i], first[j]) !=
                    distance_squared(second[permutation[i]], second[permutation[j]])
                same_distances = false
                break
            end
        end
        if same_distances
            return true
        end
    end
    return false
end

"""Group embeddings modulo congruence and graph automorphisms."""
function geometric_classes(embeddings, automorphisms)
    representatives = Vector{Vector{Point}}()
    memberships = Vector{Vector{Int}}()
    for (index, embedding) in enumerate(embeddings)
        class_index = findfirst(rep -> geometrically_equivalent(embedding, rep, automorphisms), representatives)
        if class_index === nothing
            push!(representatives, embedding)
            push!(memberships, [index])
        else
            push!(memberships[class_index], index)
        end
    end
    return representatives, memberships
end

function parse_rational(raw::AbstractString)
    pieces = split(raw, '/')
    length(pieces) == 1 && return parse(BigInt, pieces[1]) // BigInt(1)
    length(pieces) == 2 || error("invalid rational: $raw")
    return parse(BigInt, pieces[1]) // parse(BigInt, pieces[2])
end

function dataset_target(path::String)
    data = read_json_file(path)
    target_entries = filter(
        entry -> get(entry, "representation", "decimal") == "exact_qsqrt3",
        data["graphs"],
    )
    length(target_entries) == 1 || error("expected exactly one target in $path")
    entry = only(target_entries)
    n = Int(data["discos"])
    centers = Vector{Point}(undef, n)
    for vertex in 1:n
        raw = entry["centros"][string(vertex)]
        centers[vertex] = Point(
            Q3(parse_rational(raw[1]["a"]), parse_rational(raw[1]["b"])),
            Q3(parse_rational(raw[2]["a"]), parse_rational(raw[2]["b"])),
        )
    end
    return (
        n=n,
        graph6=String(entry["graph6"]),
        candidate_id=String(entry["critical_id"]),
        centers=centers,
    )
end

"""Evaluate the centroidal second moment for points in Q(sqrt(3))."""
function second_moment(points::Vector{Point})
    total = Q3(0)
    for i in 1:(length(points) - 1), j in (i + 1):length(points)
        total += distance_squared(points[i], points[j])
    end
    return total / (2 * length(points))
end

"""Enumerate and classify every Q(sqrt(3)) realization constructed for the dataset target graph."""
function analyze_target(path::String)
    target = dataset_target(path)
    graph = ContactGraph(target.graph6)
    seed, steps = if graph.n == 2
        graph.edges == [ContactEdge(1, 2)] ||
            throw(ArgumentError("the two-disk target must have exactly one contact"))
        ((1, 2), NTuple{3,Int}[])
    else
        triangle_construction(graph)
    end
    embeddings = enumerate_embeddings(graph, seed, steps)
    classifications = classify_embedding.(Ref(graph), embeddings)
    valid = embeddings[classifications .== :exact_penny_realization]
    automorphisms = graph_automorphisms(graph)
    representatives, memberships = geometric_classes(valid, automorphisms)
    target_found = any(
        geometrically_equivalent(target.centers, representative, automorphisms)
        for representative in representatives
    )

    return (
        target=target,
        graph=graph,
        seed=seed,
        steps=steps,
        total_embeddings=length(embeddings),
        classifications=classifications,
        valid=valid,
        automorphisms=automorphisms,
        representatives=representatives,
        memberships=memberships,
        target_found=target_found,
    )
end

# Public API alias.
const exact_objective = second_moment

end

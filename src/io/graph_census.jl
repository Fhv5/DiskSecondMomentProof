"""Load the supported connected graph6 census and check combinatorial integrity."""
function load_connected_graphs(n::Int)
    codes = filter(!isempty, strip.(readlines(connected_graph_path(n))))
    expected = Dict(2 => 1, 3 => 2, 4 => 5, 5 => 13, 6 => 46)[n]
    length(codes) == expected ||
        throw(ArgumentError("the n=$n census has $(length(codes)) entries, expected $expected"))
    length(unique(codes)) == length(codes) ||
        throw(ArgumentError("the n=$n census contains duplicate graph6 strings"))
    graphs = ContactGraph.(codes)
    all(graph -> graph.n == n, graphs) ||
        throw(ArgumentError("the n=$n census contains a graph of another order"))
    length(Census.filter_planar_graphs(codes, n)) == length(codes) ||
        throw(ArgumentError("the n=$n census contains a non-outerplanar graph"))
    return codes
end

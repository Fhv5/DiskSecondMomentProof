module SearchConventions

using JSON3
using ..DiskSecondMomentProof: BoxRecords
using ..CandidateSearch: LeafExclusion, LeaflessSearch, BoundaryRoots

const LE = LeafExclusion
const LS = LeaflessSearch
const BR = BoundaryRoots

edge_arrays(edges) = [collect(edge) for edge in edges]
canonical_edge(edge) = minmax(edge...)

"""Describe the deterministic chart used by one angular search."""
function search_chart(graph; leaf::Union{Nothing,Int}=nothing)
    tree = leaf === nothing ? LS.rooted_tree(graph) : LE.anchored_tree(graph, leaf)
    k = graph.n - 2
    reflection, constraints = leaf === nothing ? (1, Tuple{Int,Int}[]) :
        LE.angular_symmetry_data(graph, tree)
    box = leaf === nothing ? LS.initial_angle_box(graph.n) :
        LE.initial_box(graph.n, reflection)
    tree_edges = Set(canonical_edge.(tree))
    extra_edges = [edge for edge in graph.edges if !(edge in tree_edges)]
    permutations = Dict{String,Any}[]
    if !isempty(constraints)
        # Include swaps involving the reflection coordinate too (pattern P3).
        for i in 1:k, j in (i + 1):k
            tree[i + 1][1] == tree[j + 1][1] || continue
            u, v = tree[i + 1][2], tree[j + 1][2]
            LE.vertex_swap_is_automorphism(graph, u, v) || continue
            push!(permutations, Dict(
                "angle_indices" => [i, j],
                "vertex_permutation" => [LE.swapped_vertex(a, u, v) for a in 1:graph.n],
            ))
        end
    end
    return Dict{String,Any}(
        "selected_leaf" => leaf,
        "anchor_edge" => collect(first(tree)),
        "anchor_centers" => [[0, 0], [2, 0]],
        "oriented_tree" => edge_arrays(tree),
        "angle_edges" => edge_arrays(tree[2:end]),
        "contact_edges" => edge_arrays(graph.edges),
        "non_tree_contacts" => edge_arrays(extra_edges),
        "reduced_multiplier_edges" => leaf === nothing && length(extra_edges) < k ?
            edge_arrays(extra_edges) : Vector{Int}[],
        "root_system" => leaf !== nothing ? nothing :
            length(extra_edges) == k ? "contact_equations" : "reduced_kkt",
        "reflection_angle_index" => k == 0 ? nothing : reflection,
        "angle_order_constraints" => edge_arrays(constraints),
        "symmetry_pattern" => isempty(constraints) ? "reflection_only" :
            string(LE.VERIFIED_SYMMETRY_PATTERNS[graph.code]),
        "sibling_transpositions" => permutations,
        "initial_angle_box" => Dict(
            "lower" => BR.rational_text.(box.lower),
            "upper" => BR.rational_text.(box.upper),
        ),
    )
end

function attach_trace!(chart, output_directory, stage, code; leaf=nothing)
    path = BoxRecords.box_record_path(joinpath(output_directory, "evidence"), stage, code; leaf=leaf)
    isfile(path) || error("missing box record for $code")
    chart["box_records_file"] = replace(relpath(path, output_directory), '\\' => '/')
    return chart
end

"""Write conventions for the searches actually performed, including partial runs."""
function write_search_conventions(output_directory, n, all_codes, leaf_result,
                                  search_results, target_graph, family_entries)
    leaf_graphs = Dict(result.graph6 => result for result in leaf_result.graph_results)
    angular_graphs = Dict(result.graph6 => result for result in search_results)
    family_graphs = Set(string(entry["graph6"]) for entry in family_entries)
    graphs = Dict{String,Any}[]
    for code in all_codes
        entry = Dict{String,Any}("graph6" => code, "searches" => Dict{String,Any}[])
        if haskey(leaf_graphs, code)
            result = leaf_graphs[code]
            entry["route"] = "leaf_exclusion"
            for run in result.runs
                chart = search_chart(LE.graph_data(code); leaf=run.leaf)
                chart["status"] = string(run.status)
                attach_trace!(chart, output_directory, "leaf", code; leaf=run.leaf)
                push!(entry["searches"], chart)
            end
        elseif haskey(angular_graphs, code)
            entry["route"] = "angular_search"
            chart = search_chart(LE.graph_data(code))
            chart["status"] = string(angular_graphs[code].status)
            attach_trace!(chart, output_directory, "angular", code)
            if !isempty(BR.registered_witnesses(code))
                chart["boundary_witnesses_file"] = "search/angle_boxes/n$n-boundary-roots.tsv"
            end
            push!(entry["searches"], chart)
        elseif code == target_graph || code in family_graphs
            entry["route"] = code == target_graph ? "algebraic_target" : "algebraic_family"
            entry["input_file"] = "data/exact_inputs/$(n)disks.json"
        else
            entry["route"] = "not_searched"
        end
        push!(graphs, entry)
    end
    record = Dict(
        "schema_version" => 1,
        "n" => n,
        "index_base" => 1,
        "angle_unit" => "radians",
        "angle_definition" => "c_child - c_parent = 2*(cos(theta), sin(theta)); principal angles modulo 2*pi",
        "reflection_action" => "theta_i -> -theta_i modulo 2*pi for all i",
        "order_constraint_definition" => "[i,j] means theta_i <= theta_j",
        "permutation_action" => "p[old] = new; c'_p[i] = c_i; recover each angle from its oriented tree edge",
        "full_multiplier_order" => "contact_edges",
        "root_state_order" => "angle_edges, then reduced_multiplier_edges when root_system is reduced_kkt",
        "subdivision" => "bisect a widest interval at its rational midpoint; break ties by lowest angle index; depth-first, lower child first",
        "graphs" => graphs,
    )
    path = joinpath(output_directory, "search", "conventions.json")
    mkpath(dirname(path))
    open(path, "w") do io
        JSON3.write(io, record)
        println(io)
    end
    return path
end

end

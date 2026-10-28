using Test
using DiskSecondMomentProof

@testset "Search conventions across the two-to-six-disk census" begin
    SC = DiskSecondMomentProof.SearchConventions
    LE = CandidateSearch.LeafFilters
    patterns = Dict{String,Int}()
    for n in 2:6, code in load_connected_graphs(n)
        graph = LE.graph_data(code)
        n == 2 && continue # The two-disk graph follows the algebraic route.
        leaf = isempty(graph.leaves) ? nothing : first(graph.leaves)
        chart = SC.search_chart(graph; leaf=leaf)
        tree = Tuple.(chart["oriented_tree"])
        @test length(tree) == n - 1
        @test chart["anchor_edge"] == collect(first(tree))
        @test chart["angle_edges"] == [collect(edge) for edge in tree[2:end]]
        # Independently check that the tree places every vertex once, from its anchor.
        placed = Set(first(tree))
        for (parent, child) in tree[2:end]
            @test parent in placed && !(child in placed)
            push!(placed, child)
        end
        @test placed == Set(1:n)
        contacts = Set(Tuple.(chart["contact_edges"]))
        @test Set(minmax.(first.(tree), last.(tree))) ⊆ contacts
        @test Set(Tuple.(chart["non_tree_contacts"])) ==
            setdiff(contacts, Set(minmax.(first.(tree), last.(tree))))
        reflection = chart["reflection_angle_index"]
        @test chart["initial_angle_box"]["lower"][reflection] == "0"
        @test all(==("22/7"), chart["initial_angle_box"]["upper"])
        for swap in chart["sibling_transpositions"]
            p = swap["vertex_permutation"]
            @test sort(p) == collect(1:n)
            @test p[chart["anchor_edge"]] == chart["anchor_edge"]
            @test Set(minmax(p[u], p[v]) for (u, v) in contacts) == contacts
            i, j = swap["angle_indices"]
            @test tree[i + 1][1] == tree[j + 1][1]
            @test p[tree[i + 1][2]] == tree[j + 1][2]
        end
        for (i, j) in chart["angle_order_constraints"]
            @test any(s -> s["angle_indices"] == [i, j], chart["sibling_transpositions"])
        end
        if leaf !== nothing
            key = chart["symmetry_pattern"]
            patterns[key] = get(patterns, key, 0) + 1
            if key == "P3"
                @test any(s -> reflection in s["angle_indices"], chart["sibling_transpositions"])
                @test all(pair -> !(reflection in pair), chart["angle_order_constraints"])
            end
        end
    end
    @test patterns == Dict("reflection_only" => 22, "P0" => 12, "P1" => 4, "P2" => 4, "P3" => 4)
end

@testset "EPVW witness coordinate and multiplier order" begin
    SC = DiskSecondMomentProof.SearchConventions
    graph = CandidateSearch.LeafFilters.graph_data("EPVW")
    chart = SC.search_chart(graph)
    @test chart["anchor_edge"] == [1, 3]
    @test chart["angle_edges"] == [[1, 6], [3, 4], [6, 2], [6, 5]]
    @test chart["reduced_multiplier_edges"] == [[2, 5], [4, 5], [4, 6]]
    @test chart["root_system"] == "reduced_kkt"
    @test isempty(chart["angle_order_constraints"])
    witness = only(CandidateSearch.BoundaryRoots.registered_witnesses("EPVW"))
    verified = CandidateSearch.BoundaryRoots.verify_witness(graph, Tuple.(chart["oriented_tree"]), witness)
    @test verified.determinant == Qsqrt3(-4608)
    # A non-involutive relabeling distinguishes old-to-new from its inverse.
    p = witness.permutation
    inverse = invperm(p)
    @test p != inverse
    augmented = vcat(graph.edges, witness.extra_contacts)
    target = Set(CandidateSearch.LeafFilters.graph_data("EQlw").edges)
    @test Set(minmax(p[u], p[v]) for (u, v) in augmented) == target
    @test Set(minmax(inverse[u], inverse[v]) for (u, v) in augmented) != target
    @test minmax(p[1], p[4]) == (4, 6)
end

@testset "Convention files describe actual searches and are inventoried" begin
    for n in (2, 3)
        mktempdir() do output
            result = DiskSecondMomentProof.ProofWorkflow.execute([
                "--n", string(n), "--precision", "128", "--record-evidence", "--output", output,
            ])
            directory = joinpath(output, "n$n")
            run = DiskSecondMomentProof.read_json_file(joinpath(directory, "run.json"))
            path = joinpath(directory, run["search_conventions_file"])
            @test run["files_sha256"]["search/conventions.json"] ==
                DiskSecondMomentProof.EvidenceArchive.sha256_file(path)
            conventions = DiskSecondMomentProof.read_json_file(path)
            @test conventions["schema_version"] == 2
            @test occursin("p[old] = new; c'_p[i] = c_i", conventions["permutation_action"])
            @test length(conventions["graphs"]) == result.connected_census
            @test count(g -> g["route"] == "algebraic_target", conventions["graphs"]) == 1
            searches = [s for g in conventions["graphs"] for s in g["searches"]]
            @test length(searches) == n - 2
            for chart in searches
                @test isfile(joinpath(directory, chart["box_records_file"]))
                @test chart["selected_leaf"] == only(result.search_results).chart.selected_leaf
                @test chart["status"] == "complete"
            end
        end
    end
end

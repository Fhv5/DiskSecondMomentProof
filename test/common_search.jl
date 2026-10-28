@testset "Common search limits and leaf access" begin
    D = DiskSecondMomentProof
    LS = D.CandidateSearch.IntervalSearch
    options = LS.SearchOptions()
    @test options.max_depth == 140
    @test options.max_boxes == 3_000_000
    @test options.minimum_width == 1 // big(2)^60
    @test options.krawczyk_depth == 12
    @test_throws ArgumentError LS.SearchOptions(max_boxes=0)
    @test_throws ArgumentError LS.SearchOptions(max_depth=-1)
    @test_throws ArgumentError LS.SearchOptions(minimum_width_power=-1)
    @test_throws ArgumentError LS.SearchOptions(krawczyk_depth=-1)
    for n in (5, 6)
        data = D.read_json_file(proof_input_path(n))
        target = string(D._announced_target_entry(n)["graph6"])
        families = Set(string(e["graph6"]) for e in data["graphs"]
            if get(e, "representation", "") == "exact_family")
        codes = D.ProofWorkflow.regular_search_graphs(load_connected_graphs(n), target, families)
        @test length(codes) == (n == 5 ? 12 : 44)
        @test count(code -> !isempty(CandidateSearch.LeafFilters.graph_data(code).leaves), codes) ==
            (n == 5 ? 9 : 33)
    end
    mktempdir() do output
        @test_throws "unresolved graphs" D.ProofWorkflow.execute([
            "--n", "3", "--max-boxes", "1", "--precision", "128", "--output", output])
        @test !isfile(joinpath(output, "n3", "run.json"))
    end
end

@testset "Leaf-chart root isolation and Cartesian transfer" begin
    D = DiskSecondMomentProof
    LS = CandidateSearch.IntervalSearch
    LE = CandidateSearch.LeafFilters
    MCS = D.MinimumCandidateSearch
    for precision in (128, 512)
        T = Core.apply_type(ArbReal, precision)
        graph = LE.graph_data("BW")
        # Use the other leaf: reconstructing the default tree would swap endpoints.
        chart = LS.angular_chart(graph; leaf=last(graph.leaves))
        @test chart.tree != LS.rooted_tree(graph)
        Q = Rational{BigInt}
        box = LE.AngleBox(Q[-1//100], Q[1//100], 12)
        status, root = MCS.isolate_candidate_box(T, graph, chart, [], box, T[1, 1])
        @test status == :isolated
        @test root.kind == :flexible_critical_root
        @test root.chart === chart
        refined = Refinement.refine_candidate("BW", root; precision=precision)
        @test refined.chart === chart
        selected = chart.selected_leaf
        other = first(graph.leaves)
        @test upperbound(refined.center_box.centers[1, selected]) < -1
        @test lowerbound(refined.center_box.centers[1, other]) > 1
        @test lowerbound(refined.moment_box) <= 4 <= upperbound(refined.moment_box)
        analysis = verify_refined_candidate(refined, "leaf-root"; precision=precision)
        @test analysis.rank_status == :full_row_rank_proved
        @test analysis.conclusion == :excluded_by_negative_intrinsic_direction
        # The normal filter sequence rejects this nonminimum before isolation.
        reason, candidate = MCS.classify_candidate_box(T, graph, chart, [], box, LS.SearchOptions())
        @test reason == :leaf_equilibrium
        @test candidate === nothing
        mktempdir() do output
            path = joinpath(output, "roots.jsonl")
            Refinement.write_refinement_report([("leaf-root", refined)], path)
            record = D.read_json_file(path)
            @test Tuple.(record["oriented_tree"]) == chart.tree
            @test record["selected_leaf"] == chart.selected_leaf
        end
    end
end

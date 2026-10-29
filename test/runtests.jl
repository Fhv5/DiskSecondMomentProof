using Test
using ArbNumerics
using DiskSecondMomentProof

@testset "Interval safety and serialized enclosures" begin
    safety = DiskSecondMomentProof.IntervalSafety
    for precision in (128, 512)
        T = Core.apply_type(ArbReal, precision)
        for value in (T(1) / T(3), -T(17) / T(7), setinterval(T(-1), T(2)))
            saved = safety.interval_record(value)
            restored = safety.restore_interval(T, saved)
            @test lowerbound(restored) <= lowerbound(value)
            @test upperbound(value) <= upperbound(restored)
        end
        @test_throws ArgumentError safety.interval_record(T("NaN"))
        @test_throws ArgumentError safety.interval_record(T("Inf"))
        @test !safety.zero_excluded(T("NaN"))
        @test !DiskSecondMomentProof.CandidateSearch.LeafFilters.strictly_outside_zero(T("NaN"))
        @test safety.verified_inverse(T[0 1; 1 0]) !== nothing
        @test safety.verified_inverse(T[1 1; 1 1]) === nothing
        @test safety.verified_inverse(T[1 0; 0 T("NaN")]) === nothing
    end
    @test_throws ArgumentError safety.restore_interval(ArbReal{128}, ["1", "0", "2", "1"])
end

@testset "Rational-to-Arb enclosure" begin
    safety = DiskSecondMomentProof.IntervalSafety
    leaf = DiskSecondMomentProof.CandidateSearch.LeafFilters
    interval = DiskSecondMomentProof.CandidateSearch.IntervalSearch
    Q = Rational{BigInt}
    near_one = Q(big(2)^600 + 1, big(2)^600)
    tiny_negative = -Q(big(2)^203 + 1, big(2)^707)
    values = Q[0, 1//2, 1//3, -22//7, 22//7,
               near_one, -near_one, tiny_negative]
    for precision in (24, 128, 512)
        T = Core.apply_type(ArbReal, precision)
        for value in values
            enclosure = safety.enclosing_rational(T, value)
            @test safety.rational_bound(lowerbound(enclosure, ArbFloat)) <= value
            @test value <= safety.rational_bound(upperbound(enclosure, ArbFloat))
        end
        for depth in (0, 24, 80, 140)
            lower = Q(-22, 7) + Q(1, big(2)^depth)
            upper = Q(22, 7) - Q(1, big(2)^(depth + 1))
            box = leaf.AngleBox(Q[lower], Q[upper], depth)
            for enclosure in (only(leaf.interval_angles(T, box)),
                              only(interval.interval_angles(T, box)))
                @test safety.rational_bound(lowerbound(enclosure, ArbFloat)) <= lower
                @test upper <= safety.rational_bound(upperbound(enclosure, ArbFloat))
            end
        end
        @test_throws ArgumentError safety.enclosing_rational_interval(T, 1//2, 1//3)
    end
end

@testset "Box transcript format" begin
    leaf = DiskSecondMomentProof.CandidateSearch.LeafFilters
    records = DiskSecondMomentProof.BoxRecords
    box = leaf.AngleBox(Rational{BigInt}[-22//7, 0],
                            Rational{BigInt}[22//7, 22//7], 0)
    io = IOBuffer()
    records.record_box(io, box, :split)
    @test String(take!(io)) == "0\tsplit\n"
end

@testset "Active workflow modules" begin
    @test_throws ErrorException DiskSecondMomentProof.ProofWorkflow.execute(["--n", "8"])
    @test_throws ArgumentError second_moment_target(8)
    leaf_result = DiskSecondMomentProof.MinimumCandidateSearch.search_candidate_graph(
        "DDW", CandidateSearch.IntervalSearch.SearchOptions(max_boxes=1),
    )
    @test leaf_result.status == :unresolved
    @test isempty(leaf_result.candidates)
    result = DiskSecondMomentProof.MinimumCandidateSearch.search_candidate_graph(
        "D`{", CandidateSearch.IntervalSearch.SearchOptions(max_boxes=1),
    )
    @test result.status == :unresolved
    @test result.stats.visited == 1
    @test !hasproperty(result.stats, :objective)
end

@testset "Census" begin
    @test load_connected_graphs(3) == ["BW", "Bw"]
    @test length(load_connected_graphs(5)) == 13
    @test length(load_connected_graphs(6)) == 46
    @test_throws ArgumentError load_connected_graphs(8)
    for n in (2, 3, 4, 5, 6)
        codes = load_connected_graphs(n)
        @test all(code -> Census.is_outerplanar(ContactGraph(code)), codes)
        @test Census.compare_with_stored_census(codes, n)
        planar = readlines(joinpath(DiskSecondMomentProof.package_root(),
            "data", "planar_graphs", "planar$n.g6"))
        @test length(planar) == Dict(2 => 1, 3 => 2, 4 => 6, 5 => 20, 6 => 99)[n]
        @test Census.compare_with_stored_census(
            Census.filter_planar_graphs(planar, n), n)
    end
    k4_with_tail = ContactGraph(5, vcat(
        [ContactEdge(i, j) for i in 1:3 for j in (i + 1):4],
        [ContactEdge(1, 5)],
    ))
    @test !Census.is_outerplanar(k4_with_tail)
    original = ContactGraph("DqK")
    permutation = [3, 2, 1, 4, 5]
    relabeled = ContactGraph(5, [
        ContactEdge(permutation[edge.src], permutation[edge.dst])
        for edge in original.edges
    ])
    @test Census._canonical_signature(original) == Census._canonical_signature(relabeled)
end

@testset "Census partition" begin
    expected_leaf_counts = Dict(2 => 0, 3 => 1, 4 => 3, 5 => 9, 6 => 33)
    for n in (2, 3, 4, 5, 6)
        codes = load_connected_graphs(n)
        leaf_codes = filter(codes) do code
            graph = ContactGraph(code)
            n != 2 && any(vertex -> degree(graph, vertex) == 1, 1:n)
        end
        data = DiskSecondMomentProof.read_json_file(proof_input_path(n))
        target = only(filter(
            entry -> get(entry, "representation", "decimal") == "exact_qsqrt3",
            data["graphs"],
        ))
        family_codes = Set(string(entry["graph6"]) for entry in data["graphs"]
                           if get(entry, "representation", "") == "exact_family")
        selected = DiskSecondMomentProof.ProofWorkflow.regular_search_graphs(
            codes, string(target["graph6"]), family_codes)
        partition = vcat([string(target["graph6"])],
                         collect(family_codes), selected)
        @test length(leaf_codes) == expected_leaf_counts[n]
        @test Set(partition) == Set(codes)
        @test length(partition) == length(codes)
    end
end

@testset "Leaf symmetry coverage guard" begin
    leaf_search = CandidateSearch.LeafFilters
    for (n, expected_ordered) in ((5, 6), (6, 18))
        ordered = 0
        for code in load_connected_graphs(n)
            graph = leaf_search.graph_data(code)
            isempty(graph.leaves) && continue
            tree = leaf_search.anchored_tree(graph, first(graph.leaves))
            _, constraints = leaf_search.angular_symmetry_data(graph, tree)
            if haskey(leaf_search.VERIFIED_SYMMETRY_PATTERNS, code)
                @test !isempty(constraints)
                ordered += 1
            else
                @test isempty(constraints)
            end
        end
        @test ordered == expected_ordered
    end

    graph = leaf_search.graph_data("DBw")
    tree = leaf_search.anchored_tree(graph, first(graph.leaves))
    reflection_index, constraints = leaf_search.angular_symmetry_data(graph, tree)
    @test reflection_index == 3
    @test constraints == [(1, 2)]
    unregistered = leaf_search.GraphData(
        "unregistered", graph.n, graph.edges, graph.adjacency,
        graph.nonedges, graph.leaves, graph.bridge_cuts,
    )
    _, fallback_constraints = leaf_search.angular_symmetry_data(unregistered, tree)
    @test isempty(fallback_constraints)
    alternative_tree = [(1, 5), (5, 2), (5, 3), (3, 4)]
    _, altered_tree_constraints =
        leaf_search.angular_symmetry_data(graph, alternative_tree)
    @test isempty(altered_tree_constraints)
end

@testset "Contact normalization" begin
    centers = Rational{BigInt}[0 2; 0 0]
    configuration = Configuration(centers)
    graph = ContactGraph(2, [ContactEdge(1, 2)])
    @test contact_equations(configuration, graph) == [0]
    @test contact_matrix(configuration, graph) == Rational{BigInt}[-2 0 2 0]
    @test contact_jacobian(configuration, graph) == contact_matrix(configuration, graph)
    hessian = only(constraint_hessians(configuration, graph))
    @test hessian == Rational{BigInt}[
        1 0 -1 0;
        0 1 0 -1;
        -1 0 1 0;
        0 -1 0 1;
    ]
end

@testset "Candidate dataset integrity" begin
    for n in (2, 3, 4, 5, 6)
        data = DiskSecondMomentProof.read_json_file(proof_input_path(n))
        @test validate_candidate_dataset(data)
    end
end

@testset "Announced targets" begin
    for n in (2, 3, 5, 6)
        data = DiskSecondMomentProof.read_json_file(proof_input_path(n))
        target = only(filter(
            entry -> get(entry, "representation", "decimal") == "exact_qsqrt3",
            data["graphs"],
        ))
        result = verify_target(target)
        @test result.conclusion == :target_minimizer_exactly_verified
        @test result.objective == result.target
        @test result.objective == 2 * sum(result.multipliers)
        @test result.rolling_dimension_modulo_rigid_motions == 0
    end
end

@testset "Target geometry" begin
    for n in (2, 3, 4, 5, 6)
        result = CandidateSearch.TargetRealizations.analyze_target(proof_input_path(n))
        @test result.target_found
        @test length(result.representatives) == 1
        @test CandidateSearch.TargetRealizations.second_moment(only(result.representatives)) ==
            Qsqrt3(second_moment_target(n))
    end
end

@testset "Negative direction for the E`dg family" begin
    data = DiskSecondMomentProof.read_json_file(proof_input_path(6))
    family = only(filter(entry -> string(entry["graph6"]) == "E`dg", data["graphs"]))
    @test FamilyCoverage.verify_edg_coverage(family) == :complete_critical_family
    result = verify_family(family)
    @test result.model == :quadratic_arc
    @test result.rank_status == :full_row_rank_proved
    @test result.rank_witness_columns == [1, 2, 4, 5, 6, 8, 9]
    @test result.second_order_status == :exact_negative_direction
    @test result.negative_second_variation == -8 // 3
    @test result.conclusion == :excluded_by_exact_negative_intrinsic_direction
    zero_rank_witness = deepcopy(family)
    zero_rank_witness["family"]["rank_minor"]["columns"] = [1, 3, 5, 7, 9, 11, 2]
    zero_rank_witness["family"]["rank_minor"]["constant"] = "0"
    @test_throws "arc-family rank-minor constant must be nonzero" verify_family(zero_rank_witness)
    analysis, _ = DiskSecondMomentProof.analyze_search_results(6, [], [family])
    family_result = only(filter(item -> item.route == :exact_family, analysis.results))
    @test family_result.conclusion == :excluded_by_exact_negative_intrinsic_direction
end

@testset "EPVW boundary root" begin
    graph = CandidateSearch.LeafFilters.graph_data("EPVW")
    tree = CandidateSearch.IntervalSearch.rooted_tree(graph)
    witnesses = CandidateSearch.BoundaryRoots.registered_witnesses("EPVW")
    @test length(witnesses) == 1
    witness = only(witnesses)
    verified = CandidateSearch.BoundaryRoots.verify_witness(graph, tree, witness)
    @test verified.determinant == Qsqrt3(-4608)
    @test witness.angle_pi_multiples == [2 // 3, 2 // 3, 2 // 3, 1 // 3]
    @test witness.extra_contacts == [(1, 4)]
    @test witness.target_graph6 == "EQlw"
    @test witness.expected_moment == Qsqrt3(11)
end

@testset "Refinement and direct transfer of isolated roots" begin
    options = CandidateSearch.IntervalSearch.SearchOptions(max_boxes=500000)
    result = DiskSecondMomentProof.MinimumCandidateSearch.search_candidate_graph("D`{", options)
    @test result.status == :complete
    @test length(result.candidates) == 2
    for (index, candidate) in enumerate(result.candidates)
        refined = Refinement.refine_candidate("D`{", candidate)
        @test length(refined.lambda_box) == result.contacts
        @test length(refined.angle_box) == result.n - 2
        analysis = DiskSecondMomentProof.verify_refined_candidate(
            refined, "test-$index",
        )
        @test analysis.isolation_status == :angle_unique_transferred
        @test !hasproperty(analysis, :krawczyk_result)
        @test analysis.multiplier_status == :sign_undetermined
        @test DiskSecondMomentProof._candidate_table_fields(
            CandidateProofResult("test-$index", :transferred_angular_root,
                analysis.conclusion, analysis))[2] == :sign_undetermined
        @test analysis.hessian !== nothing
        @test analysis.hessian.status == :negative_direction
        @test analysis.objective_status == :strictly_above_target
        @test analysis.conclusion == :excluded_by_negative_intrinsic_direction
    end
end

@testset "Closure guard" begin
    target = CandidateProofResult(
        "target",
        :exact_target,
        :target_minimizer_exactly_verified,
        nothing,
    )

    excluded = CandidateProofResult(
        "candidate",
        :transferred_angular_root,
        :excluded_by_objective_bound,
        nothing,
    )

    complete = ComputationalProofResult(
        "test",
        6,
        2,
        2,
        [target, excluded],
        Dict(:exact_target => 1, :transferred_angular_root => 1),
        Dict(
            :target_minimizer_exactly_verified => 1,
            :excluded_by_objective_bound => 1,
        ),
    )

    @test isnothing(
        DiskSecondMomentProof.ProofWorkflow.assert_closure_complete(complete),
    )

    inconclusive = CandidateProofResult(
        "candidate",
        :transferred_angular_root,
        :special_treatment_required,
        nothing,
    )

    incomplete = ComputationalProofResult(
        "test",
        6,
        2,
        2,
        [target, inconclusive],
        Dict(:exact_target => 1, :transferred_angular_root => 1),
        Dict(
            :target_minimizer_exactly_verified => 1,
            :special_treatment_required => 1,
        ),
    )

    @test_throws ErrorException DiskSecondMomentProof.ProofWorkflow.assert_closure_complete(incomplete)
end

@testset "Three-disk complete workflow and evidence" begin
    data = DiskSecondMomentProof.read_json_file(proof_input_path(3))
    target = only(data["graphs"])
    target_result = verify_target(target)
    @test target_result.contact_rank == 3
    @test target_result.rolling_dimension_modulo_rigid_motions == 0
    @test target_result.minimum_noncontact_gap === nothing
    @test target_result.multipliers == fill(Qsqrt3(1//3), 3)
    @test target_result.objective == Qsqrt3(2)
    invalid = deepcopy(target)
    invalid["centros"]["1"][1]["a"] = "-2"
    invalid["centros"]["2"][1]["a"] = "2"
    @test_throws "contact equation failed" verify_target(invalid)

    for n in (2, 3, 5, 6)
        T = ArbReal{128}
        interval_target = DiskSecondMomentProof._target_value(T, n)
        @test DiskSecondMomentProof.IntervalSafety.rational_bound(lowerbound(interval_target, ArbFloat)) <= second_moment_target(n)
        @test second_moment_target(n) <= DiskSecondMomentProof.IntervalSafety.rational_bound(upperbound(interval_target, ArbFloat))
    end
    @test_throws ArgumentError DiskSecondMomentProof._target_value(ArbReal{128}, 8)

    mktempdir() do output
        result = DiskSecondMomentProof.ProofWorkflow.execute([
            "--n", "3", "--precision", "128", "--record-evidence", "--output", output,
        ])
        @test result.connected_census == 2
        @test length(result.search_results) == 1
        @test only(result.search_results).graph6 == "BW"
        @test only(result.search_results).status == :complete
        @test isempty(only(result.search_results).unresolved_boxes)
        @test isempty(only(result.search_results).candidates)
        @test result.analysis.processed_candidates == result.analysis.declared_candidates == 1
        @test only(result.analysis.results).conclusion == :target_minimizer_exactly_verified
        run = DiskSecondMomentProof.read_json_file(joinpath(output, "n3", "run.json"))
        @test run["schema_version"] == 3
        @test haskey(run["source_sha256"], ".gitattributes")
        @test haskey(run["source_text_sha256_lf"], "src/evidence/archive.jl")
        @test run["project_text_sha256_lf"] ==
            DiskSecondMomentProof.EvidenceArchive.sha256_text_lf(joinpath(dirname(@__DIR__), "Project.toml"))
        @test haskey(run["files_sha256"], "analysis/n3-interval-decisions.jsonl")
        @test run["connected_census"] == 2
        @test run["leaf_excluded"] == 1
        @test run["interval_graphs"] == 1
        @test run["isolated_boxes"] == 0
        @test run["analysis_objects"] == 1
        @test haskey(run["source_sha256"], "data/planar_graphs/planar3.g6")
        @test haskey(run["source_sha256"], "data/exact_inputs/3disks.json")
        @test haskey(run["files_sha256"], "refinement/candidates.jsonl")
        @test filesize(joinpath(output, "n3", "refinement", "candidates.jsonl")) == 0
        summary = readlines(joinpath(output, "n3", "search", "angle_boxes", "n3-summary.tsv"))
        @test length(summary) == 2
        report = read(joinpath(output, "n3", "analysis", "evidence", "n3", "001-C-n3-g6-4277-001.txt"), String)
        @test occursin("minimum_noncontact_gap = nothing", report)
    end
end

include("four_disks.jl")
include("two_disks.jl")
include("evidence_records.jl")

include("interval_minmax.jl")
include("full_equilibrium.jl")
include("search_conventions.jl")
include("public_api.jl")
include("common_search.jl")
include("execution_summary.jl")

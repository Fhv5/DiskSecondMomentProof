using Test
using ArbNumerics
using DiskSecondMomentProof

@testset "Two-disk target and reservation guard" begin
    Q = Qsqrt3
    @test second_moment_target(2) == 1//1
    @test load_connected_graphs(2) == ["A_"]
    @test Census.filter_planar_graphs(["A_"], 2) == ["A_"]
    @test_throws ArgumentError Census.filter_planar_graphs(["A?"], 2)
    @test_throws ArgumentError Census.filter_planar_graphs(["A_", "A_"], 2)
    data = DiskSecondMomentProof.read_json_file(proof_input_path(2))
    entry = only(data["graphs"])
    @test validate_candidate_dataset(data)
    target_result = verify_target(entry)
    @test target_result.contact_rank == 1
    @test target_result.rolling_dimension_modulo_rigid_motions == 0
    @test target_result.multipliers == [Q(1//2)]
    @test target_result.objective == target_result.target == Q(1)
    @test target_result.minimum_noncontact_gap === nothing
    audit = CandidateSearch.TargetRealizations.analyze_target(proof_input_path(2))
    @test audit.seed == (1, 2)
    @test isempty(audit.steps)
    @test audit.total_embeddings == 1
    @test audit.target_found && length(audit.representatives) == 1
    @test audit.classifications == [:exact_penny_realization]
    invalid = deepcopy(entry)
    invalid["centros"]["2"][1]["a"] = "2"
    @test_throws ArgumentError verify_target(invalid)
    invalid = deepcopy(entry)
    invalid["multipliers"]["1-2"] = "0"
    @test_throws ArgumentError verify_target(invalid)
    @test isempty(DiskSecondMomentProof.ProofWorkflow.regular_search_graphs(["A_"], "A_", Set{String}()))
    @test_throws ArgumentError DiskSecondMomentProof.MinimumCandidateSearch.search_candidate_graph(
        "A_", CandidateSearch.IntervalSearch.SearchOptions())

end

@testset "Two-disk complete workflow and evidence" begin
    mktempdir() do output
        result = DiskSecondMomentProof.ProofWorkflow.execute([
            "--n", "2", "--precision", "128", "--record-evidence", "--output", output,
        ])
        @test result.connected_census == 1
        @test isempty(result.search_results)
        @test result.analysis.processed_candidates == result.analysis.declared_candidates == 1
        @test only(result.analysis.results).conclusion == :target_minimizer_exactly_verified
        @test isnothing(DiskSecondMomentProof.ProofWorkflow.assert_closure_complete(result.analysis))
        incomplete = (n=2, processed_candidates=1, declared_candidates=1, results=typeof(result.analysis.results)())
        @test_throws ErrorException DiskSecondMomentProof.ProofWorkflow.assert_closure_complete(incomplete)
        run = DiskSecondMomentProof.read_json_file(joinpath(output, "n2", "run.json"))
        @test run["connected_census"] == 1
        @test run["leaf_graphs"] == run["leaf_excluded"] == 0
        @test run["interval_graphs"] == run["isolated_boxes"] == 0
        @test run["analysis_objects"] == 1
        @test haskey(run["source_sha256"], "data/planar_graphs/planar2.g6")
        @test haskey(run["source_sha256"], "docs/proof.md")
        @test filesize(joinpath(output, "n2", "refinement", "candidates.jsonl")) == 0
        @test length(readlines(joinpath(output, "n2", "search", "angle_boxes", "n2-summary.tsv"))) == 1
    end
end

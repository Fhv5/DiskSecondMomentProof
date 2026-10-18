using Test
using ArbNumerics
using DiskSecondMomentProof

@testset "Four-disk exhaustive census and planar obstruction" begin
    raw = readlines(joinpath(DiskSecondMomentProof.package_root(), "data", "planar_graphs", "planar4.g6"))
    @test raw == ["CF", "CN", "CR", "C^", "Cr", "C~"]
    @test load_connected_graphs(4) == raw[1:5]
    coverage = Census.verify_four_disk_planar_census(raw)
    @test coverage.status == :complete_connected_census
    @test coverage.connected_classes == 6
    @test coverage.k4_graph6 == "C~"
    @test coverage.gram_rank == 3 > coverage.ambient_dimension == 2
    @test !Census.is_outerplanar(ContactGraph("C~"))
    @test Census.filter_planar_graphs(raw, 4) == raw[1:5]
    @test_throws ArgumentError Census.verify_four_disk_planar_census(raw[1:5])
    @test_throws ArgumentError Census.verify_four_disk_planar_census(vcat(raw[1:5], [raw[1]]))
    @test_throws ArgumentError Census.verify_four_disk_planar_census(vcat(raw[1:5], ["C?"]))
end

@testset "Four-disk family and zero-multiplier endpoints" begin
    D = DiskSecondMomentProof
    data = D.read_json_file(proof_input_path(4))
    endpoint_entry, family_entry = data["graphs"]
    @test validate_candidate_dataset(data)
    @test verify_four_disk_rhombus_coverage(family_entry) == :complete_penny_family
    @test D.FamilyCoverage.verify_family_coverage(family_entry) == :complete_penny_family
    family = verify_family(family_entry)
    @test family.objective == family.target == 4//1
    @test family.model == :four_disk_rhombus
    @test family.noncontact_count == 2
    @test family.multipliers == fill(1//2, 4)
    @test family.rank_status == :full_row_rank_proved
    @test family.second_order_status == :exact_flat_intrinsic_direction
    endpoint = verify_four_disk_endpoint(endpoint_entry, family_entry)
    @test endpoint.objective == endpoint.target == Qsqrt3(4)
    @test endpoint.family_id == family.candidate_id
    @test endpoint.contact_rank == 5
    @test Set(keys(endpoint.boundary_to_endpoint_permutations)) == Set(["1", "3"])
    @test endpoint.minimum_noncontact_gap == Qsqrt3(8)
    @test count(iszero, endpoint.multipliers) == 1
    @test endpoint.multipliers[end] == Qsqrt3(0)
    @test_throws "multipliers are not strictly positive" verify_target(endpoint_entry)

    for (a, b) in ((6//5, 8//5), (8//5, 6//5))
        _, _, symbolic, _, _ = D._four_disk_rhombus_terms(family_entry)
        centers = [D._rhombus_evaluate(symbolic[row, i], Qsqrt3(a), Qsqrt3(b)) for row in 1:2, i in 1:4]
        @test D._qsqrt3_second_moment(centers) == Qsqrt3(4)
        @test D._four_disk_distance(centers, 1, 4) != Qsqrt3(4)
        @test D._four_disk_distance(centers, 2, 3) != Qsqrt3(4)
    end
    rational_text(value) = string(numerator(value), "/", denominator(value))
    for boundary in (1, 3)
        variant = deepcopy(endpoint_entry)
        variant["graph6"] = boundary == 1 ? "Cv" : "Cz"
        graph = ContactGraph(variant["graph6"])
        centers = D._four_disk_boundary_centers(family_entry, boundary)
        variant["contactos"] = [[edge.src, edge.dst] for edge in graph.edges]
        variant["centros"] = Dict{String,Any}(string(i) => [Dict{String,Any}("a" => rational_text(centers[row, i].a),
            "b" => rational_text(centers[row, i].b)) for row in 1:2] for i in 1:4)
        variant["multipliers"] = Dict{String,Any}("$(edge.src)-$(edge.dst)" =>
            ((edge.src, edge.dst) in D._RHOMBUS_EDGES ? "1/2" : "0") for edge in graph.edges)
        verified = verify_four_disk_endpoint(variant, family_entry)
        @test verified.objective == Qsqrt3(4)
        @test count(iszero, verified.multipliers) == 1
    end

    for mutation in (
        entry -> (entry["family"]["domain"]["a_squared_interval"]["upper"] = "4"),
        entry -> (entry["family"]["domain"]["a_squared_interval"]["lower_closed"] = true),
        entry -> (entry["family"]["domain"]["positive_axes"] = false),
        entry -> (entry["family"]["domain"]["relation"] = "b^2=5-a^2"),
        entry -> (entry["family"]["centers"]["1"][1]["p"] = ["0", "-2"]),
        entry -> (entry["family"]["multipliers"]["1-2"] = "1/3"),
        entry -> (entry["family"]["rank_minor_columns"] = [1, 2, 3, 4]),
        entry -> (entry["family"]["flat_direction"]["1"][1] = Dict{String,Any}("p" => ["0"], "q" => ["0"])),
        entry -> (entry["graph6"] = "C^"),
    )
        invalid = deepcopy(family_entry)
        mutation(invalid)
        @test_throws ArgumentError verify_family(invalid)
    end
    invalid_endpoint = deepcopy(endpoint_entry)
    invalid_endpoint["multipliers"]["3-4"] = "1/2"
    @test_throws ArgumentError verify_four_disk_endpoint(invalid_endpoint, family_entry)
end

@testset "Four-disk complete workflow, equality closure and evidence" begin
    D = DiskSecondMomentProof
    @test second_moment_target(4) == 4//1
    @test D._target_value(ArbReal{128}, 4) == ArbReal{128}(4)
    mktempdir() do output
        result = D.ProofWorkflow.execute([
            "--n", "4", "--precision", "128", "--record-evidence", "--output", output,
        ])
        @test result.leaf_filter.status == :proved_excluded
        @test result.leaf_filter.total_graphs == 5
        @test result.leaf_filter.leaf_graphs == result.leaf_filter.proved_graphs == 3
        @test Set(graph.graph6 for graph in result.leaf_filter.graph_results) == Set(["CF", "CN", "CR"])
        @test all(isempty(run.unresolved) for graph in result.leaf_filter.graph_results for run in graph.runs)
        @test isempty(result.search_results)
        @test result.analysis.processed_candidates == result.analysis.declared_candidates == 2
        @test result.analysis.conclusion_counts == Dict(:target_endpoint_exactly_verified => 1, :target_family_exactly_verified => 1)
        @test isnothing(D.ProofWorkflow.assert_closure_complete(result.analysis))
        endpoint, family = result.analysis.results
        proof_with(results) = ComputationalProofResult("test", 4, length(results), length(results),
            results, Dict{Symbol,Int}(), Dict{Symbol,Int}())
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([endpoint]))
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([family]))
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([endpoint, family, family]))
        bad_link = CandidateProofResult("wrong", family.route, family.conclusion, family.evidence)
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([endpoint, bad_link]))
        fake = CandidateProofResult(endpoint.candidate_id, endpoint.route, endpoint.conclusion, nothing)
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([fake, family]))
        partial_fields = [getfield(endpoint.evidence, name) for name in fieldnames(FourDiskEndpointResult)]
        partial_fields[10] = Dict("1" => endpoint.evidence.boundary_to_endpoint_permutations["1"])
        partial = CandidateProofResult(endpoint.candidate_id, endpoint.route, endpoint.conclusion,
            FourDiskEndpointResult(partial_fields...))
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([partial, family]))
        pending = CandidateProofResult("pending", :transferred_angular_root, :special_treatment_required, nothing)
        @test_throws ErrorException D.ProofWorkflow.assert_closure_complete(proof_with([endpoint, family, pending]))
        run = D.read_json_file(joinpath(output, "n4", "run.json"))
        @test run["connected_census"] == 5
        @test run["leaf_excluded"] == 3
        @test run["angular_graphs"] == run["angular_candidates"] == 0
        @test run["analysis_objects"] == 2
        @test run["equality_family_objects"] == run["equality_endpoint_objects"] == 1
        @test run["four_disk_scope"]["excluded_graph6"] == "C~"
        @test run["four_disk_scope"]["family_boundaries"] == [1, 3]
        @test haskey(run["source_sha256"], "docs/proof.md")
        @test filesize(joinpath(output, "n4", "refinement", "candidates.jsonl")) == 0
        summary = read(joinpath(output, "n4", "analysis", "n4-candidate-summary.tsv"), String)
        @test occursin("nonnegative_exact_with_zero_contact", summary)
        @test occursin("exact_flat_intrinsic_direction", summary)
        @test occursin("target_family_exactly_verified", summary)
    end
end

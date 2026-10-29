@testset "Execution summary derives counts and minima from results" begin
    D = DiskSecondMomentProof
    candidate(code, id, negative) = CandidateProofResult(id, :transferred_angular_root,
        :excluded_by_objective_bound, (graph6=code, objective_status=:strictly_above_target,
            hessian=negative ? (status=:negative_direction,) : nothing))
    graph(code, boxes) = (graph6=code, status=:complete, candidates=fill(nothing, boxes), unresolved_boxes=[])
    target = CandidateProofResult("target", :exact_target, :target_minimizer_exactly_verified,
        (graph6="target-g6", objective=Qsqrt3(17//5)))
    family = CandidateProofResult("family", :exact_family, :target_family_exactly_verified,
        (graph6="family-g6", objective=17//5))
    objects = [target, family, candidate("root-g6", "box-1", true), candidate("root-g6", "box-2", false)]
    analysis = (processed_candidates=length(objects), declared_candidates=length(objects), results=objects)
    searches = [graph("empty-g6", 0), graph("root-g6", 2)]
    census = ["empty-g6", "root-g6", "target-g6", "family-g6"]
    algebraic = ["target-g6", "family-g6"]
    data = D.execution_summary(analysis; census_codes=census, graph_results=searches,
        algebraic_codes=algebraic, complete=true)
    @test data.census == data.accounted == length(census)
    @test data.interval == data.completed == length(searches)
    @test data.with_candidates == data.without_candidates == 1
    @test data.retained == data.analysed == data.excluded == data.above_target == 2
    @test data.negative == 1
    @test data.minimising_graphs == sort(algebraic)
    @test data.minimum == Qsqrt3(17//5)
    buffer = IOBuffer()
    D.print_proof_summary(analysis; n=8, census_codes=census, graph_results=searches,
        algebraic_codes=algebraic, complete=true, io=buffer)
    output = String(take!(buffer))
    @test occursin("Execution for 8 disks", output)
    @test occursin(r"Graph classes accounted for\s+4 / 4", output)
    @test occursin(r"Minimising contact graphs \(graph6\)\s+family-g6, target-g6", output)
    @test occursin(r"Minimum second moment\s+17/5", output)
    @test occursin(r"Second moment strictly above 17/5\s+2", output)
    @test !occursin("Proof status", output) && !occursin("All minimisers", output)
    @test !occursin("Processed:", output) && !occursin("Routes:", output)

    pending = (graph6="root-g6", status=:unresolved, candidates=[nothing], unresolved_boxes=[nothing])
    incomplete = D.print_proof_summary(nothing; n=8, census_codes=census,
        graph_results=[graph("empty-g6", 0), pending], algebraic_codes=algebraic, io=buffer)
    output = String(take!(buffer))
    @test !incomplete.complete && incomplete.unresolved_boxes == 1
    @test occursin(r"Minimum second moment\s+not established", output)
    @test !occursin("17/5", output)
    @test_throws ArgumentError D.execution_summary(analysis; census_codes=census,
        graph_results=[graph("empty-g6", 0)], algebraic_codes=algebraic, complete=true)
    partial = D.print_proof_summary(analysis; n=8, census_codes=census,
        graph_results=searches, algebraic_codes=algebraic, complete=false, io=buffer)
    @test !partial.complete
    @test occursin(r"Minimum second moment\s+not established", String(take!(buffer)))
end

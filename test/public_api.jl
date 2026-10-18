@testset "Public API aliases" begin
    D = DiskSecondMomentProof
    @test D.TargetResult === D.ExactMinimizerResult
    @test D.FamilyResult === D.ExactFamilyResult
    @test D.verify_target === D.verify_exact_minimizer
    @test D.verify_family === D.verify_exact_family
    @test D.proof_input_path === D.exact_input_path
    @test D.FamilyCoverage === D.ExactFamilyCoverage
    @test D.IntervalSafety.rational_bound === D.IntervalSafety.exact_rational_bound
    R = D.CandidateSearch.TargetRealizations
    @test R.Point === R.ExactPoint
    @test R.second_moment === R.exact_objective

    # Evidence reports identify these types by name.
    @test nameof(D.TargetResult) == :ExactMinimizerResult
    @test nameof(D.FamilyResult) == :ExactFamilyResult
    @test nameof(R.Point) == :ExactPoint

    # Existing inputs and their terminal conclusions use the stored format.
    data = D.read_json_file(D.exact_input_path(3))
    entry = only(filter(item -> item["representation"] == "exact_qsqrt3", data["graphs"]))
    @test D.verify_exact_minimizer(entry).conclusion == :target_minimizer_exactly_verified
    @test D.verify_target(entry) isa D.ExactMinimizerResult
    family_data = D.read_json_file(D.exact_input_path(4))
    family_entry = only(filter(item -> item["representation"] == "exact_family", family_data["graphs"]))
    @test D.verify_exact_family(family_entry).conclusion == :target_family_exactly_verified
    @test D.verify_family(family_entry) isa D.ExactFamilyResult
end

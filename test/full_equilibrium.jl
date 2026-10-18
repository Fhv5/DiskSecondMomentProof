@testset "Rectangular interval equilibrium solver" begin
    solver = DiskSecondMomentProof.FullEquilibrium
    D = DiskSecondMomentProof
    LE = D.CandidateSearch.LeafExclusion
    LS = D.CandidateSearch.LeaflessSearch
    QB = Rational{BigInt}
    for precision in (128, 512)
        T = Core.apply_type(ArbReal, precision)
        # An overdetermined compatible system: eight equations, four unknowns.
        B = T[-1 0 0 0; 0 0 0 -1; 1 0 0 0; 0 -1 0 0;
              0 0 1 0; 0 1 0 0; 0 0 -1 0; 0 0 0 1]
        gradient = T[-1,-1,1,-1,1,1,-1,1]
        status, multipliers = solver.enclose_multipliers(T, permutedims(B), gradient)
        @test status == :ok
        @test all(v -> lowerbound(v) <= 1 <= upperbound(v), multipliers)

        status, _ = solver.enclose_multipliers(T, reshape(T[1,1], 1, 2), T[1,2])
        @test status == :stationarity

        status, multipliers = solver.enclose_multipliers(T, reshape(T[0,2], 1, 2), T[0,2])
        @test status == :ok
        @test lowerbound(only(multipliers)) <= 1 <= upperbound(only(multipliers))

        status, _ = solver.enclose_multipliers(T, T[1 2 3; 1 2 3], T[2,4,6])
        @test status == :rank

        uncertain = D.IntervalSafety.enclosing_rational_interval(T, -1//1, 1//1)
        status, _ = solver.enclose_multipliers(T, reshape(T[uncertain], 1, 1), T[1])
        @test status == :rank

        graph = LE.graph_data("Bw")
        tree = LS.rooted_tree(graph)
        box = LE.LeafAngleBox(QB[1], QB[11//10], 0)
        terms = LS.reduced_terms(T, graph, tree, LS.interval_angles(T, box))
        A = LE.contact_matrix_interval(T, terms.raw, graph)
        status, multipliers = solver.enclose_multipliers(T, A, vec(terms.centers))
        @test status == :ok
        @test all(v -> lowerbound(v) <= T(1)/3 <= upperbound(v), multipliers)

        # Even all-row interval elimination need not resolve a dependency that
        # ordinary enclosures have forgotten: lambda=2+cos(theta), lambda=3+cos(theta).
        cosine = D.IntervalSafety.enclosing_rational_interval(T, -1//1, 1//1)
        status, _ = solver.enclose_multipliers(T, reshape(T[1,1], 1, 2), T[2+cosine,3+cosine])
        @test status == :ok

        # A contradictory leftover row after two pivots must exclude the box.
        status, _ = solver.enclose_multipliers(T, T[1 0 1; 0 1 1], T[1,2,4])
        @test status == :stationarity
        # Negative solutions are enclosed here; the caller applies sign constraints.
        status, multipliers = solver.enclose_multipliers(T, reshape(T[1,2], 1, 2), T[-1,-2])
        @test status == :ok
        @test upperbound(only(multipliers)) < 0
        @test solver.enclose_multipliers(T, T[1 0; 0 1; 1 1], T[1,1])[1] == :rank
        @test solver.enclose_multipliers(T, reshape(T[NaN], 1, 1), T[1])[1] == :rank
        @test solver.enclose_multipliers(T, reshape(T[1], 1, 1), T[Inf])[1] == :rank
        @test_throws DimensionMismatch solver.enclose_multipliers(T, T[1 2], T[1])

        # Interval matrices covering whole families, including uncertainty in pivots.
        # Every sampled rational solution must lie in the returned enclosure.
        for unknowns in 1:4
            equations = unknowns + 2
            rational_matrix = [QB(i == j ? 8 * unknowns : (i + 2*j) % 3 - 1)
                            for i in 1:equations, j in 1:unknowns]
            radius = QB(1,100)
            coefficients = [D.IntervalSafety.enclosing_rational_interval(
                T, x-radius, x+radius) for x in rational_matrix]
            lambda_intervals = [D.IntervalSafety.enclosing_rational_interval(
                T, QB(j)-radius, QB(j)+radius) for j in 1:unknowns]
            gradient_intervals = coefficients * lambda_intervals
            original_coefficients = deepcopy(coefficients)
            original_gradient = deepcopy(gradient_intervals)
            status, multipliers = solver.enclose_multipliers(
                T, permutedims(coefficients), gradient_intervals)
            @test status == :ok
            @test coefficients == original_coefficients
            @test gradient_intervals == original_gradient
            for sample in -1:1
                rational_coefficients = rational_matrix .+ sample * radius
                rational_lambda = [QB(j) - sample * radius for j in 1:unknowns]
                rational_gradient = rational_coefficients * rational_lambda
                for i in 1:equations
                    @test D.IntervalSafety.rational_bound(lowerbound(gradient_intervals[i], ArbFloat)) <= rational_gradient[i] <=
                          D.IntervalSafety.rational_bound(upperbound(gradient_intervals[i], ArbFloat))
                end
                for j in 1:unknowns
                    @test D.IntervalSafety.rational_bound(lowerbound(multipliers[j], ArbFloat)) <= rational_lambda[j] <=
                          D.IntervalSafety.rational_bound(upperbound(multipliers[j], ArbFloat))
                end
            end
        end
    end
end

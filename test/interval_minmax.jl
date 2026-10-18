@testset "FLINT min/max precision and interval containment" begin
    safety = DiskSecondMomentProof.IntervalSafety
    cases = (
        ((1//1, 2//1), (3//1, 4//1)),
        ((-4//1, -2//1), (-3//1, 5//1)),
        ((1//3, 2//3), (1//2, 3//4)),
        ((0//1, 0//1), (-1//1, 1//1)),
        ((2//1, 2//1), (2//1, 2//1)),
    )
    for bits in (128, 512)
        T = Core.apply_type(ArbReal, bits)
        for ((a, b), (c, d)) in cases
            x = safety.enclosing_rational_interval(T, a, b)
            y = safety.enclosing_rational_interval(T, c, d)
            before_x, before_y = copy(x), copy(y)
            for (operation, rational_operation) in
                    ((safety.interval_min, min), (safety.interval_max, max))
                result = operation(x, y)
                @test safety.finite_interval(result)
                @test safety.rational_bound(lowerbound(result, ArbFloat)) <=
                      rational_operation(a, c)
                @test safety.rational_bound(upperbound(result, ArbFloat)) >=
                      rational_operation(b, d)
            end
            @test x == before_x && y == before_y
        end
        # Exercise repeated min/max calls under GC.
        for repetition in 1:200
            x = T(repetition) / T(3)
            y = T(repetition + 1) / T(3)
            lo = safety.interval_min(x, y)
            hi = safety.interval_max(x, y)
            @test lowerbound(lo) <= upperbound(x)
            @test upperbound(hi) >= lowerbound(y)
            repetition % 20 == 0 && GC.gc()
        end
    end
end

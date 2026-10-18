using SHA

@testset "Raw and LF-normalized source hashes" begin
    archive = DiskSecondMomentProof.EvidenceArchive
    mktempdir() do directory
        lf = joinpath(directory, "lf.jl")
        crlf = joinpath(directory, "crlf.jl")
        mixed = joinpath(directory, "mixed.jl")
        changed = joinpath(directory, "changed.jl")
        binary = joinpath(directory, "binary.bin")
        write(lf, "α=1\nβ=2\n")
        write(crlf, "α=1\r\nβ=2\r\n")
        write(mixed, "α=1\r\nβ=2\n")
        write(changed, "α=1\nβ=3\n")
        write(binary, UInt8[0xff,0x00,0x0d,0x0a])
        @test archive.sha256_file(lf) != archive.sha256_file(crlf)
        @test archive.sha256_text_lf(lf) == archive.sha256_text_lf(crlf) == archive.sha256_text_lf(mixed)
        @test archive.sha256_text_lf(lf) != archive.sha256_text_lf(changed)
        @test_throws ArgumentError archive.sha256_text_lf(binary)
        inventory = archive.text_file_hashes(directory)
        @test Set(keys(inventory)) == Set(["lf.jl","crlf.jl","mixed.jl","changed.jl"])
        @test inventory["lf.jl"] == bytes2hex(sha256(codeunits("α=1\nβ=2\n")))
    end
end

@testset "Interval decisions are exported with witnesses and root references" begin
    T = ArbReal{128}
    hessian = IntrinsicHessianVerification(
        :negative_direction, 4, 1, [1], [1,2,3], zeros(T,4,4),
        reshape(T[-2],1,1), T[1,0], T(-2),
    )
    moment = setinterval(T(7), T(71)/T(10))
    regular = RegularCandidateResult(
        "sample-root", "sample", "D`{", 5, 6, 128,
        :angle_unique_transferred, :exact_graph_proved, T(8),
        :full_row_rank_proved, 6, :original_balance_proved,
        :not_excluded, nothing, nothing, :sign_undetermined,
        T[0,1], hessian, moment, :strictly_above_target,
        :excluded_by_negative_intrinsic_direction, String[],
    )
    candidate = CandidateProofResult("sample-root", :transferred_angular_root,
                                     regular.conclusion, regular)
    result = ComputationalProofResult("fixture",5,1,1,[candidate],
              Dict(:transferred_angular_root=>1), Dict(regular.conclusion=>1))
    safety = DiskSecondMomentProof.IntervalSafety
    mktempdir() do directory
        write_candidate_reports(result,directory)
        data = DiskSecondMomentProof.read_json_file(
            joinpath(directory,"n5-interval-decisions.jsonl"))
        @test data["candidate_id"] == "sample-root"
        @test data["target_rational"] == "34/5"
        @test data["selected_conclusion"] == "excluded_by_negative_intrinsic_direction"
        @test data["root_enclosure_candidate_id"] == "sample-root"
        @test data["root_enclosure_file"] == "../refinement/candidates.jsonl"
        restored = safety.restore_interval(T,data["objective_interval"])
        @test lowerbound(restored) <= lowerbound(moment)
        @test upperbound(moment) <= upperbound(restored)
        margin = safety.restore_interval(T,data["objective_margin_interval"])
        @test lowerbound(margin) > 0
        negative = safety.restore_interval(T,data["negative_intrinsic_value"])
        @test upperbound(negative) < 0
        @test length(data["negative_intrinsic_direction"]) == 2
        @test data["multiplier_status"] == "sign_undetermined"
        @test length(data["multipliers"]) == 2
    end
end

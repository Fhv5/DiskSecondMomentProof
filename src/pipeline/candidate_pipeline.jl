"""Route, conclusion, and mathematical evidence for one candidate."""
struct CandidateProofResult
    candidate_id::String
    route::Symbol
    conclusion::Symbol
    evidence::Any
end

"""Aggregate the analysis results for one disk count."""
struct ComputationalProofResult
    filepath::String
    n::Int
    declared_candidates::Int
    processed_candidates::Int
    results::Vector{CandidateProofResult}
    route_counts::Dict{Symbol,Int}
    conclusion_counts::Dict{Symbol,Int}
end

function _announced_target_entry(n::Int)
    data = read_json_file(proof_input_path(n))
    targets = filter(
        entry -> get(entry, "representation", "decimal") == "exact_qsqrt3",
        data["graphs"],
    )
    length(targets) == 1 || error("expected exactly one target for n=$n")
    return only(targets)::Dict{String,Any}
end

"""Analyze the roots transferred directly from the angular search."""
function analyze_search_results(n::Int, graph_results, family_entries; precision=512)
    results = CandidateProofResult[]
    refined_records = Tuple{String,Any}[]
    target_entry = _announced_target_entry(n)
    if n == 4
        family_entry = only(filter(entry -> entry["family"]["model"] == "four_disk_rhombus", family_entries))
        target_evidence = verify_four_disk_endpoint(target_entry, family_entry)
        target_route = :exact_target_endpoint
    else
        target_evidence = verify_target(target_entry)
        target_route = :exact_target
    end
    push!(results, CandidateProofResult(
        string(target_entry["critical_id"]), target_route,
        target_evidence.conclusion, target_evidence,
    ))

    for entry in family_entries
        evidence = verify_family(entry)
        push!(results, CandidateProofResult(
            string(entry["critical_id"]), :exact_family,
            evidence.conclusion, evidence,
        ))
    end

    for graph_result in graph_results
        for (index, candidate) in enumerate(graph_result.candidates)
            identifier = "S-n$(n)-g6-$(uppercase(bytes2hex(codeunits(graph_result.graph6))))-$index"
            refined = Refinement.refine_candidate(
                graph_result.graph6, candidate; precision=precision,
            )
            push!(refined_records, (identifier, refined))
            evidence = verify_refined_candidate(refined, identifier; precision=precision)
            push!(results, CandidateProofResult(
                identifier, :transferred_angular_root,
                evidence.conclusion, evidence,
            ))
        end
    end

    route_counts = Dict{Symbol,Int}()
    conclusion_counts = Dict{Symbol,Int}()
    for result in results
        route_counts[result.route] = get(route_counts, result.route, 0) + 1
        conclusion_counts[result.conclusion] = get(conclusion_counts, result.conclusion, 0) + 1
    end
    analysis = ComputationalProofResult(
        "Angular search and direct interval transfer", n,
        length(results), length(results), results, route_counts, conclusion_counts,
    )
    return analysis, refined_records
end

"""Convert a candidate identifier into a portable report filename."""
function _report_filename(candidate_id::AbstractString)
    safe = replace(String(candidate_id), r"[^A-Za-z0-9._-]" => "_")
    return isempty(safe) ? "candidate" : safe
end

"""Summarize proof statuses without rounding an interval to an apparent zero."""
function _candidate_table_fields(candidate::CandidateProofResult)
    evidence = candidate.evidence
    if evidence isa RegularCandidateResult
        hessian_status = evidence.hessian === nothing ? :not_computed : evidence.hessian.status
        return (evidence.graph6, evidence.multiplier_status, hessian_status,
            evidence.objective_status)
    elseif evidence isa TargetResult
        return (evidence.graph6, :strictly_positive_exact, :rolling_rigid,
            :equal_to_target_exact)
    elseif evidence isa FourDiskEndpointResult
        return (evidence.graph6, :nonnegative_exact_with_zero_contact,
            :exact_equality_family_endpoint, :equal_to_target_exact)
    elseif evidence isa FamilyResult
        multiplier_status = all(value -> value > 0, evidence.multipliers) ?
            :strictly_positive_exact :
            all(value -> value >= 0, evidence.multipliers) ?
                :nonnegative_exact : :negative_present_exact
        return (evidence.graph6, multiplier_status, evidence.second_order_status,
            evidence.objective > evidence.target ? :strictly_above_target_exact :
            evidence.objective == evidence.target ? :equal_to_target_exact :
            :strictly_below_target_exact)
    end
    throw(ArgumentError("unsupported candidate evidence type: $(typeof(evidence))"))
end

"""Return an explicit endpoint representation for interval proof data."""
function _evidence_scalar_text(value)
    value isa ArbReal && return "[$(lowerbound(value)), $(upperbound(value))]"
    return repr(value)
end

"""Write a recursively structured plain-text representation of proof evidence."""
function _write_evidence(io::IO, value, indentation::Int=0)
    prefix = repeat(" ", indentation)
    if value === nothing || value isa Number || value isa Symbol || value isa AbstractString || value isa Bool
        println(io, prefix, _evidence_scalar_text(value))
    elseif value isa AbstractArray{<:ArbReal}
        println(io, prefix, "interval_array size=", size(value))
        for index in CartesianIndices(value)
            println(
                io,
                prefix,
                "  ",
                Tuple(index),
                " = ",
                _evidence_scalar_text(value[index]),
            )
        end
    elseif value isa AbstractArray
        println(io, prefix, repr(value))
    elseif isstructtype(typeof(value))
        println(io, prefix, nameof(typeof(value)))
        for field in fieldnames(typeof(value))
            print(io, prefix, "  ", field, " = ")
            field_value = getfield(value, field)
            if field_value isa AbstractArray || field_value === nothing ||
                    field_value isa Number || field_value isa Symbol ||
                    field_value isa AbstractString || field_value isa Bool
                if field_value isa AbstractArray{<:ArbReal}
                    println(io)
                    _write_evidence(io, field_value, indentation + 4)
                else
                    println(io, _evidence_scalar_text(field_value))
                end
            else
                println(io)
                _write_evidence(io, field_value, indentation + 4)
            end
        end
    else
        println(io, prefix, repr(value))
    end
end

"""
    write_candidate_reports(result, output_directory)

Write a tabular index, plain-text evidence, and rational interval decision
records for regular candidates. Report generation has no role in any mathematical decision.
"""
function write_candidate_reports(
    result::ComputationalProofResult,
    output_directory::AbstractString,
)
    mkpath(output_directory)
    summary_path = joinpath(output_directory, "n$(result.n)-candidate-summary.tsv")
    evidence_directory = joinpath(output_directory, "evidence", "n$(result.n)")
    mkpath(evidence_directory)
    open(summary_path, "w") do io
        println(io, "candidate_id\tgraph6\troute\tmultiplier_status\thessian_status\tobjective_status\tconclusion\tevidence_file")
        for (index, candidate) in enumerate(result.results)
            filename = lpad(string(index), 3, '0') * "-" *
                _report_filename(candidate.candidate_id) * ".txt"
            path = joinpath(evidence_directory, filename)
            open(path, "w") do evidence_io
                println(evidence_io, "candidate_id = ", repr(candidate.candidate_id))
                println(evidence_io, "route = ", repr(candidate.route))
                println(evidence_io, "conclusion = ", repr(candidate.conclusion))
                println(evidence_io, "evidence =")
                _write_evidence(evidence_io, candidate.evidence, 2)
            end
            relative_path = joinpath("evidence", "n$(result.n)", filename)
            graph6, multiplier_status, hessian_status, objective_status =
                _candidate_table_fields(candidate)
            println(
                io,
                join((candidate.candidate_id, graph6, candidate.route,
                    multiplier_status, hessian_status, objective_status,
                    candidate.conclusion, relative_path), '\t'),
            )
        end
    end
    decision_path = joinpath(output_directory, "n$(result.n)-interval-decisions.jsonl")
    open(decision_path, "w") do io
        for (index, candidate) in enumerate(result.results)
            evidence = candidate.evidence
            evidence isa RegularCandidateResult || continue
            T = eltype(evidence.lambda_box)
            moment = evidence.objective_interval
            margin = moment === nothing ? nothing : moment - _target_value(T, result.n)
            curvature = evidence.hessian === nothing ? nothing : evidence.hessian.negative_value
            direction = evidence.hessian === nothing ? nothing :
                IntervalSafety.interval_record.(evidence.hessian.negative_direction)
            target = second_moment_target(result.n)
            encode(value) = value === nothing ? nothing : IntervalSafety.interval_record(value)
            filename = lpad(string(index), 3, '0') * "-" *
                _report_filename(candidate.candidate_id) * ".txt"
            JSON3.write(io, Dict(
                "schema_version" => 1,
                "candidate_id" => candidate.candidate_id,
                "graph6" => evidence.graph6,
                "n" => result.n,
                "selected_conclusion" => string(candidate.conclusion),
                "objective_status" => string(evidence.objective_status),
                "target_rational" => string(numerator(target)) * "/" * string(denominator(target)),
                "objective_interval" => encode(moment),
                "objective_margin_interval" => encode(margin),
                "global_exclusion_condition" => "inf(M-U_n)>0",
                "negative_intrinsic_value" => encode(curvature),
                "negative_intrinsic_direction" => direction,
                "local_exclusion_condition" => "sup(v^T H_intrinsic v)<0",
                "minimum_noncontact_gap" => encode(evidence.minimum_noncontact_gap),
                "rank_status" => string(evidence.rank_status),
                "multiplier_status" => string(evidence.multiplier_status),
                "multipliers" => IntervalSafety.interval_record.(evidence.lambda_box),
                "root_enclosure_file" => "../refinement/candidates.jsonl",
                "root_enclosure_candidate_id" => candidate.candidate_id,
                "analysis_evidence_file" => replace(joinpath("evidence", "n$(result.n)", filename), '\\' => '/'),
                "endpoint_encoding" => "lower numerator, lower denominator, upper numerator, upper denominator",
            ))
            println(io)
        end
    end
    return summary_path, evidence_directory
end

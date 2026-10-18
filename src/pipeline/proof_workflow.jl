module ProofWorkflow

using ArbNumerics
using ..DiskSecondMomentProof
const LS = DiskSecondMomentProof.CandidateSearch.LeaflessSearch
const LE = DiskSecondMomentProof.CandidateSearch.LeafExclusion
const BR = DiskSecondMomentProof.CandidateSearch.BoundaryRoots
const CD = DiskSecondMomentProof.CandidateDeduplication
const OFLE = DiskSecondMomentProof.ObjectiveFreeLeafExclusion
using ..MinimumCandidateSearch: search_candidate_graph

const TERMINAL_NON_TARGET_CONCLUSIONS = Set([
    :excluded_by_leaf_condition,
    :excluded_by_first_order,
    :excluded_by_negative_intrinsic_direction,
    :excluded_by_objective_bound,
    :excluded_by_exact_negative_intrinsic_direction,
    :excluded_by_exact_whole_family_objective,
])

"""Read an optional name/value argument with a supplied default."""
function option_value(arguments, name, default)
    index = findfirst(==(name), arguments)
    index === nothing && return default
    index < length(arguments) || error("missing value after $name")
    return arguments[index + 1]
end

"""Check whether a boolean command-line flag is present."""
has_flag(arguments, name) = name in arguments

"""Persist rational bounds for every unresolved search box."""
function write_unresolved_box_report(results, n, output_directory)
    path = joinpath(output_directory, "n$n-unresolved-boxes.tsv")
    open(path, "w") do io
        println(io, "graph6\tbox\tdepth\tlower\tupper")
        for result in results
            for (index, box) in enumerate(result.unresolved_boxes)
                println(
                    io,
                    join((
                        result.graph6,
                        index,
                        box.depth,
                        join(string.(box.lower), ','),
                        join(string.(box.upper), ','),
                    ), '\t'),
                )
            end
        end
    end
    return path
end

"""
Require analysis to give a terminal mathematical conclusion for every object
used by the complete proof.

The two-, three-, five- and six-disk targets must be verified in Q(sqrt(3)). The
four-disk case requires a linked equality family and its endpoint.
Every other analysis object must have a terminal exclusion conclusion. Inconclusive, diagnostic, inconsistent,
or target-overlap statuses are fatal for a complete proof run.
"""
function assert_closure_complete(analysis)
    analysis.processed_candidates == analysis.declared_candidates || error(
        "Analysis processed $(analysis.processed_candidates) of " *
        "$(analysis.declared_candidates) declared objects",
    )

    if analysis.n == 4
        equality_ids = DiskSecondMomentProof.assert_four_disk_closure(analysis)
        all(result -> result.candidate_id in equality_ids ||
            result.conclusion in TERMINAL_NON_TARGET_CONCLUSIONS, analysis.results) ||
            error("Four-disk closure has nonterminal objects")
        length(unique(result.candidate_id for result in analysis.results)) == length(analysis.results) ||
            error("Four-disk closure has duplicate object identifiers")
        return nothing
    end

    target_results = [
        result for result in analysis.results
        if result.route == :exact_target
    ]

    length(target_results) == 1 || error(
        "Analysis must contain exactly one target; found " *
        "$(length(target_results))",
    )

    target = only(target_results)
    target.conclusion == :target_minimizer_exactly_verified || error(
        "Target verification did not terminate successfully: " *
        "$(target.candidate_id) => $(target.conclusion)",
    )

    nonterminal = [
        result for result in analysis.results
        if result.route != :exact_target &&
           !(result.conclusion in TERMINAL_NON_TARGET_CONCLUSIONS)
    ]

    if !isempty(nonterminal)
        details = join(
            (
                "$(result.candidate_id) => $(result.conclusion)"
                for result in nonterminal
            ),
            ", ",
        )
        error("Closure is incomplete; nonterminal conclusions: $details")
    end

    return nothing
end

"""Select regular search graphs after the leaf, target, and family routes."""
function regular_search_graphs(all_codes, target_graph, family_graphs)
    return [
        code for code in all_codes
        if isempty(LE.graph_data(code).leaves) && code != target_graph &&
           !(code in family_graphs)
    ]
end

"""Execute the two-, three-, four-, five- or six-disk proof workflow and write all stage reports."""
function execute(arguments)
    started = time()
    value_options = Set(("--n", "--precision", "--max-boxes", "--max-depth",
        "--minimum-width-power", "--krawczyk-depth", "--leaf-max-boxes",
        "--leaf-max-depth", "--leaf-minimum-width-power", "--graph6", "--output"))
    index = 1
    while index <= length(arguments)
        option = arguments[index]
        if option == "--skip-analysis" || option == "--record-evidence"
            index += 1
            continue
        end
        option in value_options || throw(ArgumentError("unknown option: $option"))
        index < length(arguments) || throw(ArgumentError("missing value after $option"))
        index += 2
    end
    n = parse(Int, option_value(arguments, "--n", "5"))
    precision = parse(Int, option_value(arguments, "--precision", "512"))
    max_boxes = parse(Int, option_value(arguments, "--max-boxes", "3000000"))
    max_depth = parse(Int, option_value(arguments, "--max-depth", "140"))
    minimum_width_power = parse(Int, option_value(arguments, "--minimum-width-power", "60"))
    krawczyk_depth = parse(Int, option_value(arguments, "--krawczyk-depth", "12"))
    leaf_max_boxes = parse(Int, option_value(arguments, "--leaf-max-boxes", "250000"))
    leaf_max_depth = parse(Int, option_value(arguments, "--leaf-max-depth", "80"))
    leaf_minimum_width_power = parse(
        Int,
        option_value(arguments, "--leaf-minimum-width-power", "24"),
    )
    graph_filter = option_value(arguments, "--graph6", "")
    run_analysis = !has_flag(arguments, "--skip-analysis")
    record_evidence = has_flag(arguments, "--record-evidence")
    options = LS.SearchOptions(
        precision=precision,
        max_boxes=max_boxes,
        max_depth=max_depth,
        minimum_width_power=minimum_width_power,
        krawczyk_depth=krawczyk_depth,
    )
    n in (2, 3, 4, 5, 6) || error(
        "the proof supports n=2, n=3, n=4, n=5 and n=6",
    )
    all_codes = DiskSecondMomentProof.load_connected_graphs(n)
    if n == 4
        raw = readlines(joinpath(DiskSecondMomentProof.package_root(), "data", "planar_graphs", "planar4.g6"))
        coverage = DiskSecondMomentProof.Census.verify_four_disk_planar_census(raw)
        DiskSecondMomentProof.Census.compare_with_stored_census(
            DiskSecondMomentProof.Census.filter_planar_graphs(raw, 4), 4) ||
            error("Four-disk filtered census differs from the stored classes")
        println("Four-disk census coverage: $(coverage.connected_classes) classes; " *
                "$(coverage.k4_graph6) excluded by Gram rank $(coverage.gram_rank)")
    end
    target_graph = string(DiskSecondMomentProof._announced_target_entry(n)["graph6"])
    target_geometry =
    DiskSecondMomentProof.CandidateSearch.TargetRealizations.analyze_target(
        DiskSecondMomentProof.proof_input_path(n),
    )
    target_geometry.target_found || error(
        "The announced target was not found among the penny realizations",
    )

    length(target_geometry.representatives) == 1 || error(
        "Target contact graph has $(length(target_geometry.representatives)) " *
        "noncongruent penny realizations",
    )

    println(
        n == 4 ? "Endpoint graph audit: one geometric penny-realization class" :
        "Target realization audit: one geometric penny-realization class",
    )
    candidate_data = DiskSecondMomentProof.read_json_file(
        DiskSecondMomentProof.proof_input_path(n),
    )
    DiskSecondMomentProof.validate_candidate_dataset(candidate_data)
    family_entries = Dict{String,Any}[]
    for raw_entry in candidate_data["graphs"]
        entry = raw_entry::Dict{String,Any}
        get(entry, "representation", "decimal") == "exact_family" || continue
        graph6 = string(entry["graph6"])
        isempty(LE.graph_data(graph6).leaves) || continue
        push!(family_entries, entry)
    end
    for entry in family_entries
        coverage = DiskSecondMomentProof.FamilyCoverage.verify_family_coverage(entry)
        coverage in (:complete_critical_family, :complete_penny_family) || error(
            "family candidate coverage remains unresolved",
        )
    end
    family_graphs = Set(string(entry["graph6"]) for entry in family_entries)
    length(family_graphs) == length(family_entries) ||
        error("duplicate family graph entries")
    all(code -> code in all_codes, family_graphs) ||
        error("a family graph is absent from the census")
    codes = regular_search_graphs(all_codes, target_graph, family_graphs)
    if !isempty(graph_filter)
        graph_filter in codes || error(
            "graph6 $graph_filter is not a regular non-target search graph for n=$n",
        )
        codes = [graph_filter]
    end
    output_root = abspath(option_value(arguments, "--output", joinpath(DiskSecondMomentProof.package_root(), "output")))
    output_directory = isempty(graph_filter) ? joinpath(output_root, "n$n") :
        joinpath(output_root, "n$n", "single-$(bytes2hex(Vector{UInt8}(codeunits(graph_filter))))")
    leaf_output_directory = joinpath(output_directory, "search", "leaf_exclusion")
    search_output_directory = joinpath(output_directory, "search", "angle_boxes")
    evidence_directory = record_evidence ? joinpath(output_directory, "evidence") : nothing
    results = LS.GraphSearchResult[]

    println("Computational proof for n=$n")
    println("Connected census: $(length(all_codes)); precision: $precision bits")
    println("Candidate search: leaf exclusion running")
    println(
        "Leaf limits: max_boxes=$leaf_max_boxes max_depth=$leaf_max_depth " *
        "minimum_width=2^(-$leaf_minimum_width_power)",
    )
    leaf_result = OFLE.prove_leaf_census(
        all_codes;
        reserved_exact_targets=n == 2 ? [target_graph] : String[],
        options=LE.LeafExclusionOptions(
            precision=precision,
            max_boxes=leaf_max_boxes,
            max_depth=leaf_max_depth,
            minimum_width_power=leaf_minimum_width_power,
        ),
        verbose=false,
        trace_directory=evidence_directory,
    )
    LE.write_leaf_exclusion_reports(leaf_result, leaf_output_directory)
    println(
        "Leaf exclusion: $(leaf_result.proved_graphs)/$(leaf_result.leaf_graphs) " *
        "graphs proved excluded",
    )
    leaf_result.status == :proved_excluded || error(
        "leaf exclusion is incomplete; refusing to omit unresolved leaf graphs",
    )
    println(n == 4 ? "Equality family endpoint reserved for algebraic analysis: $target_graph" :
            "Announced target reserved for algebraic analysis: $target_graph")
    println(
        "Leafless families reserved for their analytic route: " *
        (isempty(family_graphs) ? "none" : join(sort!(collect(family_graphs)), ", ")),
    )
    println("Angular search graphs: $(length(codes))")
    println("Objective pruning: disabled; multiplier sign: nonnegative")
    for (index, code) in enumerate(codes)
        println("[$index/$(length(codes))] Starting graph $code")
        flush(stdout)
        result = if evidence_directory === nothing
            search_candidate_graph(code, options)
        else
            path = DiskSecondMomentProof.BoxRecords.box_record_path(
                evidence_directory, "angular", code)
            mkpath(dirname(path))
            open(path, "w") do io
                search_candidate_graph(code, options;
                    trace=(box, reason) ->
                        DiskSecondMomentProof.BoxRecords.record_box(io, box, reason))
            end
        end
        push!(results, result)
        println(
            "[$index/$(length(codes))] $code status=$(result.status) " *
            "boxes=$(result.stats.visited) candidates=$(length(result.candidates)) " *
            "unresolved=$(length(result.unresolved_boxes)) " *
            "time=$(round(result.elapsed_seconds; digits=2))s",
        )
    end

    LS.write_search_reports(results, n, search_output_directory)
    if record_evidence
        conventions_path = DiskSecondMomentProof.SearchConventions.write_search_conventions(
            output_directory, n, all_codes, leaf_result, results, target_graph,
            family_entries)
        println("Search conventions: $conventions_path")
    end
    unresolved_path = write_unresolved_box_report(results, n, search_output_directory)
    classes_by_result = CD.deduplicate_results(results)
    deduplication_path = CD.write_deduplication_report(
        results,
        classes_by_result,
        search_output_directory,
    )
    total_candidates = sum((length(result.candidates) for result in results); init=0)
    deduplicated_candidates = sum((length(classes) for classes in classes_by_result); init=0)
    unresolved_graphs = count(result -> result.status == :unresolved, results)
    println("Completed graphs: $(length(results) - unresolved_graphs)/$(length(results))")
    println("Raw isolated candidate boxes: $total_candidates")
    println("Diagnostic geometric classes: $deduplicated_candidates")
    println("Unresolved graphs: $unresolved_graphs")
    println("Search reports: $search_output_directory")
    println("Unresolved-box report: $unresolved_path")
    println("Deduplication report: $deduplication_path")

    unresolved_graphs == 0 || error(
        "Candidate search has unresolved graphs; see the written reports",
    )
    println("Candidate search completeness guard: passed")

    analysis = nothing
    if run_analysis
        println("Refinement and transfer: $total_candidates isolated angular roots")
        analysis, refined_records = DiskSecondMomentProof.analyze_search_results(
            n,
            results,
            family_entries;
            precision=precision,
        )
        refinement_path = DiskSecondMomentProof.Refinement.write_refinement_report(
            refined_records,
            joinpath(output_directory, "refinement", "candidates.jsonl"),
        )
        println("Refinement report: $refinement_path")
        DiskSecondMomentProof.print_proof_summary(analysis)
        report_directory = joinpath(output_directory, "analysis")
        DiskSecondMomentProof.write_candidate_reports(analysis, report_directory)
        println("Analysis reports: $report_directory")
        assert_closure_complete(analysis)
        if isempty(graph_filter)
            println("Closure: all required checks passed for the full census")
        else
            println(
                "Analysis terminal-status guard: passed for the selected graph; " *
                "this is a partial diagnostic run, not a complete global proof",
            )
        end
    end
    println("Total elapsed time: $(round(time() - started; digits=2))s")
    if !run_analysis
        println(
            "Analysis was skipped; this execution is a search diagnostic run, " *
            "not a complete global proof",
        )
    elseif isempty(graph_filter)
        println("Complete proof workflow: all required checks passed")
        if record_evidence
            parameters = Dict(
                "precision" => precision,
                "max_boxes" => max_boxes,
                "max_depth" => max_depth,
                "minimum_width_power" => minimum_width_power,
                "krawczyk_depth" => krawczyk_depth,
                "leaf_max_boxes" => leaf_max_boxes,
                "leaf_max_depth" => leaf_max_depth,
                "leaf_minimum_width_power" => leaf_minimum_width_power,
                "objective_pruning" => false,
            )
            record = DiskSecondMomentProof.EvidenceArchive.write_run_record(
                output_directory, n, parameters, leaf_result, results, analysis)
            println("Evidence inventory: $record")
        end
    end
    return (n=n, leaf_filter=leaf_result, search_results=results, analysis=analysis)
end

end

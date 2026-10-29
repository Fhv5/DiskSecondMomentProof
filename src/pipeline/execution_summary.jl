const MINIMISING_CONCLUSIONS = Set([
    :target_minimizer_exactly_verified,
    :target_family_exactly_verified,
    :target_endpoint_exactly_verified,
])

"""Compute reporting totals from search results and verified analysis objects."""
function execution_summary(analysis; census_codes, graph_results, algebraic_codes,
                           complete=false)
    census = Set(String.(census_codes))
    algebraic = Set(String.(algebraic_codes))
    interval = setdiff(census, algebraic)
    objects = analysis === nothing ? CandidateProofResult[] : analysis.results
    candidates = filter(item -> item.route == :transferred_angular_root, objects)
    algebraic_objects = filter(item -> item.route != :transferred_angular_root, objects)
    terminal = ProofWorkflow.TERMINAL_NON_TARGET_CONCLUSIONS
    retained = sum((length(result.candidates) for result in graph_results); init=0)
    excluded = count(item -> item.conclusion in terminal, candidates)
    completed = Set(result.graph6 for result in graph_results
        if result.status == :complete && isempty(result.unresolved_boxes))
    accounted = Set{String}()
    for result in graph_results
        result.graph6 in completed || continue
        items = filter(item -> item.evidence.graph6 == result.graph6, candidates)
        if length(items) == length(result.candidates) &&
           all(item -> item.conclusion in terminal, items)
            push!(accounted, result.graph6)
        end
    end
    for item in algebraic_objects
        if item.conclusion in union(terminal, MINIMISING_CONCLUSIONS)
            push!(accounted, item.evidence.graph6)
        end
    end
    minimisers = filter(item -> item.conclusion in MINIMISING_CONCLUSIONS, objects)
    values = [convert(Qsqrt3, item.evidence.objective) for item in minimisers]
    same_value = !isempty(values) && all(value -> value == first(values), values)
    ids = [item.candidate_id for item in objects]
    search_codes = [result.graph6 for result in graph_results]
    ready = analysis !== nothing && accounted == census &&
        Set(search_codes) == interval && length(search_codes) == length(interval) &&
        length(completed) == length(interval) && length(candidates) == retained == excluded &&
        analysis.processed_candidates == analysis.declared_candidates == length(objects) &&
        length(unique(ids)) == length(ids) && same_value
    complete && !ready && throw(ArgumentError("execution summary has incomplete census or analysis coverage"))
    return (
        census=length(census), interval=length(interval), algebraic=length(algebraic),
        completed=length(intersect(completed, interval)),
        without_candidates=count(result -> result.graph6 in completed && isempty(result.candidates), graph_results),
        with_candidates=count(result -> !isempty(result.candidates), graph_results),
        retained=retained,
        unresolved_boxes=sum((length(result.unresolved_boxes) for result in graph_results); init=0),
        analysed=length(candidates), excluded=excluded,
        above_target=count(item -> item.evidence.objective_status == :strictly_above_target, candidates),
        negative=count(item -> item.evidence.hessian !== nothing &&
            item.evidence.hessian.status == :negative_direction, candidates),
        unresolved_candidates=max(retained - excluded, 0),
        targets=filter(item -> item.route in (:exact_target, :exact_target_endpoint), algebraic_objects),
        families=filter(item -> item.route == :exact_family, algebraic_objects),
        accounted=length(intersect(accounted, census)),
        minimising_graphs=sort!(unique([item.evidence.graph6 for item in minimisers])),
        minimum=same_value ? first(values) : nothing,
        complete=complete && ready,
    )
end

_summary_number(value::Qsqrt3) = CandidateSearch.BoundaryRoots.qsqrt3_text(value)
_summary_row(io, label, value) = println(io, "  ", rpad(label, 49), lpad(string(value), 12))

"""Print census, search and analysis totals without treating boxes as graph classes."""
function print_proof_summary(analysis; n, census_codes, graph_results, algebraic_codes,
                             complete=false, elapsed_seconds=nothing,
                             evidence_directory=nothing, io::IO=stdout)
    summary = execution_summary(analysis; census_codes, graph_results, algebraic_codes, complete)
    row(label, value) = _summary_row(io, label, value)
    println(io, "\nExecution for $n disks\n")
    println(io, "GRAPH CENSUS")
    row("Connected outerplanar graph classes", summary.census)
    row("Assigned to interval search", summary.interval)
    row("Assigned to algebraic verification", summary.algebraic)
    println(io, "\nINTERVAL SEARCH")
    row("Graph searches completed", "$(summary.completed) / $(summary.interval)")
    row("Graphs excluded without retained candidates", summary.without_candidates)
    row("Graphs with retained candidates", summary.with_candidates)
    row("Candidate boxes retained", summary.retained)
    row("Unresolved boxes", summary.unresolved_boxes)
    println(io, "\nCANDIDATE ANALYSIS")
    row("Candidate boxes analysed", "$(summary.analysed) / $(summary.retained)")
    if summary.minimum !== nothing
        row("Second moment strictly above $(_summary_number(summary.minimum))", summary.above_target)
    end
    row("Negative intrinsic direction verified", summary.negative)
    row("Candidate boxes excluded", "$(summary.excluded) / $(summary.retained)")
    row("Unresolved candidates", summary.unresolved_candidates)
    println(io, "\nALGEBRAIC VERIFICATION")
    if analysis === nothing
        row("Algebraic analysis", "not performed")
    else
        row("Target graphs verified", "$(count(item -> item.conclusion in MINIMISING_CONCLUSIONS, summary.targets)) / $(length(summary.targets))")
    end
    critical = filter(item -> !(item.conclusion in MINIMISING_CONCLUSIONS), summary.families)
    minimising = filter(item -> item.conclusion in MINIMISING_CONCLUSIONS, summary.families)
    if !isempty(critical)
        row("Critical families excluded", "$(count(item -> item.conclusion in ProofWorkflow.TERMINAL_NON_TARGET_CONCLUSIONS, critical)) / $(length(critical))")
    end
    if !isempty(minimising)
        row("Minimising families verified", "$(length(minimising)) / $(length(minimising))")
    end
    println(io, "\nFINAL RESULT")
    row("Graph classes accounted for", "$(summary.accounted) / $(summary.census)")
    row("Minimising contact graphs (graph6)", summary.complete ? join(summary.minimising_graphs, ", ") : "not established")
    row("Minimum second moment", summary.complete ? _summary_number(summary.minimum) : "not established")
    if !summary.complete
        println(io, "\nGlobal conclusion not established by this execution.")
    end
    elapsed_seconds === nothing || println(io, "\nExecution time: $(round(elapsed_seconds; digits=2)) s")
    evidence_directory === nothing || println(io, "Evidence directory: $evidence_directory")
    return summary
end

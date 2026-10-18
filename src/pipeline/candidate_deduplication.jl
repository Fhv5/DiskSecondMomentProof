module CandidateDeduplication

using ArbNumerics
using ..DiskSecondMomentProof

const TR = DiskSecondMomentProof.CandidateSearch.TargetRealizations

export CandidateClass,
       deduplicate_result,
       deduplicate_results,
       write_deduplication_report

struct CandidateClass
    representative_index::Int
    member_indices::Vector{Int}
    matching_automorphisms::Vector{Vector{Int}}
end

function squared_distance_matrix(candidate)
    centers = candidate.center_box
    n = size(centers, 2)
    distances = Matrix{ArbReal}(undef, n, n)
    for i in 1:n
        distances[i, i] = zero(centers[1, i])
        for j in (i + 1):n
            dx = centers[1, i] - centers[1, j]
            dy = centers[2, i] - centers[2, j]
            value = dx * dx + dy * dy
            distances[i, j] = value
            distances[j, i] = value
        end
    end
    return distances
end

intervals_overlap(first, second) =
    DiskSecondMomentProof.IntervalSafety.interval_max(lowerbound(first), lowerbound(second)) <=
    DiskSecondMomentProof.IntervalSafety.interval_min(upperbound(first), upperbound(second))

function compatible_automorphism(first, second, automorphisms)
    n = size(first, 1)
    for permutation in automorphisms
        compatible = true
        for i in 1:(n - 1), j in (i + 1):n
            if !intervals_overlap(first[i, j], second[permutation[i], permutation[j]])
                compatible = false
                break
            end
        end
        compatible && return permutation
    end
    return nothing
end

function deduplicate_result(result)
    isempty(result.candidates) && return CandidateClass[]
    graph = ContactGraph(result.graph6)
    automorphisms = TR.graph_automorphisms(graph)
    distances = squared_distance_matrix.(result.candidates)
    classes = CandidateClass[]
    for candidate_index in eachindex(result.candidates)
        match = nothing
        matching_permutation = nothing
        for (class_index, candidate_class) in enumerate(classes)
            permutation = compatible_automorphism(
                distances[candidate_class.representative_index],
                distances[candidate_index],
                automorphisms,
            )
            if permutation !== nothing
                match = class_index
                matching_permutation = permutation
                break
            end
        end
        if match === nothing
            identity = collect(1:result.n)
            push!(classes, CandidateClass(candidate_index, [candidate_index], [identity]))
        else
            push!(classes[match].member_indices, candidate_index)
            push!(classes[match].matching_automorphisms, matching_permutation)
        end
    end
    return classes
end

deduplicate_results(results) = [deduplicate_result(result) for result in results]

function write_deduplication_report(results, classes_by_result, output_directory)
    mkpath(output_directory)
    path = joinpath(output_directory, "deduplicated-candidates.tsv")
    open(path, "w") do io
        println(io, "graph6\tclass\trepresentative\tmembers\tautomorphisms\tmethod")
        for (result, classes) in zip(results, classes_by_result)
            for (class_index, candidate_class) in enumerate(classes)
                members = join(candidate_class.member_indices, ',')
                automorphisms = join(
                    (join(permutation, ',') for permutation in candidate_class.matching_automorphisms),
                    ';',
                )
                println(
                    io,
                    join((
                        result.graph6,
                        class_index,
                        candidate_class.representative_index,
                        members,
                        automorphisms,
                        "interval_distance_compatibility",
                    ), '\t'),
                )
            end
        end
    end
    return path
end

end

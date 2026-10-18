module DiskSecondMomentProof

using ArbNumerics
using LinearAlgebra

include("core/configuration.jl")
include("core/contact_graph.jl")
include("core/exact_numbers.jl")
include("objective/base.jl")
include("objective/second_moment.jl")
include("core/contact_constraints.jl")
include("numerics/intrinsic_hessian.jl")
include("numerics/interval_safety.jl")
include("numerics/full_equilibrium.jl")
include("io/repository_paths.jl")
include("io/json_parser.jl")
include("evidence/box_records.jl")
include("evidence/archive.jl")
include("census/Census.jl")
include("analysis/candidate_dataset.jl")
include("analysis/regular_candidates.jl")
include("analysis/exact_candidates.jl")
include("analysis/exact_families.jl")
include("analysis/four_disk_rhombus.jl")
include("search/CandidateSearch.jl")
include("search/ExactFamilyCoverage.jl")
const FamilyCoverage = ExactFamilyCoverage
include("refinement/Refinement.jl")
include("pipeline/candidate_pipeline.jl")
include("io/graph_census.jl")
include("search/objective_free_leaf_exclusion.jl")
include("search/minimum_candidate_search.jl")
include("pipeline/candidate_deduplication.jl")
include("evidence/search_conventions.jl")
include("pipeline/proof_workflow.jl")

export Configuration,
       num_disks,
       as_vector,
       centroid,
       center_at_centroid,
       ensure_centered,
       squared_distance,
       ContactEdge,
       ContactGraph,
       decode_graph6,
       num_contacts,
       degree,
       contact_matrix,
       euclidean_motions,
       AbstractObjective,
       evaluate,
       gradient,
       ambient_hessian,
       SecondMomentObjective,
       second_moment_target,
       constraint_hessians,
       lagrangian_hessian,
       contact_equations,
       contact_jacobian,
       Qsqrt3,
       IntrinsicHessianVerification,
       RegularCandidateResult,
       TargetResult,
       FamilyResult,
       FourDiskEndpointResult,
       verify_four_disk_endpoint,
       verify_four_disk_rhombus_coverage,
       validate_candidate_dataset,
       verify_refined_candidate,
       verify_target,
       verify_family,
       CandidateProofResult,
       ComputationalProofResult,
       analyze_search_results,
       print_proof_summary,
       write_candidate_reports,
       load_connected_graphs,
       connected_graph_path,
       proof_input_path,
       Census,
       FamilyCoverage,
       Refinement,
       CandidateSearch,
       ExactMinimizerResult,
       ExactFamilyResult,
       verify_exact_minimizer,
       verify_exact_family,
       exact_input_path,
       ExactFamilyCoverage

end

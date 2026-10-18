"""Return the package root directory."""
package_root() = normpath(joinpath(@__DIR__, "..", ".."))

"""Return the path of the connected graph6 census for `n`."""
function connected_graph_path(n::Int)
    n in (2, 3, 4, 5, 6) || throw(ArgumentError("n must be 2, 3, 4, 5 or 6"))
    return joinpath(package_root(), "data", "connected_graphs", "graphs$n.txt")
end

"""Return the path of the algebraic proof inputs for `n`."""
function proof_input_path(n::Int)
    n in (2, 3, 4, 5, 6) || throw(ArgumentError("n must be 2, 3, 4, 5 or 6"))
    return joinpath(package_root(), "data", "exact_inputs", "$(n)disks.json")
end

# Public API alias.
const exact_input_path = proof_input_path

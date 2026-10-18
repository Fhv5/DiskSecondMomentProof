"""
    constraint_hessians(c, g; contact_distance=2)

Builds the constant Hessian of each contact equation g_ij = (||c_i-c_j||^2-4)/2

    ||c_i-c_j||^2/2-contact_distance^2/2.

The returned matrices use the same normalization as `contact_equations` and
`contact_matrix`.
"""
function constraint_hessians(c::Configuration{T}, g::ContactGraph;
                             contact_distance=2) where {T}
    n = num_disks(c)
    scale = one(T)
    result = Matrix{T}[]

    for edge in g.edges
        H = zeros(T, 2 * n, 2 * n)
        i, j = edge.src, edge.dst
        for coordinate in 0:1
            ii = 2 * i - 1 + coordinate
            jj = 2 * j - 1 + coordinate
            H[ii, ii] = scale
            H[jj, jj] = scale
            H[ii, jj] = -scale
            H[jj, ii] = -scale
        end
        push!(result, H)
    end
    return result
end

"""
    lagrangian_hessian(H_objective, lambda, c, g; contact_distance=2)

Builds the ambient Hessian of the Lagrangian using exactly the same constraint
normalization as the KKT system.
"""
function lagrangian_hessian(H_objective::AbstractMatrix{T},
                            lambda::AbstractVector{T},
                            c::Configuration{T}, g::ContactGraph;
                            contact_distance=2) where {T}
    length(lambda) == num_contacts(g) || throw(DimensionMismatch("one multiplier is required per contact"))
    H = Matrix{T}(H_objective)
    for (multiplier, constraint_hessian) in
            zip(lambda, constraint_hessians(c, g; contact_distance=contact_distance))
        H .-= multiplier .* constraint_hessian
    end
    return H
end

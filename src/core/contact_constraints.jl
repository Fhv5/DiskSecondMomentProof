"""
    contact_matrix(c::Configuration{T}, g::ContactGraph; contact_distance=2) where {T}

Construct the Jacobian of the contact constraints g_ij = (||c_i-c_j||^2-4)/2
`g_ij(c) = norm(c_i - c_j)^2 / 2 - contact_distance^2 / 2`.
Each row contains `c_i - c_j` at vertex `i` and its negative at vertex `j`.
"""
function contact_matrix(c::Configuration{T}, g::ContactGraph; contact_distance=2) where {T}
    n = num_disks(c)
    m = num_contacts(g)
    A = zeros(T, m, 2 * n)

    for (k, edge) in enumerate(g.edges)
        i, j = edge.src, edge.dst
        u_x = c.centers[1, i] - c.centers[1, j]
        u_y = c.centers[2, i] - c.centers[2, j]

        A[k, 2 * i - 1] = u_x
        A[k, 2 * i]     = u_y
        A[k, 2 * j - 1] = -u_x
        A[k, 2 * j]     = -u_y
    end

    return A
end

"""
    contact_equations(c::Configuration{T}, g::ContactGraph; contact_distance=2) where {T}

Return the contact equations `norm(c_i - c_j)^2 / 2 - contact_distance^2 / 2`.
"""
function contact_equations(c::Configuration{T}, g::ContactGraph; contact_distance=2) where {T}
    values = zeros(T, num_contacts(g))
    distance_squared = T(contact_distance)^2
    for (k, edge) in enumerate(g.edges)
        dx = c.centers[1, edge.src] - c.centers[1, edge.dst]
        dy = c.centers[2, edge.src] - c.centers[2, edge.dst]
        values[k] = (dx^2 + dy^2 - distance_squared) / T(2)
    end
    return values
end

"""
    contact_jacobian(c::Configuration, g::ContactGraph; contact_distance=2)

Return the Jacobian of `contact_equations`, equal to `contact_matrix`.
"""
contact_jacobian(c::Configuration, g::ContactGraph; contact_distance=2) =
    contact_matrix(c, g; contact_distance=contact_distance)

"""
    euclidean_motions(c::Configuration{T}) where {T}

Return the (2n x 3) matrix of planar rigid motions: two translations and
infinitesimal rotation around the origin.
"""
function euclidean_motions(c::Configuration{T}) where {T}
    n = num_disks(c)
    E = zeros(T, 2 * n, 3)

    for i in 1:n
        E[2 * i - 1, 1] = one(T)
        E[2 * i, 2]     = one(T)
        E[2 * i - 1, 3] = -c.centers[2, i]
        E[2 * i, 3]     = c.centers[1, i]
    end

    return E
end

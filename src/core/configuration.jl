"""
    Configuration{T}

Stores the 2D coordinates of n disk centers as a (2  x  n) matrix of type T.
Column j corresponds to [x_j, y_j]^T.
"""
struct Configuration{T}
    centers::Matrix{T}

    function Configuration(centers::Matrix{T}) where {T}
        @assert size(centers, 1) == 2 "Configuration matrix must have 2 rows (x and y coordinates)."
        new{T}(centers)
    end
end

"""
    Configuration(vec::AbstractVector{T})

Constructs a Configuration by reshaping a flat 2n vector [x1, y1, x2, y2, ...] into a (2  x  n) matrix.
"""
function Configuration(vec::AbstractVector{T}) where {T}
    @assert iseven(length(vec)) "Vector length must be even (2n)."
    n = div(length(vec), 2)
    centers = reshape(copy(vec), 2, n)
    return Configuration(centers)
end

"""
    num_disks(c::Configuration)

Returns the number of disks n in the configuration.
"""
function num_disks(c::Configuration)
    return size(c.centers, 2)
end

"""
    as_vector(c::Configuration)

Returns the configuration centers flattened into a 1D vector of length 2n.
"""
function as_vector(c::Configuration)
    return vec(c.centers)
end

"""
    centroid(c::Configuration{T}) where {T}

Calculates the centroid [x_bar, y_bar] of the disk cluster as a 1D vector of length 2.
"""
function centroid(c::Configuration{T}) where {T}
    n = num_disks(c)
    return vec(sum(c.centers, dims=2) ./ n)
end

"""
    center_at_centroid(c::Configuration{T}) where {T}

Translates the entire configuration so that its centroid is at the origin (0, 0).
"""
function center_at_centroid(c::Configuration{T}) where {T}
    c_bar = centroid(c)
    return Configuration(c.centers .- c_bar)
end

"""
    ensure_centered(c::Configuration{T}) where {T}

Verifies if the configuration is centered at the origin.
If the centroid is strictly zero, returns the configuration as-is.
Otherwise, translates the configuration so that its centroid is at the origin.
"""
function ensure_centered(c::Configuration{T}) where {T}
    c_bar = centroid(c)
    if iszero(c_bar[1]) && iszero(c_bar[2])
        return c
    end
    return center_at_centroid(c)
end

"""
    squared_distance(c::Configuration, i::Int, j::Int)

Computes the squared Euclidean distance (Deltax^2 + Deltay^2) between disk i and disk j.
"""
function squared_distance(c::Configuration, i::Int, j::Int)
    dx = c.centers[1, i] - c.centers[1, j]
    dy = c.centers[2, i] - c.centers[2, j]
    return dx^2 + dy^2
end

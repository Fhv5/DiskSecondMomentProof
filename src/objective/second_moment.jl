"""
    SecondMomentObjective <: AbstractObjective

Represents the centroidal second moment objective function:
    M(c) = (1/2) * sum_{i=1}^n ||c_i - c_bar||^2 = (1 / 2n) * sum_{1<=i<j<=n} ||c_i - c_j||^2
"""
struct SecondMomentObjective <: AbstractObjective end

"""Return the announced second moment target as a rational number."""
function second_moment_target(n::Int)
    n == 2 && return BigInt(1) // BigInt(1)
    n == 3 && return BigInt(2) // BigInt(1)
    n == 4 && return BigInt(4) // BigInt(1)
    n == 5 && return BigInt(34) // BigInt(5)
    n == 6 && return BigInt(29) // BigInt(3)
    throw(ArgumentError("no announced second-moment target is configured for n=$n"))
end

"""
    evaluate(::SecondMomentObjective, c::Configuration{T}) where {T}

Calculates the scalar second moment M(c) = (1/2) * sum ||c_i - c_bar||^2.
"""
function evaluate(::SecondMomentObjective, c::Configuration{T}) where {T}
    if T <: AbstractFloat
        @warn "[DiskSecondMomentProof: objective] `evaluate` is operating with floating-point arithmetic."
    end

    n = num_disks(c)
    c_bar = centroid(c)
    total = zero(T)

    for i in 1:n
        dx = c.centers[1, i] - c_bar[1]
        dy = c.centers[2, i] - c_bar[2]
        total += dx^2 + dy^2
    end

    return total / 2
end

"""
    gradient(::SecondMomentObjective, c::Configuration{T}) where {T}

Computes the 2n-dimensional gradient vector gradientM(c):
    gradient_{c_i} M(c) = c_i - c_bar
"""
function gradient(::SecondMomentObjective, c::Configuration{T}) where {T}
    if T <: AbstractFloat
        @warn "[DiskSecondMomentProof: objective] `gradient` is operating with floating-point arithmetic."
    end

    n = num_disks(c)
    c_bar = centroid(c)
    grad = zeros(T, 2 * n)

    for i in 1:n
        grad[2 * i - 1] = c.centers[1, i] - c_bar[1]
        grad[2 * i]     = c.centers[2, i] - c_bar[2]
    end

    return grad
end

"""
    ambient_hessian(::SecondMomentObjective, c::Configuration{T}) where {T}

Computes the constant (2n  x  2n) ambient Hessian matrix gradient^2M(c):
    d^2M / dc_i^2 = (1 - 1/n) * I_2
    d^2M / dc_i dc_j = - (1/n) * I_2   (for i != j)
"""
function ambient_hessian(::SecondMomentObjective, c::Configuration{T}) where {T}
    if T <: AbstractFloat
        @warn "[DiskSecondMomentProof: objective] `ambient_hessian` is operating with floating-point arithmetic."
    end

    n = num_disks(c)
    H = zeros(T, 2 * n, 2 * n)
    diag_val = (T(n) - one(T)) / T(n)
    offdiag_val = -one(T) / T(n)

    for i in 1:n
        for j in 1:n
            val = (i == j) ? diag_val : offdiag_val
            H[2 * i - 1, 2 * j - 1] = val
            H[2 * i,     2 * j]     = val
        end
    end

    return H
end

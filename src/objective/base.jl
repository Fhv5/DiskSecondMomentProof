"""
    AbstractObjective

Base type for objective functions supported by the analysis. A concrete
objective implements `evaluate`, `gradient`, and `ambient_hessian`.
"""
abstract type AbstractObjective end

"""Evaluate the objective at a configuration."""
function evaluate(obj::AbstractObjective, c::Configuration)
    error("`evaluate` is not implemented for objective type $(typeof(obj)).")
end

"""Return the objective gradient in Cartesian coordinates."""
function gradient(obj::AbstractObjective, c::Configuration)
    error("`gradient` is not implemented for objective type $(typeof(obj)).")
end

"""Return the unconstrained Cartesian Hessian of the objective."""
function ambient_hessian(obj::AbstractObjective, c::Configuration)
    error("`ambient_hessian` is not implemented for objective type $(typeof(obj)).")
end

"""
    IceColumn

Solution to the ice-column heat equation (Moreno-Parada et al., 2024).

Always carries both dimensional (K) and dimensionless (θ = T/T_air) temperatures,
as well as the equilibrium (stationary) solution.

# Fields
- `z`       : depth coordinate (m), length nz; z=0 at base, z=L at surface
- `zeta`    : normalised coordinate ζ = z/L ∈ [0,1], length nz
- `t`       : time vector (yr), length nt
- `T`       : temperature (K), size nz × nt
- `theta`   : dimensionless temperature θ = T/T_air, size nz × nt
- `T_eq`    : equilibrium (stationary) temperature (K), length nz
- `theta_eq`: dimensionless equilibrium temperature, length nz
- `par`     : `IceColumnPar` used to compute this solution
"""
struct IceColumn
    z::Vector{Float64}
    zeta::Vector{Float64}
    t::Vector{Float64}
    T::Matrix{Float64}
    theta::Matrix{Float64}
    T_eq::Vector{Float64}
    theta_eq::Vector{Float64}
    par::IceColumnPar
end

function Base.show(io::IO, sol::IceColumn)
    nz, nt = size(sol.T)
    print(io, "IceColumn(nz=$nz, nt=$nt, t=$(extrema(sol.t)) yr | $(sol.par))")
end

"""
    to_celsius(sol::IceColumn) -> IceColumn

Return a copy of `sol` with `T` and `T_eq` converted to degrees Celsius.
`theta` and `theta_eq` are unchanged (always dimensionless).
"""
function to_celsius(sol::IceColumn)
    IceColumn(sol.z, sol.zeta, sol.t,
              sol.T .- 273.15,
              sol.theta,
              sol.T_eq .- 273.15,
              sol.theta_eq,
              sol.par)
end

"""
    to_celsius!(sol::IceColumn)

Convert `sol.T` and `sol.T_eq` to degrees Celsius in-place.
Returns `sol`.
"""
function to_celsius!(sol::IceColumn)
    sol.T    .-= 273.15
    sol.T_eq .-= 273.15
    sol
end

# ---- internal helpers -------------------------------------------------------

"""Build a uniform ζ grid on [0,1] of length nz."""
_zeta_grid(nz::Int) = range(0.0, 1.0, length=nz) |> collect

"""
Resolve the vertical grid from keyword arguments `nz`, `zeta`, and `z`.

`nz == 0` is the sentinel for "not explicitly set"; when no custom grid is
given it defaults to a uniform grid of 100 points. Passing both `zeta` and
`z` is an error. Passing `nz > 0` alongside a custom grid whose length
differs from `nz` is an error.
"""
function _resolve_zeta_grid(nz::Int,
                             zeta::Union{Nothing, AbstractVector{<:Real}},
                             z::Union{Nothing, AbstractVector{<:Real}},
                             L::Real)
    !isnothing(zeta) && !isnothing(z) &&
        error("provide `zeta` or `z`, not both")
    if !isnothing(zeta)
        nz > 0 && nz != length(zeta) &&
            error("`nz=$nz` conflicts with length(zeta)=$(length(zeta))")
        return collect(Float64, zeta)
    end
    if !isnothing(z)
        nz > 0 && nz != length(z) &&
            error("`nz=$nz` conflicts with length(z)=$(length(z))")
        return collect(Float64, z) ./ L
    end
    _zeta_grid(nz == 0 ? 100 : nz)
end

"""Assemble an IceColumn from dimensionless θ arrays and par."""
function _make_solution(zeta::Vector{Float64},
                        theta_eq::Vector{Float64},
                        theta::Matrix{Float64},
                        t::Vector{Float64},
                        par::IceColumnPar;
                        celsius::Bool=false)
    z    = zeta .* par.L
    T_eq = theta_eq .* par.T_air
    T    = theta    .* par.T_air
    if celsius
        T_eq .-= 273.15
        T    .-= 273.15
    end
    IceColumn(z, zeta, t, T, theta, T_eq, theta_eq, par)
end

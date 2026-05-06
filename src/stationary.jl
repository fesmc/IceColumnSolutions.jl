"""
Stationary (equilibrium) solution — Appendix B, Moreno-Parada et al. (2024).

The stationary temperature profile ϑ(ξ) satisfies (Eq. B1 in non-dimensional form):
    ϑ_ξξ - Pe·ξ·ϑ_ξ = Ω,   ξ ∈ [0,1]
    ϑ_ξ = γ,               ξ = 0
    β'·ϑ_ξ + ϑ = 1,        ξ = 1

where Ω = Br + Λ is the total dimensionless heat source.

Solution (Eq. B2):
    ϑ(ξ) = Ω·(ξ²/2)·₂F₂(1,1; 3/2, 2; -a²ξ²) + A·erf(aξ) + B

with a = sqrt(Pe/2), A = -γ·sqrt(π/(4a)), B from the top BC.

Special case Pe=0: purely diffusive + source gives a quadratic polynomial.

All computations use complex arithmetic where Pe < 0 (downward ice flow) makes
a imaginary; the physical result is real and recovered via `real()`.
"""

using HypergeometricFunctions: pFq
using SpecialFunctions: erf

# ---- helpers ----------------------------------------------------------------

"""₂F₂(1,1; 3/2, 2; z) — generalised hypergeometric function."""
@inline _2F2_1122(z) = pFq((1.0, 1.0), (1.5, 2.0), z)

"""₂F₂(2,2; 5/2, 3; z)"""
@inline _2F2_2253(z) = pFq((2.0, 2.0), (2.5, 3.0), z)

"""
Stationary solution at a vector of ζ values for given dimensionless parameters.
Returns a real vector of θ values.
"""
function _stationary_theta(zeta::AbstractVector{Float64}, par::IceColumnPar)
    Pe = par.Pe
    Ω  = par.Br + par.Lambda
    γ  = par.gamma
    β  = par.beta_prime

    if abs(Pe) < 1e-12
        return _stationary_theta_Pe0(zeta, Ω, γ, β)
    end

    a  = sqrt(complex(Pe / 2.0))   # imaginary for Pe < 0 (downward flow)
    a2 = a^2                        # = Pe/2

    # Constants A and B from Appendix B (Eq. B2)
    # d/dξ[erf(aξ)]|_{ξ=0} = 2a/√π, so BC ϑ_ξ(0)=γ gives A = γ√π/(2a)
    A = γ * sqrt(complex(π)) / (2 * a)

    # ϑ(1) contribution from ₂F₂ term:  Ω/2 · ₂F₂(1,1;3/2,2; -a²)
    F1_at1 = _2F2_1122(-a2)

    # ϑ_ξ(1) = Ω·[₂F₂(1,1;3/2,2;-a²) - a²/3·₂F₂(2,2;5/2,3;-a²)] + A·(2a/√π)·e^(-a²)
    F2_at1 = _2F2_2253(-a2)
    dF_at1 = F1_at1 - (a2 / 3.0) * F2_at1

    erf_a  = erf(a)
    exp_a2 = exp(-a2)

    B = (1.0
         - A * (2a * exp_a2 / sqrt(π) * β + erf_a)
         - Ω * ((β + 0.5) * F1_at1 - β * a2 / 3.0 * F2_at1))

    # Evaluate ϑ at each ζ
    θ = similar(zeta, ComplexF64)
    for i in eachindex(zeta)
        ξ  = zeta[i]
        ζi = a2 * ξ^2      # = Pe/2 · ξ²
        θ[i] = Ω * (ξ^2 / 2) * _2F2_1122(-ζi) + A * erf(a * ξ) + B
    end

    real.(θ)
end

"""Stationary solution for Pe = 0 (pure diffusion + source)."""
function _stationary_theta_Pe0(zeta, Ω, γ, β)
    # ϑ(ξ) = -Ω·ξ²/2 + γ·ξ + C
    # BC at ξ=1: β·(-Ω + γ) + (-Ω/2 + γ + C) = 1  →  C = 1 + (β+0.5)·Ω - (β+1)·γ
    C = 1.0 + (β + 0.5) * Ω - (β + 1.0) * γ
    [-Ω * ξ^2 / 2 + γ * ξ + C for ξ in zeta]
end

# ---- public API -------------------------------------------------------------

"""
    solve_stationary(par::IceColumnPar; nz=100, celsius=false) -> IceColumn

Compute the equilibrium (stationary) temperature profile ϑ(ξ).

Returns an `IceColumn` where `T` == `T_eq` (nz × 1 matrix, repeated once so the
struct is always consistent). Time vector `t` contains a single value `Inf`.
"""
function solve_stationary(par::IceColumnPar; nz::Int=100, celsius::Bool=false)
    zeta     = _zeta_grid(nz)
    theta_eq = _stationary_theta(zeta, par)
    theta    = reshape(theta_eq, :, 1)   # nz × 1
    _make_solution(zeta, theta_eq, theta, [Inf], par; celsius=celsius)
end

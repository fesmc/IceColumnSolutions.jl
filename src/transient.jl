"""
Transient solution — Appendix A, Moreno-Parada et al. (2024).

The transient perturbation μ(ξ,τ) = θ(ξ,τ) - ϑ(ξ) satisfies:
    μ_τ = μ_ξξ - Pe·ξ·μ_ξ,   ξ ∈ [0,1], τ > 0
    μ_ξ = 0,                   ξ = 0
    β'·μ_ξ + μ = 0,            ξ = 1
    μ(ξ,0) = θ₀(ξ) - ϑ(ξ)

Pe ≠ 0 solution (Eq. A7):
    μ(ξ,τ) = Σ_n  Aₙ · M(αₙ; 1/2; Pe·ξ²/2) · exp(λₙ·τ)

where M = Kummer ₁F₁, λₙ = 2·Pe·αₙ < 0 (decay).

Eigenvalues αₙ from BC at ξ=1 (Eq. A8, corrected sign):
    β'·Pe·2αₙ·M(αₙ+1, 3/2, Pe/2) + M(αₙ, 1/2, Pe/2) = 0

The Kummer ODE gives X'' - Pe·ξ·X' = 2Pe·α·X, so λₙ = 2Pe·αₙ.
For decay (λₙ < 0):  αₙ·Pe < 0  →  αₙ < 0 for Pe > 0, αₙ > 0 for Pe < 0.

Pe = 0 solution (pure diffusion):
    μ(ξ,τ) = Σ_n  Aₙ · cos(kₙ·ξ) · exp(-kₙ²·τ)

where kₙ satisfies: -β'·kₙ·sin(kₙ) + cos(kₙ) = 0.
For β'=0 (Dirichlet top): kₙ = (n-1/2)π exactly.
"""

using HypergeometricFunctions: M as Kummer
using QuadGK: quadgk
using Roots: find_zeros

# ---- initial condition helpers -----------------------------------------------

"""
    uniform(T_val) -> Function

Initial condition: uniform temperature `T_val` (Kelvin).
Returns a function `θ₀(ξ, par)` = T_val / T_air.
"""
uniform(T_val::Real) = (ξ, par) -> T_val / par.T_air

"""
    stationary_init(par) -> Function

Initial condition: the stationary (equilibrium) temperature profile ϑ(ξ).
Produces no transient modes (μ = 0 everywhere).
"""
stationary_init(par::IceColumnPar) = (ξ, par2) -> _stationary_theta_scalar(ξ, par2)

# ---- eigenvalue machinery ---------------------------------------------------

"""
Eigenvalue equation residual at ξ=1 for a given α (Eq. A8, Pe ≠ 0).
Returns zero when α is a valid eigenvalue parameter.
"""
function _eigen_residual(α, par::IceColumnPar)
    Pe = par.Pe
    β  = par.beta_prime
    z  = complex(Pe / 2)
    real(β * Pe * 2α * Kummer(α + 1, 1.5, z) + Kummer(α, 0.5, z))
end

"""Pe=0 eigenvalue equation residual: -β'·k·sin(k) + cos(k) = 0."""
_eigen_residual_Pe0(k, β) = -β * k * sin(k) + cos(k)

"""
Find the first `n` wavenumbers kₙ for Pe=0.
For β'=0: kₙ = (n-1/2)π. For β'>0: root-find per interval.
"""
function _find_k_Pe0(β, n::Int)
    if β < 1e-12
        return [(i - 0.5) * π for i in 1:n]
    end
    ks = Float64[]
    for i in 1:n
        lo = (i - 1) * π + 1e-6
        hi = i * π - 1e-6
        append!(ks, find_zeros(k -> _eigen_residual_Pe0(k, β), lo, hi))
    end
    sort!(ks)[1:n]
end

"""
Find the first `n` eigenvalue parameters αₙ for Pe ≠ 0.

Decay requires λₙ = 2Pe·αₙ < 0:
  Pe > 0  →  αₙ < 0  (search negative range)
  Pe < 0  →  αₙ > 0  (search positive range)

The scan range is set from physical estimates: eigenvalues σₙ ≈ -(nπ)²
at large n, so |αₙ| ≈ n²π²/(2|Pe|). A scan grid of ~50·n points
suffices to detect all sign changes with adequate spacing.
"""
function _find_alpha_values(par::IceColumnPar, n::Int)
    Pe      = par.Pe
    scan_max = max(50.0, 1.5 * n^2 * π^2 / (2 * abs(Pe)))
    n_grid   = max(1000, 50 * n)

    # αₙ has opposite sign to Pe for decay
    if Pe > 0
        grid = range(-scan_max, -1e-6, length=n_grid)
    else
        grid = range(1e-6, scan_max, length=n_grid)
    end

    alphas = find_zeros(α -> _eigen_residual(α, par), grid)
    if length(alphas) < n
        error("Only found $(length(alphas)) eigenvalues; requested $n. " *
              "Try increasing scan range.")
    end
    sort!(alphas, by=abs)[1:n]
end

"""
    eigenvalues(par, n) -> (alphas_or_ks, lambdas)

Return the first `n` eigenvalue parameters and decay rates λₙ < 0.

For Pe ≠ 0: `alphas_or_ks[n] = αₙ` (Kummer parameter), `lambdas[n] = 2·Pe·αₙ < 0`.
For Pe = 0: `alphas_or_ks[n] = kₙ` (wavenumber), `lambdas[n] = -kₙ²`.
"""
function eigenvalues(par::IceColumnPar, n::Int)
    if abs(par.Pe) < 1e-12
        ks      = _find_k_Pe0(par.beta_prime, n)
        lambdas = [-k^2 for k in ks]
        return ks, lambdas
    end
    alphas  = _find_alpha_values(par, n)
    lambdas = [2.0 * par.Pe * α for α in alphas]   # λₙ = 2Pe·αₙ < 0
    alphas, lambdas
end

# ---- eigenfunctions & coefficients ------------------------------------------

"""
Eigenfunction at ξ given the stored eigenvalue parameter α_or_k.
  Pe ≠ 0: M(α_or_k, 1/2, Pe·ξ²/2)
  Pe = 0: cos(α_or_k · ξ)
"""
function _eigenfunction(ξ::Real, α_or_k::Real, par::IceColumnPar)
    if abs(par.Pe) < 1e-12
        return cos(α_or_k * ξ)
    end
    real(Kummer(complex(α_or_k), 0.5, complex(par.Pe * ξ^2 / 2)))
end

"""Weight function ϱ(ξ) = exp(-Pe·ξ²/2)."""
_weight(ξ::Real, Pe::Real) = exp(-Pe * ξ^2 / 2)

"""
Compute series coefficient Aₙ (Eq. A9) by weighted inner product.
"""
function _series_coefficient(α_or_k, par, delta0_fn)
    Φ(ξ)  = _eigenfunction(ξ, α_or_k, par)
    ϱ(ξ)  = _weight(ξ, par.Pe)
    num, _ = quadgk(ξ -> delta0_fn(ξ) * ϱ(ξ) * Φ(ξ), 0.0, 1.0; rtol=1e-8)
    den, _ = quadgk(ξ -> Φ(ξ)^2 * ϱ(ξ), 0.0, 1.0; rtol=1e-8)
    num / den
end

# ---- scalar stationary helper -----------------------------------------------

"""Stationary θ at a single ξ value."""
function _stationary_theta_scalar(ξ::Real, par::IceColumnPar)
    _stationary_theta([Float64(ξ)], par)[1]
end

# ---- main solver ------------------------------------------------------------

"""
    solve(par, ts; init, n_modes=50, nz=100, celsius=false) -> IceColumn

Compute the transient temperature evolution θ(ξ,τ) at times `ts` (years).

# Arguments
- `par`     : `IceColumnPar` with all parameters
- `ts`      : time vector (yr) at which to evaluate the solution
- `init`    : initial condition function `θ₀(ξ, par)` (dimensionless);
              use `uniform(T_val)` or `stationary_init(par)`.
- `n_modes` : number of series terms (default 50)
- `nz`      : number of vertical grid points (default 100)
- `celsius` : if true, return T and T_eq in °C (default false)
"""
function solve(par::IceColumnPar, ts::AbstractVector{<:Real};
               init,
               n_modes::Int = 50,
               nz::Int      = 100,
               celsius::Bool = false)

    zeta     = _zeta_grid(nz)
    theta_eq = _stationary_theta(zeta, par)

    # Non-dimensional times τ = κ·t / L²
    kappa_yr = par.kappa * YR_TO_S
    tau_vec  = [kappa_yr * t / par.L^2 for t in ts]

    # Initial perturbation δ₀(ξ) = θ₀(ξ) - ϑ(ξ)
    delta0(ξ) = init(ξ, par) - _stationary_theta_scalar(ξ, par)

    # Eigenvalues
    params_eig, lambdas = eigenvalues(par, n_modes)

    # Series coefficients
    An = [_series_coefficient(params_eig[n], par, delta0) for n in 1:n_modes]

    # Precompute eigenfunction matrix Φ[i,n] = Φₙ(ξᵢ)
    Phi = Matrix{Float64}(undef, nz, n_modes)
    for n in 1:n_modes, i in eachindex(zeta)
        Phi[i, n] = _eigenfunction(zeta[i], params_eig[n], par)
    end

    # Evaluate θ(ξ, τ) = ϑ(ξ) + Phi * (An .* exp.(λ·τ))
    nt    = length(ts)
    theta = Matrix{Float64}(undef, nz, nt)
    An_v  = collect(Float64, An)

    for j in eachindex(tau_vec)
        decay        = exp.(lambdas .* tau_vec[j])   # n_modes vector
        theta[:, j]  = theta_eq .+ Phi * (An_v .* decay)
    end

    t_vec = collect(Float64, ts)
    _make_solution(zeta, theta_eq, theta, t_vec, par; celsius=celsius)
end

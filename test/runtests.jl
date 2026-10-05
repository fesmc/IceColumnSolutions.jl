using IceColumnSolutions
using Test

# Shared dimensional anchors
const L     = 1000.0
const T_air = 250.0
const kappa = 1e-6
const k     = 2.0

# ---- stationary solution ----------------------------------------------------

@testset "Stationary — Pe=0 pure diffusion" begin
    # With Pe=0, Br=0, Λ=0 the stationary solution is the linear profile
    # ϑ(ξ) = γ·ξ + C. With β'=0 (Dirichlet): C = 1 - γ, so ϑ(0)=1-γ, ϑ(1)=1.
    γ = 2.0
    par = IceColumnPar(L, T_air, kappa, k, 0.0, γ, 0.0)
    sol = solve_stationary(par; nz=101)

    # ϑ(1) = 1 (Dirichlet top BC)
    @test sol.theta_eq[end] ≈ 1.0  atol=1e-12
    # ϑ(0) = 1 - γ (base value from linear profile)
    @test sol.theta_eq[1] ≈ 1.0 - γ  atol=1e-12
    # Profile is linear: slope = γ
    dtheta = diff(sol.theta_eq) ./ diff(sol.zeta)
    @test all(isapprox.(dtheta, γ; atol=1e-10))
end

@testset "Stationary — Neumann base BC" begin
    # θ_ξ(0) = γ should hold for all parameter sets
    γ = 1.5
    for Pe in [-3.0, 0.0, 3.0]
        par = IceColumnPar(L, T_air, kappa, k, 0.2, γ, Pe; Br=0.5)
        sol = solve_stationary(par; nz=201)
        dxi = sol.zeta[2] - sol.zeta[1]
        slope_base = (sol.theta_eq[2] - sol.theta_eq[1]) / dxi
        @test slope_base ≈ γ  atol=5e-3   # O(h) forward-difference, h=1/200
    end
end

@testset "Stationary — Robin top BC" begin
    # β'·θ_ξ(1) + θ(1) = 1
    for (Pe, β) in [(-5.0, 0.5), (5.0, 1.0), (0.0, 0.3)]
        par = IceColumnPar(L, T_air, kappa, k, β, 2.0, Pe)
        sol = solve_stationary(par; nz=501)
        n   = length(sol.zeta)
        h   = sol.zeta[n] - sol.zeta[n-1]
        # 3-point O(h²) backward difference to reduce truncation error
        slope_top = (3sol.theta_eq[n] - 4sol.theta_eq[n-1] + sol.theta_eq[n-2]) / (2h)
        bc_val = β * slope_top + sol.theta_eq[n]
        @test bc_val ≈ 1.0  atol=1e-2
    end
end

@testset "Stationary — benchmark exp1 (Pe=0,Br=0)" begin
    par = benchmark(:exp1)
    sol = solve_stationary(par)
    @test sol.theta_eq[end] ≈ 1.0  atol=1e-12
    @test sol.theta_eq[1]   ≈ 1.0 - par.gamma  atol=1e-12
end

# ---- transient solution ------------------------------------------------------

@testset "Transient — convergence to stationary" begin
    # Starting from uniform T, solution should approach stationary at large τ
    par = benchmark(:exp1)
    T0  = 0.9 * T_air   # slightly off from equilibrium
    # t_final chosen so λ₁·τ ≈ -(π/2)²·κ_yr·t/L² << -10: need t ≈ 200_000 yr
    ts  = [0.0, 100.0, 1000.0, 300_000.0]
    sol = solve(par, ts; init=uniform(T0), n_modes=50, nz=51)

    # At t=0 the series approximates θ₀ in the interior.  The uniform IC violates
    # the Dirichlet BC (θ₀(1)=0.9 ≠ 1=ϑ(1)), so eigenfunctions (all zero at ξ=1)
    # can never reproduce the boundary value — skip the last grid point.
    interior = 1:length(sol.zeta)-1
    @test all(isapprox.(sol.theta[interior, 1], T0 / T_air; atol=5e-2))

    # At large t the solution should converge to the equilibrium profile
    theta_last = sol.theta[:, end]
    @test maximum(abs.(theta_last .- sol.theta_eq)) < 1e-4
end

@testset "Transient — stationary init gives no change" begin
    # Starting exactly at equilibrium: all Aₙ ≈ 0, solution stays constant
    par = benchmark(:exp2)
    ts  = [0.0, 1000.0, 10_000.0]
    sol = solve(par, ts; init=stationary_init(par), n_modes=5, nz=51)

    for j in 1:length(ts)
        @test maximum(abs.(sol.theta[:, j] .- sol.theta_eq)) < 1e-4
    end
end

@testset "Transient — BC satisfied at all times" begin
    # Robin BC: β'·θ_ξ(1,t) + θ(1,t) ≈ 1
    # exp4 has Pe=5 (downward flow) with steep gradients at ξ=1; use 3-point
    # backward difference to reduce finite-difference truncation error.
    par = benchmark(:exp4)
    ts  = [10.0, 500.0, 5000.0]
    sol = solve(par, ts; init=uniform(0.8 * T_air), n_modes=5, nz=201)
    β   = par.beta_prime
    n   = size(sol.theta, 1)
    h   = sol.zeta[n] - sol.zeta[n-1]

    for j in eachindex(ts)
        slope = (3sol.theta[n,j] - 4sol.theta[n-1,j] + sol.theta[n-2,j]) / (2h)
        bc    = β * slope + sol.theta[n, j]
        @test bc ≈ 1.0  atol=5e-2
    end
end

# ---- comparison with a finite-difference solution ----------------------------
#
# Independent reference: θ_τ = θ_ξξ + Pe·ξ·θ_ξ + Ω (Pe > 0 for downward flow),
# θ_ξ(0) = γ, β'·θ_ξ(1) + θ(1) = 1, second-order finite differences on a uniform
# grid with ghost points for both boundary conditions, backward Euler in time.

"Solve a tridiagonal system (sub-, main and super-diagonal a, b, c) by the Thomas algorithm."
function thomas(a, b, c, d)
    n = length(b); cp = zeros(n); dp = zeros(n); x = zeros(n)
    cp[1] = c[1] / b[1]; dp[1] = d[1] / b[1]
    for i in 2:n
        m = b[i] - a[i] * cp[i-1]
        cp[i] = i < n ? c[i] / m : 0.0
        dp[i] = (d[i] - a[i] * dp[i-1]) / m
    end
    x[n] = dp[n]
    for i in n-1:-1:1
        x[i] = dp[i] - cp[i] * x[i+1]
    end
    return x
end

"""
Finite-difference operator A·θ + r for the interior and boundary nodes
ξ_i = i/N, i = 0…N (Dirichlet top for β' = 0: θ_N = 1 is kept fixed).
Returns the tridiagonal coefficients (a, b, c) and the source r.
"""
function fd_operator(par, N)
    h = 1.0 / N; Pe = par.Pe; Ω = par.Br + par.Lambda; γ = par.gamma; β = par.beta_prime
    a = zeros(N+1); b = zeros(N+1); c = zeros(N+1); r = fill(Ω, N+1)
    for i in 0:N
        ξ = i * h
        lo = 1/h^2 - Pe*ξ/(2h); hi = 1/h^2 + Pe*ξ/(2h); j = i + 1
        a[j] = lo; b[j] = -2/h^2; c[j] = hi
        if i == 0                       # ghost θ_{-1} = θ_1 - 2hγ
            c[j] += lo; r[j] -= lo * 2h * γ; a[j] = 0.0
        elseif i == N && β > 0          # ghost θ_{N+1} = θ_{N-1} + 2h(1 - θ_N)/β
            a[j] += hi; b[j] -= hi * 2h / β; r[j] += hi * 2h / β; c[j] = 0.0
        end
    end
    return a, b, c, r
end

function fd_stationary(par, N)
    a, b, c, r = fd_operator(par, N)
    if par.beta_prime == 0              # Dirichlet: θ_N = 1
        a[end] = 0.0; b[end] = 1.0; r[end] = -1.0
    end
    return thomas(a, b, c, -r)
end

function fd_transient(par, N, θ0, τs; dτ = 1e-5)
    a, b, c, r = fd_operator(par, N)
    θ = copy(θ0); dir = par.beta_prime == 0
    dir && (θ[end] = 1.0)
    out = zeros(N+1, length(τs)); τ = 0.0
    for (n, τn) in enumerate(τs)
        while τ < τn - 1e-12
            d = θ .+ dτ .* r
            aa = -dτ .* a; bb = 1 .- dτ .* b; cc = -dτ .* c
            if dir
                aa[end] = 0.0; bb[end] = 1.0; d[end] = 1.0
            end
            θ = thomas(aa, bb, cc, d); τ += dτ
        end
        out[:, n] = θ
    end
    return out
end

@testset "Stationary — finite-difference reference" begin
    N = 800
    for (Pe, β, Br, Λ) in [(5.0, 0.0, 0.0, 0.0), (20.0, 0.0, 0.0, 0.0), (-5.0, 0.0, 0.0, 0.0),
                           (5.0, 1.0, 0.0, 3.0), (7.0, 0.0, 6.0, 0.0)]
        par = IceColumnPar(L, T_air, kappa, k, β, 2.0, Pe; Br = Br, Lambda = Λ)
        sol = solve_stationary(par; nz = N + 1)
        @test maximum(abs.(sol.theta_eq .- fd_stationary(par, N))) < 1e-4
    end
    # With geothermal heating (γ < 0), downward flow (Pe > 0) brings the base
    # closer to the surface temperature than pure diffusion, upward flow less close
    θb(Pe) = solve_stationary(IceColumnPar(L, T_air, kappa, k, 0.0, -0.5, Pe); nz = 3).theta_eq[1]
    @test θb(5.0) < θb(0.0) < θb(-5.0)
end

@testset "Transient — finite-difference reference" begin
    N  = 400
    τs = [0.02, 0.1, 0.5]
    kappa_yr = kappa * 365.25 * 24 * 3600
    ts = τs .* L^2 ./ kappa_yr
    for (Pe, β) in [(5.0, 0.0), (20.0, 0.0), (-5.0, 0.0), (5.0, 1.0)]
        par = IceColumnPar(L, T_air, kappa, k, β, 2.0, Pe)
        sol = solve(par, ts; init = uniform(0.9 * T_air), n_modes = 15, nz = N + 1)
        ref = fd_transient(par, N, fill(0.9, N + 1), τs)
        @test maximum(abs.(sol.theta .- ref)) < 2e-3
    end
end

@testset "Transient — decay toward the stationary profile" begin
    # The slowest mode must decay on the advective time scale ~1/Pe (in τ),
    # not on a spuriously long one (as with a sign mismatch between the
    # stationary and transient operators).
    for Pe in [5.0, 20.0, 60.0]
        par = IceColumnPar(L, T_air, kappa, k, 0.0, 2.0, Pe)
        _, λ = eigenvalues(par, 3)
        @test all(λ .< 0)
        @test -λ[1] > 0.5 * Pe
    end
end

# ---- eigenvalue equation -----------------------------------------------------

@testset "Eigenvalue equation satisfied (Pe≠0)" begin
    for exp in [:exp2, :exp4]   # Pe=5 for both
        par = benchmark(exp)
        alphas, lambdas = eigenvalues(par, 10)
        # All decay rates must be negative
        @test all(lambdas .< 0)
        # Each α satisfies the eigenvalue equation at ξ=1 (find_zeros accuracy ~1e-5)
        for α in alphas
            @test abs(IceColumnSolutions._eigen_residual(α, par)) < 1e-4
        end
    end
end

@testset "Eigenvalue equation satisfied (Pe=0)" begin
    par = benchmark(:exp1)   # Pe=0, β'=0
    ks, lambdas = eigenvalues(par, 10)
    for (n, k) in enumerate(ks)
        @test k ≈ (n - 0.5) * π  atol=1e-12   # exact formula for Dirichlet
        @test lambdas[n] ≈ -k^2  atol=1e-12
    end
end

# ---- unit conversion ---------------------------------------------------------

@testset "to_celsius" begin
    par = benchmark(:exp1)
    sol = solve_stationary(par)
    sol_c = to_celsius(sol)

    @test sol_c.T_eq ≈ sol.T_eq .- 273.15
    @test sol_c.theta_eq == sol.theta_eq   # dimensionless unchanged
end

@testset "celsius keyword" begin
    par   = benchmark(:exp1)
    sol_K = solve_stationary(par; celsius=false)
    sol_C = solve_stationary(par; celsius=true)

    @test sol_C.T_eq ≈ sol_K.T_eq .- 273.15
end

# ---- benchmarks --------------------------------------------------------------

@testset "benchmark constructors" begin
    # exp1-3 have β'=0 (Dirichlet BC): θ(1) = 1
    for sym in [:exp1, :exp2, :exp3]
        par = benchmark(sym)
        @test par isa IceColumnPar
        sol = solve_stationary(par)
        @test sol isa IceColumn
        @test sol.theta_eq[end] ≈ 1.0  atol=1e-3
    end
    # exp4 has β'=1 (Robin BC): θ(1) ≠ 1 in general; just check it runs
    par4 = benchmark(:exp4)
    @test par4 isa IceColumnPar
    @test solve_stationary(par4) isa IceColumn
    @test_throws ArgumentError benchmark(:exp99)
end

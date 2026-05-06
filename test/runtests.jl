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
    # exp4 has Pe=5 (upward flow) with steep gradients at ξ=1; use 3-point
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

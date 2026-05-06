"""
Benchmark parameter sets from Moreno-Parada et al. (2024).

Table 2 lists four experiments used in Section 4 (convergence study):
  exp1: Pe=0,  Br=0, γ=2, β'=0, Λ=0  (pure diffusion, Dirichlet top)
  exp2: Pe=5,  Br=0, γ=2, β'=0, Λ=0  (advection dominant)
  exp3: Pe=0,  Br=6, γ=2, β'=0, Λ=0  (strain heating)
  exp4: Pe=5,  Br=0, γ=2, β'=1, Λ=3  (advection + insulation + horizontal)

Dimensional anchors (not constrained by the paper; chosen to give round numbers):
  L=1000 m, T_air=250 K, kappa=1e-6 m²/s, k=2.0 W/m/K
"""

"""
    benchmark(exp::Symbol; L=1000.0, T_air=250.0, kappa=1e-6, k=2.0) -> IceColumnPar

Return the `IceColumnPar` for one of the benchmark experiments from Table 2 of
Moreno-Parada et al. (2024). Valid symbols: `:exp1`, `:exp2`, `:exp3`, `:exp4`.
"""
function benchmark(exp::Symbol;
                   L::Float64     = 1000.0,
                   T_air::Float64 = 250.0,
                   kappa::Float64 = 1e-6,
                   k::Float64     = 2.0)

    params = Dict(
        :exp1 => (Pe=0.0,  Br=0.0, gamma=2.0, beta_prime=0.0, Lambda=0.0),
        :exp2 => (Pe=5.0,  Br=0.0, gamma=2.0, beta_prime=0.0, Lambda=0.0),
        :exp3 => (Pe=0.0,  Br=6.0, gamma=2.0, beta_prime=0.0, Lambda=0.0),
        :exp4 => (Pe=5.0,  Br=0.0, gamma=2.0, beta_prime=1.0, Lambda=3.0),
    )

    haskey(params, exp) || throw(ArgumentError(
        "Unknown benchmark '$exp'. Choose from: $(join(keys(params), ", "))"))

    p = params[exp]
    IceColumnPar(L, T_air, kappa, k, p.beta_prime, p.gamma, p.Pe;
                 Br=p.Br, Lambda=p.Lambda)
end

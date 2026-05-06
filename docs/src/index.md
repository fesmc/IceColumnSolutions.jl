# IceColumnSolutions.jl

**Analytical solutions for transient temperature evolution in an ice column.**

This package implements the closed-form solutions derived by
[Moreno-Parada et al. (2024)](https://doi.org/10.5194/tc-18-4215-2024)
for the 1-D advective–diffusive heat equation in a glacier or ice sheet column.
It is intended as a lightweight reference tool for:

- reproducing the paper's benchmark experiments,
- verifying the thermal component of numerical ice-sheet models,
- rapid parameter-space exploration without a full numerical solver.

## Installation

```julia
using Pkg
Pkg.add(url = "https://github.com/fesmc/IceColumnSolutions.jl")
```

## Quick start

```julia
using IceColumnSolutions

# Physical parameters: 1000 m ice column, no flow, geothermal flux only
par = IceColumnPar(1000.0,   # L     — ice thickness (m)
                   250.0,    # T_air — surface temperature (K)
                   1e-6,     # κ     — thermal diffusivity (m² s⁻¹)
                   2.0,      # k     — thermal conductivity (W m⁻¹ K⁻¹)
                   0.0,      # β     — surface insulation (m); 0 = Dirichlet
                   0.06;     # G     — geothermal heat flux (W m⁻²)
                   w0 = -0.3)  # surface velocity (m yr⁻¹); negative = downward

# Equilibrium (stationary) temperature profile
sol_eq = solve_stationary(par; celsius=true)

# Transient evolution from a cold initial state
sol = solve(par, [0.0, 1_000.0, 10_000.0, 50_000.0];
            init    = uniform(220.0),   # uniform 220 K start
            n_modes = 20,
            celsius = true)

# sol.T       — temperature in °C, size nz × nt
# sol.T_eq    — equilibrium temperature in °C, length nz
# sol.z       — depth coordinate (m), 0 at base, L at surface
```

## Contents

| Section | Description |
|---------|-------------|
| [Theory](@ref) | PDE, boundary conditions, dimensionless parameters, solution formulae |
| [Benchmark Experiments](@ref "Benchmark Experiments") | Four test cases from Table 2 of the paper with figures |
| [Model Verification](@ref) | Using the analytical solution to verify numerical ice-sheet models |
| [API Reference](@ref) | Full documentation of all exported functions and types |

## Citation

If you use this package, please cite the original paper:

> Moreno-Parada, D., Robinson, A., Montoya, M., and Alvarez-Solas, J.:
> Analytical solutions for the transient evolution of temperature in an ice column,
> *The Cryosphere*, 18, 4215–4232, 2024.
> <https://doi.org/10.5194/tc-18-4215-2024>

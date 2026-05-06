# Benchmark Experiments

The four experiments in Table 2 of Moreno-Parada et al. (2024) are designed to
test different physical regimes of the solution.  They share the same
dimensional anchors (chosen to give round dimensionless numbers) but isolate
individual effects.

## Shared dimensional anchors

| Parameter | Value | Description |
|-----------|-------|-------------|
| ``L`` | 1000 m | Ice thickness |
| ``T_\text{air}`` | 250 K | Surface temperature |
| ``\kappa`` | ``10^{-6}`` m² s⁻¹ | Thermal diffusivity |
| ``k`` | 2.0 W m⁻¹ K⁻¹ | Thermal conductivity |

## Parameter table

| Experiment | Pe | Br | ``\gamma`` | ``\beta'`` | ``\Lambda`` | Key physics |
|------------|----|----|------------|------------|-------------|-------------|
| [Exp 1](@ref "Exp 1 — Pure diffusion") | 0  | 0 | 2 | 0 | 0 | Pure diffusion, linear equilibrium |
| [Exp 2](@ref "Exp 2 — Advection")      | 5  | 0 | 2 | 0 | 0 | Upward advection (ablation zone) |
| [Exp 3](@ref "Exp 3 — Strain heating") | 0  | 6 | 2 | 0 | 0 | Strong strain heating |
| [Exp 4](@ref "Exp 4 — Full case")      | 5  | 0 | 2 | 1 | 3 | All effects combined |

## Diffusion timescale

The characteristic diffusion time for these anchors is

```math
t_\text{diff} = \frac{L^2}{\kappa \cdot \text{yr}} \approx 31\,700\;\text{yr} .
```

One dimensionless time unit ``\tau = 1`` corresponds to roughly 31 700 yr.
High-mode transients decay as ``e^{\lambda_n \tau}``; the slowest mode (n=1)
has a decay time of ``1/|\lambda_1|`` dimensionless units — see each experiment
page for the specific value.

## Accessing the benchmarks

Each experiment's `IceColumnPar` is obtained with:

```julia
using IceColumnSolutions
par = benchmark(:exp1)   # :exp1, :exp2, :exp3, or :exp4
```

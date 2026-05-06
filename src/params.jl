# seconds per Julian year — used for all time-unit conversions
const YR_TO_S = 365.25 * 24 * 3600

"""
    IceColumnPar

Parameters for the 1D advective–diffusive ice-column heat equation
(Moreno-Parada et al., 2024, doi:10.5194/tc-18-4215-2024).

Holds both dimensional and dimensionless (non-dimensional) fields.
Construct with physical parameters (primary) or dimensionless parameters.

# Dimensional fields
- `L`       : ice thickness (m)
- `T_air`   : air temperature (K)
- `kappa`   : thermal diffusivity (m² s⁻¹)
- `k`       : thermal conductivity (W m⁻¹ K⁻¹)
- `beta`    : surface thermal insulation length scale (m); 0 = perfect conductor
- `G`       : geothermal heat flux (W m⁻²)
- `Q`       : basal frictional heat flux (W m⁻²)
- `w0`      : surface vertical velocity (m yr⁻¹); negative = downward (accumulation)
- `S`       : depth-uniform strain heating rate (W m⁻³)
- `H`       : depth-averaged horizontal advection (K yr⁻¹)
- `p`       : velocity exponent for power-law profile w(z)=w0*(z/L)^p; 1 = linear

# Dimensionless fields (Table 1, Moreno-Parada et al. 2024)
- `Pe`         : Péclet number  = w0*L / (kappa*YR_TO_S)
- `Br`         : Brinkman number = L²*S / (k*T_air)
- `gamma`      : basal heat parameter = -(G+Q)*L / (k*T_air)
- `beta_prime` : dimensionless surface insulation = beta/L
- `Lambda`     : dimensionless horizontal advection = L²*H / (kappa_yr*T_air)
"""
struct IceColumnPar
    # --- dimensional ---
    L::Float64
    T_air::Float64
    kappa::Float64
    k::Float64
    beta::Float64
    G::Float64
    Q::Float64
    w0::Float64
    S::Float64
    H::Float64
    p::Float64
    # --- dimensionless ---
    Pe::Float64
    Br::Float64
    gamma::Float64
    beta_prime::Float64
    Lambda::Float64
end

"""
    IceColumnPar(L, T_air, kappa, k, beta, G; w0=0.0, Q=0.0, S=0.0, H=0.0, p=1.0)

Construct from **physical** parameters. Dimensionless fields are computed automatically.

# Arguments
- `L`     : ice thickness (m)
- `T_air` : air temperature (K)
- `kappa` : thermal diffusivity (m² s⁻¹)
- `k`     : thermal conductivity (W m⁻¹ K⁻¹)
- `beta`  : surface insulation (m); β=0 → prescribed surface temperature
- `G`     : geothermal heat flux (W m⁻²)
- `w0`    : surface vertical velocity (m yr⁻¹); default 0 (pure diffusion)
- `Q`     : basal frictional heat (W m⁻²); default 0
- `S`     : strain heating rate (W m⁻³); default 0
- `H`     : horizontal advection rate (K yr⁻¹); default 0
- `p`     : velocity exponent; default 1 (linear)
"""
function IceColumnPar(L, T_air, kappa, k, beta, G;
                      w0=0.0, Q=0.0, S=0.0, H=0.0, p=1.0)
    kappa_yr = kappa * YR_TO_S          # m² yr⁻¹
    Pe         = w0 * L / kappa_yr
    Br         = L^2 * S / (k * T_air)
    gamma      = -(G + Q) * L / (k * T_air)
    beta_prime = beta / L
    Lambda     = L^2 * H / (kappa_yr * T_air)
    IceColumnPar(Float64(L), Float64(T_air), Float64(kappa), Float64(k),
                 Float64(beta), Float64(G), Float64(Q), Float64(w0),
                 Float64(S), Float64(H), Float64(p),
                 Pe, Br, gamma, beta_prime, Lambda)
end

"""
    IceColumnPar(L, T_air, kappa, k, beta_prime, gamma, Pe; Br=0.0, Lambda=0.0, p=1.0)

Construct from **dimensionless** parameters plus the four dimensional anchors
needed to recover physical fields.

# Arguments
- `L`, `T_air`, `kappa`, `k` : dimensional anchors (same units as above)
- `beta_prime` : dimensionless surface insulation β̃ = β/L
- `gamma`      : dimensionless basal heat parameter γ
- `Pe`         : Péclet number
- `Br`         : Brinkman number; default 0
- `Lambda`     : dimensionless horizontal advection; default 0
- `p`          : velocity exponent; default 1

Physical fields are back-calculated where possible; those that cannot be
uniquely recovered (Q, S, H split from their combined dimensionless form)
are set to 0 with G, w0, beta absorbing the full value.
"""
function IceColumnPar(L, T_air, kappa, k, beta_prime, gamma, Pe;
                      Br=0.0, Lambda=0.0, p=1.0)
    kappa_yr   = kappa * YR_TO_S
    beta       = beta_prime * L
    # recover combined physical quantities
    G          = -gamma * k * T_air / L   # Q absorbed into G (Q=0 by convention)
    Q          = 0.0
    w0         = Pe * kappa_yr / L
    S          = Br * k * T_air / L^2
    H          = Lambda * kappa_yr * T_air / L^2
    IceColumnPar(Float64(L), Float64(T_air), Float64(kappa), Float64(k),
                 Float64(beta), Float64(G), Float64(Q), Float64(w0),
                 Float64(S), Float64(H), Float64(p),
                 Float64(Pe), Float64(Br), Float64(gamma),
                 Float64(beta_prime), Float64(Lambda))
end

# Convenience: show compact dimensionless summary
function Base.show(io::IO, p::IceColumnPar)
    print(io, "IceColumnPar(L=$(p.L) m, T_air=$(p.T_air) K | ",
          "Pe=$(round(p.Pe,digits=3)), Br=$(round(p.Br,digits=3)), ",
          "γ=$(round(p.gamma,digits=3)), β'=$(round(p.beta_prime,digits=3)), ",
          "Λ=$(round(p.Lambda,digits=3)))")
end

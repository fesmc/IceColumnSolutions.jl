# Theory

## The governing equation

The temperature ``T(z,t)`` in a one-dimensional vertical ice column of
thickness ``L`` satisfies the advective–diffusive heat equation

```math
\rho c_p \frac{\partial T}{\partial t}
= k \frac{\partial^2 T}{\partial z^2}
- \rho c_p\, w(z)\, \frac{\partial T}{\partial z}
+ S - \rho c_p H ,
\qquad z \in [0, L],\; t > 0
```

where ``z = 0`` is the ice base, ``z = L`` is the surface,
``w(z) = w_0 (z/L)^p`` is the vertical velocity (negative = downward),
``S`` is the strain-heating rate (W m⁻³), and ``H`` is a depth-averaged
horizontal heat advection rate (K yr⁻¹).

### Boundary conditions

| Location | Condition | Physical meaning |
|----------|-----------|-----------------|
| ``z = 0`` (base) | ``k \,\partial_z T = G + Q`` | Prescribed heat flux: geothermal ``G`` + basal friction ``Q`` |
| ``z = L`` (surface) | ``\beta\,\partial_z T + T = T_\text{air}`` | Robin: ``\beta = 0`` gives Dirichlet (prescribed surface temperature), ``\beta > 0`` gives partial insulation |

## Non-dimensionalisation

Introducing ``\theta = T / T_\text{air}``, ``\xi = z/L``,
``\tau = \kappa t / L^2`` (``\kappa = k/(\rho c_p)``), the equation becomes

```math
\frac{\partial\theta}{\partial\tau}
= \frac{\partial^2\theta}{\partial\xi^2}
- \mathrm{Pe}\,\xi\,\frac{\partial\theta}{\partial\xi}
+ \Omega ,
\qquad \xi\in[0,1]
```

with boundary conditions

```math
\frac{\partial\theta}{\partial\xi}\bigg|_{\xi=0} = \gamma ,
\qquad
\beta'\,\frac{\partial\theta}{\partial\xi}\bigg|_{\xi=1}
+ \theta\big|_{\xi=1} = 1 .
```

The five dimensionless parameters are (Table 1 of Moreno-Parada et al. 2024):

| Symbol | Name | Definition | Typical range |
|--------|------|-----------|--------------|
| ``\mathrm{Pe}`` | Péclet number | ``w_0 L / (\kappa \cdot \text{yr})`` | ``[-10,\; 0]`` accumulation; ``[0,\; 10]`` ablation |
| ``\mathrm{Br}`` | Brinkman number | ``L^2 S / (k T_\text{air})`` | ``[0,\; 10]`` |
| ``\gamma`` | Basal heat parameter | ``-(G+Q)L / (k T_\text{air})`` | ``[-1,\; -0.05]`` realistic; ``[0, 5]`` benchmark |
| ``\beta'`` | Surface insulation | ``\beta / L`` | ``0`` (Dirichlet) to ``\sim 1`` |
| ``\Lambda`` | Horizontal advection | ``L^2 H / (\kappa_\text{yr} T_\text{air})`` | ``[-5,\; 5]`` |

and ``\Omega = \mathrm{Br} + \Lambda`` is the total dimensionless heat source.

## Solution decomposition

The total solution is split into a stationary part ``\vartheta(\xi)``
(which satisfies the steady equation with all source terms) and a
transient part ``\mu(\xi,\tau)`` that carries the initial condition:

```math
\theta(\xi,\tau) = \vartheta(\xi) + \mu(\xi,\tau) .
```

## Stationary solution

``\vartheta(\xi)`` satisfies the ODE (Appendix B, Eq. B1)

```math
\vartheta_{\xi\xi} - \mathrm{Pe}\,\xi\,\vartheta_\xi = \Omega ,
\qquad \vartheta_\xi(0) = \gamma ,\quad
\beta'\,\vartheta_\xi(1) + \vartheta(1) = 1 .
```

For ``\mathrm{Pe} \ne 0`` the solution is (Eq. B2, corrected sign for ``A``):

```math
\vartheta(\xi)
= \frac{\Omega\,\xi^2}{2}\;{}_2F_2\!\left(1,1;\tfrac{3}{2},2;\,-a^2\xi^2\right)
+ A\,\mathrm{erf}(a\xi) + B ,
```

where ``a = \sqrt{\mathrm{Pe}/2}`` (imaginary for ``\mathrm{Pe} < 0``),

```math
A = \frac{\gamma\sqrt{\pi}}{2a}, \qquad
B = 1 - A\!\left(\frac{2a\,e^{-a^2}}{\sqrt{\pi}}\,\beta' + \mathrm{erf}(a)\right)
    - \Omega\!\left[\!\left(\beta'+\tfrac{1}{2}\right)F_1 - \tfrac{\beta' a^2}{3}F_2\right],
```

with ``F_1 = {}_2F_2(1,1;\tfrac{3}{2},2;-a^2)`` and
``F_2 = {}_2F_2(2,2;\tfrac{5}{2},3;-a^2)``.
All arithmetic is performed in ``\mathbb{C}``; the physical result is real.

**Special case ``\mathrm{Pe} = 0``** (pure diffusion):

```math
\vartheta(\xi) = -\frac{\Omega\xi^2}{2} + \gamma\xi + C ,\qquad
C = 1 + \left(\beta'+\tfrac{1}{2}\right)\Omega - (\beta'+1)\gamma .
```

## Transient solution

``\mu(\xi,\tau)`` satisfies the homogeneous PDE with homogeneous BCs
(Appendix A, Eq. A7):

```math
\mu(\xi,\tau)
= \sum_{n=1}^{\infty}
  A_n\;M\!\left(\alpha_n;\,\tfrac{1}{2};\,\frac{\mathrm{Pe}\,\xi^2}{2}\right)
  \exp(\lambda_n\,\tau) ,
```

where ``M(a,b,z)`` is the Kummer confluent hypergeometric function
``{}_1F_1(a;b;z)``, and ``\lambda_n = 2\,\mathrm{Pe}\,\alpha_n < 0``.

**Special case ``\mathrm{Pe} = 0``**: eigenfunctions are ``\cos(k_n\xi)``
with ``\lambda_n = -k_n^2``, where ``k_n`` satisfies
``-\beta' k_n \sin k_n + \cos k_n = 0``
(for ``\beta' = 0``: ``k_n = (n-\tfrac{1}{2})\pi``).

### Eigenvalue equation

The parameters ``\alpha_n`` satisfy (Eq. A8):

```math
\beta'\,\mathrm{Pe}\,2\alpha_n\,M\!\left(\alpha_n+1;\,\tfrac{3}{2};\,\frac{\mathrm{Pe}}{2}\right)
+ M\!\left(\alpha_n;\,\tfrac{1}{2};\,\frac{\mathrm{Pe}}{2}\right) = 0 .
```

For decay (``\lambda_n < 0``): ``\alpha_n`` and ``\mathrm{Pe}`` have
**opposite signs** — for ``\mathrm{Pe} > 0`` the eigenvalues are negative,
for ``\mathrm{Pe} < 0`` they are positive.

### Series coefficients

```math
A_n = \frac{\displaystyle\int_0^1
      \bigl[\theta_0(\xi)-\vartheta(\xi)\bigr]\,\varrho(\xi)\,\Phi_n(\xi)\,d\xi}
{\displaystyle\int_0^1 \Phi_n(\xi)^2\,\varrho(\xi)\,d\xi} ,
\qquad \varrho(\xi) = e^{-\mathrm{Pe}\,\xi^2/2} ,
```

computed numerically using adaptive Gauss–Kronrod quadrature.

---

## Known errata in Moreno-Parada et al. (2024)

Two sign errors were found in the published appendices during implementation.
The corrected formulas are used throughout this package.

**Eq. B2 — constant ``A``:**
The paper states ``A = -\gamma\sqrt{\pi/(4a)}``.
The correct formula derived from the base BC ``\vartheta_\xi(0)=\gamma`` is
``A = \gamma\sqrt{\pi}/(2a)``.
The published formula has both the wrong sign and the wrong power of ``a``.

**Eq. A7 — decay rate ``\lambda_n``:**
The paper states ``\lambda_n = -2\,\mathrm{Pe}\,\alpha_n``.
The correct formula from the Kummer ODE
(``X''-\mathrm{Pe}\,\xi\,X' = 2\,\mathrm{Pe}\,\alpha\,X``) is
``\lambda_n = 2\,\mathrm{Pe}\,\alpha_n``.
With the paper's sign the modes for ``\mathrm{Pe}<0`` would grow rather than decay.

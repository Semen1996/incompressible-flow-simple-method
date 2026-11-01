# Lid-driven cavity flow: the SIMPLE method

An Octave/MATLAB program that solves the steady incompressible Navier–Stokes equations in a square cavity whose top lid moves at a constant speed. It uses:

- the SIMPLE algorithm (Patankar, Ch. 6);
- the finite volume method on a staggered (MAC) grid;
- Patankar's generalized convection–diffusion formula with a choice of scheme: power-law, hybrid, upwind or central;
- under-relaxation of velocity and pressure.

The method is derived step by step in the article [docs/physics-en.md](docs/physics-en.md): the equations, discretization, the staggered grid, the SIMPLE algorithm, boundary conditions, linear solvers, under-relaxation and the convergence criterion. Comments in the code cite the article's sections and equations.

## Parameters

All parameters are set at the top of [main.m](main.m):

| Parameter  | Default   | Meaning                                                              |
|------------|-----------|----------------------------------------------------------------------|
| `Re`       | 1000      | Reynolds number $\mathrm{Re} = \rho U L / \mu$                       |
| `nx`, `ny` | 64, 64    | number of cells in $x$ and $y$                                       |
| `scheme`   | `'power'` | scheme: `'power'`, `'hybrid'`, `'upwind'`, `'central'` (Sec. 2.3)    |
| `alpha_u`  | 0.7       | velocity under-relaxation (Sec. 7.2)                                 |
| `alpha_p`  | 0.3       | pressure under-relaxation (Sec. 7.3)                                 |
| `max_iter` | 5000      | maximum number of SIMPLE iterations                                  |
| `tol`      | 1e-6      | tolerance on the normalized residuals (Sec. 7.6)                     |

The cavity size, lid speed and density are all equal to one. If the residuals grow or oscillate, reduce `alpha_u` and `alpha_p`. The rule $\alpha_p \approx 1 - \alpha_u$ usually works well (Sec. 7.5).

## Structure

| File                                   | Purpose                                                              |
|----------------------------------------|----------------------------------------------------------------------|
| [main.m](main.m)                       | parameters, SIMPLE loop, convergence check                           |
| [solve_u.m](solve_u.m), [solve_v.m](solve_v.m) | momentum equations for $u$ and $v$ with the current pressure |
| [solve_pc.m](solve_pc.m)               | pressure correction equation for $p'$                                |
| [coefficients.m](coefficients.m)       | coefficients $a_E, a_W, a_N, a_S$ from fluxes $F$ and conductances $D$ |
| [scheme_A.m](scheme_A.m)               | function $A(\lvert P \rvert)$ for the chosen scheme                  |
| [solve_relaxed.m](solve_relaxed.m)     | residual, under-relaxation and direct solve                          |
| [assemble.m](assemble.m)               | assembly of the sparse five-point matrix                             |
| [plot_results.m](plot_results.m)       | plot of the velocity magnitude                                       |

One SIMPLE iteration:

1. Solve the momentum equations with the current pressure $p^*$ to get $u^*$, $v^*$.
2. Solve the pressure correction equation for $p'$.
3. Correct the pressure with under-relaxation: $p = p^* + \alpha_p\, p'$.
4. Correct the velocities with the full correction $p'$.
5. Check that the normalized momentum and continuity residuals are below `tol`.

## References

1. Patankar S. V. _Numerical Heat Transfer and Fluid Flow._ Hemisphere, 1980.
2. Ferziger J. H., Perić M. _Computational Methods for Fluid Dynamics._ 3rd ed., Springer, 2002.
3. [Simulating fluid with Chorin's method](https://github.com/Semen1996/incompressible-flow-chorin-method) (in Russian): the same class of problems solved with the projection method.

The full reference list is at the end of the [article](docs/physics-en.md).

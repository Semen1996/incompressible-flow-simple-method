# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A 2D lid-driven square cavity solver for the steady incompressible Navier–Stokes equations, written in Octave/MATLAB. It uses the SIMPLE algorithm (Patankar, Ch. 6) with the finite volume method on a staggered (MAC) grid. Convection is handled by Patankar's generalized coefficient formula with a selectable scheme.

`docs/physics.md` is a companion article in Russian that derives every step. Code comments cite its section numbers (`Section 4.3`, `(2.3)`, `(7.1)`). When you change the numerics, keep the code and the cited sections consistent.

## Running

There is no build, test suite or linter. Run the script from the repo root, since the functions are resolved from the current directory:

```
octave --no-gui main.m        # or: run main.m from the Octave/MATLAB prompt
```

All parameters (`Re`, `nx`, `ny`, `scheme`, `alpha_u`, `alpha_p`, `max_iter`, `tol`) are hard-coded at the top of `main.m`. The script prints whether it converged and plots the velocity magnitude.

## Architecture

**Array layout (first index = x, second index = y):**
- `p`: `nx × ny`, at cell centres.
- `u`: `(nx+1) × ny`, on vertical faces. `u(1,:)` and `u(nx+1,:)` lie on the walls; the unknowns are `u(2:nx,:)`.
- `v`: `nx × (ny+1)`, on horizontal faces. `v(:,1)` and `v(:,ny+1)` lie on the walls; the unknowns are `v(:,2:ny)`.

**One SIMPLE iteration in `main.m`:** `solve_u`/`solve_v` (momentum with the current `p`, giving `u*`, `v*` and `d`) → `solve_pc` (pressure correction `p'`) → `p += alpha_p * p'` → velocity correction with the full, unrelaxed `p'`.

**Shared pipeline for each transport equation:**
1. Build the face fluxes `F` and conductances `D` on the variable's own control volume.
2. `coefficients.m` combines them via `scheme_A.m` into `aE, aW, aN, aS` (each the size of the unknown array). Then `aP = Σ a_nb`.
3. Boundary handling is done by zeroing the link to a boundary node and moving its known value into `b`. Wall-adjacent conductances use `2*D`, because the distance to the wall is half a cell.
4. `solve_relaxed.m` computes the unrelaxed residual from the current field first, then applies under-relaxation (`aP/alpha`) and solves directly with `\`.
5. `assemble.m` builds the sparse five-point matrix with linear index `k = i + (j-1)*m` (E = k+1, W = k−1, N = k+m, S = k−m). Coefficients that point outside the array are dropped.

**Conventions to preserve:**
- `d_u = dy ./ aP_relaxed` and `d_v = dx ./ aP_relaxed` use the *relaxed* `aP`. They are zero on wall faces, so `solve_pc` needs no explicit boundary treatment for `p'`.
- The level of `p'` is pinned by replacing the equation in cell `(1,1)` with `p' = 0`.
- Convergence: the momentum residuals are normalized by `rho*U_lid^2*L` and the mass residual by `rho*U_lid*L`. All three must be below `tol`.
- The code is vectorized; there are no loops over cells. Keep new code in that style, with short English comments that cite article sections.

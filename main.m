% MAIN  Lid-driven square cavity flow, SIMPLE method.
%
%   The steady incompressible Navier–Stokes equations are solved with
%   the SIMPLE method (Patankar, Ch. 6) on a staggered grid. Discretization:
%   finite volume method; convection: Patankar's power-law scheme.
%
%   Section numbers in the comments refer to the article docs/physics.md.

clear all; close all;

%% Physical parameters
L     = 1.0;                   % cavity width
H     = 1.0;                   % cavity height
U_lid = 1.0;                   % lid velocity
rho   = 1.0;                   % density
Re    = 1000;                   % Reynolds number Re = rho*U_lid*L/mu
mu    = rho * U_lid * L / Re;  % dynamic viscosity

%% Numerical parameters
nx = 64;                       % number of cells in x
ny = 64;                       % number of cells in y
scheme      = 'power';         % convection scheme: 'power', 'hybrid', 'upwind', 'central' (Section 2.3)
alpha_u     = 0.7;             % velocity under-relaxation (Section 7.2)
alpha_p     = 0.3;             % pressure under-relaxation (Section 7.3)
max_iter    = 5000;            % maximum number of SIMPLE iterations
tol         = 1e-6;            % tolerance on normalized residuals (Section 7.6)

dx = L / nx;
dy = H / ny;

%% Fields on the staggered grid (Section 3.3); first index is x, second is y
p = zeros(nx,   ny);           % pressure at cell centres
u = zeros(nx+1, ny);           % u on vertical faces; u(1,:), u(nx+1,:) are on the walls
v = zeros(nx,   ny+1);         % v on horizontal faces; v(:,1), v(:,ny+1) are on the walls

% Scales for residual normalization (Section 7.6)
scale_mom  = rho * U_lid^2 * L;
scale_mass = rho * U_lid * L;

converged = false;

for iter = 1:max_iter
  % Step 2: momentum equations with the current pressure p* -> u*, v*
  [u_star, d_u, res_u] = solve_u(u, v, p, U_lid, rho, mu, dx, dy, alpha_u, scheme);
  [v_star, d_v, res_v] = solve_v(u, v, p, rho, mu, dx, dy, alpha_u, scheme);

  % Step 3: pressure correction equation for p'
  [pc, res_m] = solve_pc(u_star, v_star, d_u, d_v, rho, dx, dy);

  % Step 4: pressure correction with under-relaxation
  p = p + alpha_p * pc;

  % Step 5: velocity correction with the full correction (Sections 4.5, 7.3)
  u = u_star;
  v = v_star;
  u(2:nx,:) = u(2:nx,:) + d_u(2:nx,:) .* (pc(1:nx-1,:) - pc(2:nx,:));
  v(:,2:ny) = v(:,2:ny) + d_v(:,2:ny) .* (pc(:,1:ny-1) - pc(:,2:ny));

  % Convergence check
  R = [res_u / scale_mom, res_v / scale_mom, res_m / scale_mass];

  if all(R < tol)
    converged = true;
    break;
  end
end

if converged
  fprintf('Converged in %d iterations \n', iter);
else
  fprintf('NOT converged in %d iterations \n', max_iter);
end

plot_results(u, v, p, L, H, U_lid, Re);

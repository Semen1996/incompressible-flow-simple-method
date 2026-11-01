function [v_new, d_v, res] = solve_v(u, v, p, rho, mu, dx, dy, alpha, scheme)
% Momentum equation for v (Section 3.4). Unknowns: v(:, 2:ny).
  nx = size(v, 1);
  ny = size(v, 2) - 1;
  n  = ny - 1;

  % Mass fluxes through the faces of the v control volume
  Fe = rho * dy * 0.5 * (u(2:nx+1,1:ny-1) + u(2:nx+1,2:ny));
  Fw = rho * dy * 0.5 * (u(1:nx,1:ny-1)   + u(1:nx,2:ny));
  Fn = rho * dx * 0.5 * (v(:,2:ny)   + v(:,3:ny+1));
  Fs = rho * dx * 0.5 * (v(:,1:ny-1) + v(:,2:ny));

  % Face conductances; at the side walls the distance to the wall is dx/2 (Section 5.1)
  De = mu * dy / dx * ones(nx, n);
  Dw = De;
  Dn = mu * dx / dy * ones(nx, n);
  Ds = Dn;
  De(nx,:) = 2 * mu * dy / dx;
  Dw(1,:)  = 2 * mu * dy / dx;

  [aE, aW, aN, aS] = coefficients(Fe, Fw, Fn, Fs, De, Dw, Dn, Ds, scheme);
  aP = aE + aW + aN + aS;

  % Pressure force: (p_P - p_N) * A_n
  b = dx * (p(:,1:ny-1) - p(:,2:ny));

  % v = 0 on all walls: links to boundary nodes are removed, no contribution to b
  aE(nx,:) = 0;   % right wall
  aW(1,:)  = 0;   % left wall
  aN(:,n)  = 0;   % lid
  aS(:,1)  = 0;   % bottom

  [v_in, res, aP_relaxed] = solve_relaxed(aP, aE, aW, aN, aS, b, v(:,2:ny), alpha);

  v_new = v;
  v_new(:,2:ny) = v_in;

  d_v = zeros(size(v));
  d_v(:,2:ny) = dx ./ aP_relaxed;
end

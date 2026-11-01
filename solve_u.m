function [u_new, d_u, res] = solve_u(u, v, p, U_lid, rho, mu, dx, dy, alpha, scheme)
% Momentum equation for u (Section 3.4). Unknowns: u(2:nx, :).
  nx = size(u, 1) - 1;
  ny = size(u, 2);
  m  = nx - 1;

  % Mass fluxes through the faces of the u control volume
  Fe = rho * dy * 0.5 * (u(2:nx,:)   + u(3:nx+1,:));
  Fw = rho * dy * 0.5 * (u(1:nx-1,:) + u(2:nx,:));
  Fn = rho * dx * 0.5 * (v(1:nx-1,2:ny+1) + v(2:nx,2:ny+1));
  Fs = rho * dx * 0.5 * (v(1:nx-1,1:ny)   + v(2:nx,1:ny));

  % Face conductances; at the bottom and the lid the distance to the wall is dy/2 (Section 5.1)
  De = mu * dy / dx * ones(m, ny);
  Dw = De;
  Dn = mu * dx / dy * ones(m, ny);
  Ds = Dn;
  Dn(:,ny) = 2 * mu * dx / dy;
  Ds(:,1)  = 2 * mu * dx / dy;

  [aE, aW, aN, aS] = coefficients(Fe, Fw, Fn, Fs, De, Dw, Dn, Ds, scheme);
  aP = aE + aW + aN + aS;

  % Pressure force: (p_P - p_E) * A_e
  b = dy * (p(1:nx-1,:) - p(2:nx,:));

  % Known boundary values move to the right-hand side
  b(m,:)  = b(m,:)  + aE(m,:) .* u(nx+1,:);   aE(m,:)  = 0;   % right wall
  b(1,:)  = b(1,:)  + aW(1,:) .* u(1,:);      aW(1,:)  = 0;   % left wall
  b(:,ny) = b(:,ny) + aN(:,ny) * U_lid;       aN(:,ny) = 0;   % lid
  aS(:,1) = 0;                                                % bottom: u = 0, no contribution to b

  [u_in, res, aP_relaxed] = solve_relaxed(aP, aE, aW, aN, aS, b, u(2:nx,:), alpha);

  u_new = u;
  u_new(2:nx,:) = u_in;

  % d = A_e / (a_e / alpha) for the velocity correction formula (Sections 4.2, 7.2)
  d_u = zeros(size(u));
  d_u(2:nx,:) = dy ./ aP_relaxed;
end

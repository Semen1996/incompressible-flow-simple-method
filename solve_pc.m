function [pc, res] = solve_pc(u, v, d_u, d_v, rho, dx, dy)
% Pressure correction equation (Sections 4.3, 4.4).
% d_u, d_v are zero on the walls, so the link to the outside vanishes by itself (Section 5.3).
  nx = size(u, 1) - 1;
  ny = size(v, 2) - 1;

  aE = rho * dy * d_u(2:nx+1,:);
  aW = rho * dy * d_u(1:nx,:);
  aN = rho * dx * d_v(:,2:ny+1);
  aS = rho * dx * d_v(:,1:ny);
  aP = aE + aW + aN + aS;

  % Mass imbalance of the u*, v* field
  b = rho * dy * (u(1:nx,:) - u(2:nx+1,:)) + rho * dx * (v(:,1:ny) - v(:,2:ny+1));
  res = sum(abs(b(:)));

  % The level of p' is fixed in cell (1,1): its equation is replaced by p' = 0 (Section 5.4)
  aP(1,1) = 1;
  aE(1,1) = 0;
  aN(1,1) = 0;
  b(1,1)  = 0;

  pc = reshape(assemble(aP, aE, aW, aN, aS) \ b(:), nx, ny);
end

function plot_results(u, v, p, L, H, U_lid, Re)
% Velocity magnitude field
  nx = size(p, 1);
  ny = size(p, 2);
  dx = L / nx;
  dy = H / ny;
  xc = ((1:nx) - 0.5) * dx;      % cell centres
  yc = ((1:ny) - 0.5) * dy;

  % Velocity at cell centres
  uc = 0.5 * (u(1:nx,:) + u(2:nx+1,:));
  vc = 0.5 * (v(:,1:ny) + v(:,2:ny+1));
  speed = sqrt(uc.^2 + vc.^2);

  % Colormap blue -> green -> red
  cmap = interp1([0; 0.5; 1], [0 0 1; 0 1 0; 1 0 0], linspace(0, 1, 256)');

  figure(1); clf;
  imagesc(xc, yc, speed');
  axis xy equal tight;
  colormap(cmap);
  colorbar;
  xlabel('x'); ylabel('y');
  title(sprintf('|u|, Re = %g', Re));
end

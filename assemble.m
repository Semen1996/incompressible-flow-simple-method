function A = assemble(aP, aE, aW, aN, aS)
% Sparse matrix of the five-point equation
%   aP*phi_P - aE*phi_E - aW*phi_W - aN*phi_N - aS*phi_S = b.
% Node (i,j) of an m x n array gets index k = i + (j-1)*m;
% its neighbours: E = k+1, W = k-1, N = k+m, S = k-m (Section 6.2).
  [m, n] = size(aP);
  k  = reshape(1:m*n, m, n);
  kE = k(1:m-1,:);
  kW = k(2:m,:);
  kN = k(:,1:n-1);
  kS = k(:,2:n);

  rows = [k(:); kE(:);     kW(:);     kN(:);     kS(:)];
  cols = [k(:); kE(:) + 1; kW(:) - 1; kN(:) + m; kS(:) - m];
  vals = [aP(:);
          -reshape(aE(1:m-1,:), [], 1);
          -reshape(aW(2:m,:),   [], 1);
          -reshape(aN(:,1:n-1), [], 1);
          -reshape(aS(:,2:n),   [], 1)];

  A = sparse(rows, cols, vals, m*n, m*n);
end

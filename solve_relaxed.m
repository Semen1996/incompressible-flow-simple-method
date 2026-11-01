function [x, res, aP_relaxed] = solve_relaxed(aP, aE, aW, aN, aS, b, x0, alpha)
% Solves a momentum equation with under-relaxation (7.1).
% The residual is computed before solving, from the current field x0 and without relaxation (Section 7.6).
  res = sum(abs(b(:) - assemble(aP, aE, aW, aN, aS) * x0(:)));

  aP_relaxed = aP / alpha;
  b_relaxed  = b + (1 - alpha) * aP_relaxed .* x0;
  x = reshape(assemble(aP_relaxed, aE, aW, aN, aS) \ b_relaxed(:), size(b));
end

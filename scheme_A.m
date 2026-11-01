function A = scheme_A(P, scheme)
% Function A(|P|) of the cell Peclet number (table in Section 2.3)
  P = abs(P);
  switch scheme
    case 'central'
      A = 1 - 0.5 * P;
    case 'upwind'
      A = ones(size(P));
    case 'hybrid'
      A = max(0, 1 - 0.5 * P);
    case 'power'
      A = max(0, (1 - 0.1 * P).^5);
    otherwise
      error('Unknown scheme: %s', scheme);
  end
end

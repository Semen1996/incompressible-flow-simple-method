function [aE, aW, aN, aS] = coefficients(Fe, Fw, Fn, Fs, De, Dw, Dn, Ds, scheme)
% Discretization coefficients from Patankar's generalized formula (2.3)
  aE = De .* scheme_A(Fe ./ De, scheme) + max(-Fe, 0);
  aW = Dw .* scheme_A(Fw ./ Dw, scheme) + max( Fw, 0);
  aN = Dn .* scheme_A(Fn ./ Dn, scheme) + max(-Fn, 0);
  aS = Ds .* scheme_A(Fs ./ Ds, scheme) + max( Fs, 0);
end

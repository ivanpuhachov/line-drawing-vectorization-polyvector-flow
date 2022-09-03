function [ energy, grad ] = polynomial_energy_2019( X, W, tau, mask, splitReIm)
constraints = tau(mask(:));
ignoreTerms = (abs(constraints) < 1e-6);
%ignoreTerms = [];
if (size(X,2)~=2)
   X = X(1:end/2)+1i*X(end/2+1:end);
   X = reshape(X,[],2);
end
nonzeros = nnz(mask);
Xext = [X(:,1) X(:,2) ones(nonzeros,1)];
tauMonomial = [ones(nonzeros,1) tau(mask(:)).^2 tau(mask(:)).^4];
polyResidual = sum(Xext.*tauMonomial,2);
polyResidual(ignoreTerms) = 0;
energy = sum((abs(polyResidual).^2).*W(mask(:)));

if nargout <= 1 % energy computed, so stop if they didn't ask for gradient
    return;
end

r = real(tau(mask(:)));
m = imag(tau(mask(:)));

a1 = real(X(:,1)); b1 = imag(X(:,1));
a2 = real(X(:,2)); b2 = imag(X(:,2));

a1(ignoreTerms) = 0;
a2(ignoreTerms) = 0;
b1(ignoreTerms) = 0;
b2(ignoreTerms) = 0;

%partials over re(x_i) and im(x_i)
DfDa1 = 2.*a1+(-2).*a2.*m.^2+2.*m.^4+(-4).*b2.*m.*r+2.*a2.*r.^2+(-12).* ...
  m.^2.*r.^2+2.*r.^4;

DfDa2 = (-2).*a1.*m.^2+2.*a2.*m.^4+(-2).*m.^6+4.*b1.*m.*r+2.*a1.*r.^2+4.* ...
  a2.*m.^2.*r.^2+(-2).*m.^4.*r.^2+2.*a2.*r.^4+2.*m.^2.*r.^4+2.*r.^6;

DfDb1 = 2.*b1+(-2).*b2.*m.^2+4.*a2.*m.*r+(-8).*m.^3.*r+2.*b2.*r.^2+8.*m.*r.^3;

DfDb2 = (-2).*b1.*m.^2+2.*b2.*m.^4+(-4).*a1.*m.*r+4.*m.^5.*r+2.*b1.*r.^2+ ...
  4.*b2.*m.^2.*r.^2+8.*m.^3.*r.^3+2.*b2.*r.^4+4.*m.*r.^5;

%Complex Derivatives (function is real => this formula is correct)
DfDx1 = W(mask(:)).*(DfDa1 + 1i*DfDb1);
DfDx2 = W(mask(:)).*(DfDa2 + 1i*DfDb2);

%DfDx1 = zeros(size(DfDx1));
%assert(norm(DfDx1(ignoreTerms))==0);
%assert(norm(DfDx2(ignoreTerms))==0);

%now assemble into the gradient vector
grad = transpose([DfDx1; DfDx2]);
if (splitReIm)
    grad = [real(grad) imag(grad)];
end
end


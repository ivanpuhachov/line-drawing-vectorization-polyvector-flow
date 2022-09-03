function [ energy, grad,outTmp ] = totalEnergy_2019(y, D, tau, m,n,beta,mask,noisy,tmp)
outTmp = [];
y_complex = y(1:end/2)+1i*y(end/2+1:end);
x1_thinband = y_complex(1:end/2);
x2_thinband = y_complex(end/2+1:end);
x_thinband = [x1_thinband x2_thinband];
x = [real(y(1:end/2)+1i*y(end/2+1:end)); imag(y(1:end/2)+1i*y(end/2+1:end))];

Daugm = ones(m+2,n+2);
Daugm(2:end-1,2:end-1) = 1-D;
smartWeights = (circshift(Daugm,1,1)+circshift(Daugm,1,2)+circshift(Daugm,-1,1)+circshift(Daugm,-1,2))/4;
smartWeights = smartWeights(2:end-1,2:end-1);
if (noisy)
    alpha = 0.1;
else
    alpha = 0.01;
end

if (nargout>1)
    [e1,g1] = polynomial_energy_2019(x_thinband,D,tau,mask,true);
    [e2,g2] = polynomial_energy_2019(x_thinband,ones(size(D)),exp(1i*pi/2).*tau,mask,true);
    [e3,g3] = l2_regularizer_2019(x,m,n,smartWeights,true,mask);
    energy = e1 + alpha*e2 + beta*e3;
    grad = g1 + alpha*g2 + beta*g3;
    grad = transpose(grad);
    if (nargin == 9)
        outTmp = [tmp; e1 e2 e3];
    end
end
end


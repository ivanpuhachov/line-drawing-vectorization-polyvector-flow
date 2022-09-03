function [fVec, df_dnu, df_dx] = anisotropicFrameFieldFun(nu, x, field,m,n)
field = reshape(field,[],2);
field1 = reshape(field(:,1),m,n);
field2 = reshape(field(:,2),m,n);
nu_cmplx = nu(:,1)+1i*nu(:,2);
tau = 1i*nu_cmplx;

x_int = [round(x(:,1)),round(x(:,2))];
idx=sub2ind([m n],x_int(:,1),x_int(:,2));
f1 = field1(idx);
f2 = field2(idx);
f1R = real(f1); f1I = imag(f1);
f2R = real(f2); f2I = imag(f2);
fVec = abs(f1+f2.*tau.^2+tau.^4).^2;
nu1 = nu(:,1); nu2 = nu(:,2);
if nargout > 1
    df_dnu = [4.*f2I.^2.*nu1.^3 + 4.*f2I.^2.*nu1.*nu2.^2 + 20.*f2I.*nu1.^4.*nu2 + 24.*f2I.*nu1.^2.*nu2.^3 + 4.*f1I.*f2I.*nu1 + 4.*f2I.*nu2.^5 - 4.*f1R.*f2I.*nu2 + 4.*f2R.^2.*nu1.^3 + 4.*f2R.^2.*nu1.*nu2.^2 + 12.*f2R.*nu1.^5 + 8.*f2R.*nu1.^3.*nu2.^2 - 4.*f2R.*nu1.*nu2.^4 + 4.*f1R.*f2R.*nu1 + 4.*f1I.*f2R.*nu2 + 8.*nu1.^7 + 24.*nu1.^5.*nu2.^2 + 24.*nu1.^3.*nu2.^4 + 8.*f1R.*nu1.^3 + 24.*f1I.*nu1.^2.*nu2 + 8.*nu1.*nu2.^6 - 24.*f1R.*nu1.*nu2.^2 - 8.*f1I.*nu2.^3 + nu1.*((4.*nu1.^3.*nu2 + f2I.*nu1.^2 - 4.*nu1.*nu2.^3 + 2.*f2R.*nu1.*nu2 - f2I.*nu2.^2 + f1I).^2 + (nu1.^4 - 6.*nu1.^2.*nu2.^2 + f2R.*nu1.^2 - 2.*f2I.*nu1.*nu2 + nu2.^4 - f2R.*nu2.^2 + f1R).^2),...
     ((4.*nu1.^3.*nu2 + f2I.*nu1.^2 - 4.*nu1.*nu2.^3 + 2.*f2R.*nu1.*nu2 - f2I.*nu2.^2 + f1I).^2 + (nu1.^4 - 6.*nu1.^2.*nu2.^2 + f2R.*nu1.^2 - 2.*f2I.*nu1.*nu2 + nu2.^4 - f2R.*nu2.^2 + f1R).^2)./1.^(1./2) + (4.*f2I.^2.*nu1.^2.*nu2 + 4.*f2I.^2.*nu2.^3 + 4.*f2I.*nu1.^5 + 24.*f2I.*nu1.^3.*nu2.^2 + 20.*f2I.*nu1.*nu2.^4 - 4.*f1R.*f2I.*nu1 - 4.*f1I.*f2I.*nu2 + 4.*f2R.^2.*nu1.^2.*nu2 + 4.*f2R.^2.*nu2.^3 + 4.*f2R.*nu1.^4.*nu2 - 8.*f2R.*nu1.^2.*nu2.^3 + 4.*f1I.*f2R.*nu1 - 12.*f2R.*nu2.^5 - 4.*f1R.*f2R.*nu2 + 8.*nu1.^6.*nu2 + 24.*nu1.^4.*nu2.^3 + 8.*f1I.*nu1.^3 + 24.*nu1.^2.*nu2.^5 - 24.*f1R.*nu1.^2.*nu2 - 24.*f1I.*nu1.*nu2.^2 + 8.*nu2.^7 + 8.*f1R.*nu2.^3)];   

    %df/dF1R, df/dF1I, df/dF2R, df/dF2I
    df_dfield = [(nu1.^2 + 2.*nu2).^(1./2).*(2.*nu1.^4 - 12.*nu1.^2.*nu2.^2 + 2.*f2R.*nu1.^2 - 4.*f2I.*nu1.*nu2 + 2.*nu2.^4 - 2.*f2R.*nu2.^2 + 2.*f1R),...
                 (nu1.^2 + 2.*nu2).^(1./2).*(8.*nu1.^3.*nu2 + 2.*f2I.*nu1.^2 - 8.*nu1.*nu2.^3 + 4.*f2R.*nu1.*nu2 - 2.*f2I.*nu2.^2 + 2.*f1I),...
                 (nu1.^2 + 2.*nu2).^(1./2).*(2.*nu1.^6 + 2.*nu1.^4.*nu2.^2 + 2.*f2R.*nu1.^4 - 2.*nu1.^2.*nu2.^4 + 4.*f2R.*nu1.^2.*nu2.^2 + 2.*f1R.*nu1.^2 + 4.*f1I.*nu1.*nu2 - 2.*nu2.^6 + 2.*f2R.*nu2.^4 - 2.*f1R.*nu2.^2),...
                 (nu1.^2 + 2.*nu2).^(1./2).*(4.*nu1.^5.*nu2 + 2.*f2I.*nu1.^4 + 8.*nu1.^3.*nu2.^3 + 4.*f2I.*nu1.^2.*nu2.^2 + 2.*f1I.*nu1.^2 + 4.*nu1.*nu2.^5 - 4.*f1R.*nu1.*nu2 + 2.*f2I.*nu2.^4 - 2.*f1I.*nu2.^2)];
    idxXP = sub2ind([m n],x_int(:,1)+1,x_int(:,2));
    idxXM = sub2ind([m n],x_int(:,1)-1,x_int(:,2));
    idxYP = sub2ind([m n],x_int(:,1),x_int(:,2)+1);
    idxYM = sub2ind([m n],x_int(:,1),x_int(:,2)-1);
    dfield1_dx = [(field1(idxXP)-field1(idxXM))/2, (field1(idxYP)-field1(idxYM))/2];
    dfield2_dx = [(field2(idxXP)-field2(idxXM))/2, (field2(idxYP)-field2(idxYM))/2];
    df_dx = [df_dfield(:,1).*real(dfield1_dx(:,1))+df_dfield(:,2).*imag(dfield1_dx(:,1))+df_dfield(:,3).*real(dfield2_dx(:,1))+df_dfield(:,4).*imag(dfield2_dx(:,1)),...
        df_dfield(:,1).*real(dfield1_dx(:,2))+df_dfield(:,2).*imag(dfield1_dx(:,2))+df_dfield(:,3).*real(dfield2_dx(:,2))+df_dfield(:,4).*imag(dfield2_dx(:,2))];
end
end


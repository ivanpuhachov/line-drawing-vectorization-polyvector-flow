function [ energy,grad ] = l2_regularizer_2019(X, m,n, D,splitReIm,mask)
if (size(X,2)~=2)
   X = X(1:end/2)+1i*X(end/2+1:end);
   X = reshape(X,[],2);
end

energy = 0;
myGrad = [];

for k=1:2
    Xreshape = inf*ones(m,n);
    Xreshape(mask) = X(:,k);
    %Xreshape(~mask) = inf;
    DxX = Dx(Xreshape,mask);
    %DxY = transpose(Dx(transpose(Xreshape),transpose(mask)));
    DxY = Dy(Xreshape,mask);
   
    myDx = DxX(mask).*sqrt(D(mask));
    myDy = DxY(mask).*sqrt(D(mask));
        
    energy = energy+ norm(myDx,'fro')^2+norm(myDy,'fro')^2;
    
    %if (i-1) is 0 (unmasked), then add -2D_i
    %if (i-1) is 1, but (i-2) is 0, then add 2D_{i-1}
    %if (i-1) is and (i-2) is 1, then add D_{i-1} 
    
    grad_mine = zeros(m,n);
    dxmat{1} = DxX;
    dxmat{2} = DxY;
    for dir=1:2
        for leftRight = 1:2
            sign = 2*leftRight-3;
            weightShr = zeros(m+2,n+2);
            weightShr(2:end-1,2:end-1) = D;
            weightShr = circshift(weightShr,sign,3-dir);
            weightShr = weightShr(2:end-1,2:end-1);
            
            DxXshr = circshift(dxmat{dir},sign,3-dir);
            maskX = false(m+4,n+4);
            maskX(3:end-2,3:end-2) = mask;
            
            maskXshr = circshift(maskX,1*sign,3-dir);
            maskXshrr = circshift(maskX,2*sign,3-dir);
            
            maskR{1} = maskX & ~maskXshr;
            maskR{2} = maskX & maskXshr & ~maskXshrr;
            maskR{3} = maskX & maskXshr & maskXshrr;
            
            for i=1:3
                maskR{i} = maskR{i}(3:end-2,3:end-2);
             end

            grad_mine(maskR{1}) = grad_mine(maskR{1}) - sign*2*dxmat{dir}(maskR{1}).*D(maskR{1});
            grad_mine(maskR{2}) = grad_mine(maskR{2}) + sign*2*DxXshr(maskR{2}).*weightShr(maskR{2});
            grad_mine(maskR{3}) = grad_mine(maskR{3}) + sign*DxXshr(maskR{3}).*weightShr(maskR{3});
        end
    end
    grad_unwrapped = grad_mine(mask);
    myGrad = [myGrad; grad_unwrapped(:)];
end

grad = transpose(myGrad);

if (splitReIm)
    grad = [real(grad) imag(grad)];
end

end


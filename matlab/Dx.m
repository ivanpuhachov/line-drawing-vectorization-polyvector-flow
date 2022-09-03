function [ D ] = Dx( X, mask )
if (nargin == 1)
    mask = true(size(X));
end

 shr = circshift(X,1,2);
 shl = circshift(X,-1,2);
D = zeros*ones(size(X));
maskX = false(size(X,1),size(X,2)+2);
maskX(:,2:end-1) = mask;
maskShr = circshift(maskX,1,2);
maskShl = circshift(maskX,-1,2);
allThree = maskX & maskShr & maskShl;
meAndLeft = maskX & ~maskShr & maskShl;
meAndRight = maskX & maskShr & ~maskShl;
allThree = allThree(:,2:end-1);
meAndLeft = meAndLeft(:,2:end-1);
meAndRight = meAndRight(:,2:end-1);
D(allThree) = 0.5*(shl(allThree)-shr(allThree));
D(meAndLeft) = shl(meAndLeft)-X(meAndLeft);
D(meAndRight) = X(meAndRight)-shr(meAndRight);
end


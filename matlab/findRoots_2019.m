function [ roots_out ] = findRoots_2019( X, m,n )
nPixels = m*n;
if (numel(X)==nPixels*4)
    X = X(1:2*nPixels)+1i*X(2*nPixels+1:end);
end

X = reshape(X,[],2);
poly_roots = [];

for i=1:nPixels
    poly_roots = [poly_roots findAndSortRoots_2019(X(i,:))];
end

roots_out = cell(1,4);
roots_out{1} = reshape(poly_roots(1,:),m,n);
roots_out{2} = reshape(poly_roots(2,:),m,n);
roots_out{3} = reshape(poly_roots(3,:),m,n);
roots_out{4} = reshape(poly_roots(4,:),m,n);
end


function [ ] = drawNPolyDeg4_2019(xGrid, yGrid, X, roots_out, drawSubset, tau)
if (nargin==4)
    drawSubset = false;
    tau = [];
elseif (nargin==5)
    tau = [];
end

tfHold = ishold;

m = size(xGrid,1);
n = size(xGrid,2);
nPixels = m*n;
sparse = false;
%if (max(m,n) > 500)
%    sparse = true;
%end

if (numel(X)==nPixels*4)
    X = X(1:2*nPixels)+1i*X(2*nPixels+1:end);
end

if (sparse && (nargin>5))
     roots_out{1} = roots_out{1}(tau~=0);
     roots_out{2} = roots_out{2}(tau~=0);
     roots_out{3} = roots_out{3}(tau~=0);
     roots_out{4} = roots_out{4}(tau~=0);
     xGrid = xGrid(tau~=0);
     yGrid = yGrid(tau~=0);
end

%figure;
hold on;

q = cell(1,4);
scale = .5;
%if (sparse)
    scale = 0.2;
%end

if (drawSubset)
    subset = find(abs(roots_out{1})>1e-6);
    subset = subset(1:2:numel(subset));
else
    subset = 1:numel(xGrid);
    subset = find(abs(roots_out{1})>1e-6);
end

xGrid = xGrid(subset);
yGrid = yGrid(subset);

for i=1:4
    %nn = abs(roots_out{i});
    %nn(abs(nn)<1e-6) = 1;
    %roots_out{i} = roots_out{i}./nn;
    roots_out{i} =  roots_out{i}(subset);
q{i} = quiver(xGrid,yGrid,real(roots_out{i}),imag(roots_out{i}),scale,'Color',[0.85 0.32 0.09],'ShowArrowHead','off','linewidth',2);
%q{i} = quiver(xGrid,yGrid,real(roots_out{i}),imag(roots_out{i}));
end

%// Get the current colormap
currentColormap = colormap(autumn); %bone?

if (~isempty(X))
    %// Now determine the color to make each arrow using a colormap
    X = reshape(X,[],2);
    mags = abs(X(:,1));
    if (sparse && (nargin>4))
        mags = mags(tau~=0);
    end
    
    mags = max(mags)-mags;
    mags = mags(subset);
    [~, ~, ind] = histcounts(mags, size(currentColormap, 1));
    ind = uint8(ind*0.7);
    %// Now map this to a colormap to get RGB
    cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
    cmap(:,:,4) = 255;
    cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);
    L = reshape(cmap(1:3,:,:), [], 4).';
    
    for i=1:4
        %set(q{i}.Head,'ColorBinding', 'interpolated', ...
        %'ColorData', L);
        set(q{i}.Tail, ...
            'ColorBinding', 'interpolated', ...
            'ColorData', reshape(cmap(1:2,:,:), [], 4).');
        
        if sparse
            set(q{i},'linewidth',0.01);
        end
    end
    
    %if (nargin>4)&&(~sparse)
    %    quiver(xGrid,yGrid,real(tau),imag(tau),'Color',[0.5 0.5 0.5],'ShowArrowHead','off');
    %end
    
end
axis equal
colormap(gray);
if tfHold
    hold on;
else
    hold off;
end
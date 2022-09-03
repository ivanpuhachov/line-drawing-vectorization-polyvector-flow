%exit
close all; clear;
img = imread('tapir_medium.png');
%img = imread('my_inputs/stool_1_microcut.png');
%img = imread('my_inputs/stool_1_microcut2.png');
%img = imread('my_inputs/cartoon-elephant07_ear.png');
%img = imread('my_inputs/fig2.jpg');

tic;
bwImg = 255-sum(img,3)/3;
m = size(bwImg,1);
n = size(bwImg,2);

%for tapir
%pointA = [21,31];
%pointB = [34,77];
%pointB = [26,53];
pointA = [18, 24];
pointB = [22, 30];

%pointA = [40, 25];
%pointB = [51, 31];
%pointB = [48, 75];

noisy = false;
origMask = true(m,n);
extMask = origMask;
[gMag, gDir] = imgradient(bwImg);
alpha = 0; %dot smoothness
beta = 10; %L2 regularization
tau = exp(1i*(gDir*pi/180+pi/2));
g = exp(1i*(gDir*pi/180));
tau(gMag < max(max(gMag))./10) = 0; 
gMag(gMag < max(max(gMag))/10) = 0;

[x_grid,y_grid] = meshgrid(1:size(bwImg,2), 1:size(bwImg,1));
nPixels = numel(bwImg);

%% OPTIMIZATION
mse = conv2((tau.*gMag).^2,[1 1 1; 1 0 1; 1 1 1],'same');
%div = conv2(gMag.^2,[1 1 1; 1 0 1; 1 1 1],'same');
div = abs(mse);
div(div==0) = 1;
mse = mse./div;
mse = mse-tau.^2;
mse(gMag==0) = 0;
weight = abs(mse);
weight = 1-weight/max(max(weight));
weight(gMag==0) = 0;
X = [-ones(1,nPixels); zeros(1,nPixels)].'; %initialize coefficients to Ts
rng(1);
X = zeros(2*nnz(extMask),1);
fun = @(y)totalEnergy_2019(y,weight,tau,m,n,beta,extMask,noisy);
fun2 = @(y,outTmp)totalEnergy_2019(y,weight,tau,m,n,beta,extMask,noisy,outTmp);
options = [];
%options.display = 'full';
options.maxFunEvals = 1e6;
options.Method = 'lbfgs';
options.maxIter = 2000;
options.optTol = 1e-6;
X0 = [real(X(:)); imag(X(:))];
X_new_narrowband = minFunc(fun,X0,options);

toc

X_new1 = zeros(m,n);
X_new2 = X_new1;
X_new1(extMask) = X_new_narrowband(1:end/4)+X_new_narrowband(end/2+1:end*3/4)*1i;
X_new2(extMask) = X_new_narrowband(end/4+1:end/2)+1i*X_new_narrowband(end*3/4+1:end);

X_new = [real(X_new1(:)); real(X_new2(:)); imag(X_new1(:)); imag(X_new2(:))];
X_cmplx = X_new(1:2*nPixels)+1i*X_new(2*nPixels+1:end);
roots_out = findRoots_2019(X_new,m,n);

%%
figure;
colormap(gray);
imagesc([1 n],[-1 -m],255*3-bwImg);
set(gca,'YDir','normal');
axis equal
hold on;

drawNPolyDeg4_2019(x_grid,-y_grid,X_new,roots_out,false)
%%
%create a graph from the narrow band
origMask = bwImg > 10;
[maskI, maskJ] = find(origMask);
idx = ones(m,n)*-1;
idx(sub2ind([m n],maskI,maskJ))=1:numel(maskI);

A = sparse(m*n, m*n);
for i=1:m
    for j=1:n
        if origMask(i,j)
            for i1=-1:1
                for j1=-1:1
                 if i+i1>0 && i+i1<=m  && j+j1>0 && j+j1<=n  && origMask(i+i1,j+j1)
                      A(idx(i,j),idx(i+i1,j+j1))=sqrt(i1^2+j1^2); 
                      A(idx(i+i1,j+j1),idx(i,j))=sqrt(i1^2+j1^2); 
                 end
                end
            end
        end
    end
end

myGraph = graph(A);

%%
imgM = m;
imgN = n;
%t = 0:0.01:2;


%elephant ear
%pointA = [25, 31];
%pointB = [78,143];
%pointA = [4,56];
%pointB = [28,57];
%pointB = [35,92];

%fox
%pointA = [80, 70];
%pointB = [127,144];

%stool_1_microcut
%pointA = [122, 28];
%pointA = [35 85];
%pointB = [80, 46];

%stool_1_microcut2
%pointA = [18,21];
%pointA = [31,40];
%pointB = [96, 188];

myPath = shortestpath(myGraph,idx(pointA(1),pointA(2)),idx(pointB(1),pointB(2)));

pathPts = [];
for ii=1:numel(myPath)
pathPts = [pathPts; maskI(myPath(ii)), maskJ(myPath(ii))];
end

%newPathPts = pathPts(1,:);
%for ii=2:numel(myPath)
%    newPathPts = [newPathPts; (newPathPts(end,:)+pathPts(ii,:))*0.5; pathPts(ii,:)];
%end

%pathPts = newPathPts;

x = pathPts;
%x = A+t'*(B-A)/2;
testAlignCurveToFrameField
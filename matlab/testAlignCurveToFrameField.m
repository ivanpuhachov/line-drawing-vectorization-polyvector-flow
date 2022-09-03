%clear
%close all

wReg = 0.1; %regularization with regular mean curvature flow
lengthReg = 1e-8;
wCenter = 0;
%initial curve
%x(:,1) = 0:0.01:2;
%x(:,2) = sin(x(:,1));

%two endpoints
%A = x(1,:);
%B = x(end,:);
% A = [21,31];
% B = [26,53];
% 
%x = A+x(:,1).*(B-A)/2;
hold on
plot(x(:,2),-x(:,1));
%plot(x(:,1),x(:,2));
tmpX{1} = x;
[lengths,tangents,~,normals,~,allTangents,next,curveSize] = computeTangentsEtc(tmpX,true);
lengths = lengths{1}; 
%tangents = [tangents{1}; tangents{1}(end,:)];
n = curveSize+1;
%lineIntegral = @(f)trapz([0;cumsum(lengths{1})],f);

%multiplies nx2x2 by nx2, giving nx2
matMult = @(A,x) sum(A.*reshape(x,[size(A,1),1,2]),3);

%laplacian
%E = [(1:n-1)' (2:n)'];
%A = adjacency_matrix(E);
%L = A-diag(sum(A,2));
%b = find(sum(A,1)==1);
%L(b,:) = 0;

%grad operator
%D = sparse(size(A,1), size(A,2));
%D(sub2ind(size(D),(1:n-1)',E(:,1)))=-1;
%D(sub2ind(size(D),(1:n-1)',E(:,2)))=1;
%D(b,:)=0;

v = 1;
%figure
%axis equal
%axis([0    2.3016   -0.0767    1.7386]);
hold on
%plot(x(:,2),-x(:,1));
N = 100;
alpha = 1e-2;
xInit = x;

%plot(x(:,1),sin(x(:,1)),'LineWidth',5);
xHistory = [];
for i=1:N
    if mod(i,100)==0
        fprintf('Iter %d\n', i);
    end
    %if mod(i,100)==0
        %plot(x(:,1),x(:,2),'Color',[double(i)/N 0 1.0-double(i)/N]);
    xHistory{i} = x;
       %plot(x(:,2),-x(:,1),'Color',[double(i)/N 0 1.0-double(i)/N]);
        %fprintf('It %d/%d: Energy: %f, step size: %f\n',i,N,energy,alpha);
    %end
    
    prev = circshift(x,1,1);
    midpoints = 0.5*(x+prev);
    midpoints(1,:) = x(1,:);
    tau = x-prev;
    tau(1,:)=tau(2,:);
    
    X = x(:,1)-1i*x(:,2);
    q = sqrt(sum(tau.^2,2)+lengthReg); q(1) = q(2);
    tau_cmplx = tau(:,1)-1i*tau(:,2);
    
    qNext = circshift(q,-1);  qNext(end) = inf;
    normals = [-tau(:,2) tau(:,1)]./q;
    
    [g,dg_dn,dg_dx] = anisotropicFrameFieldFun(normals,midpoints,X_cmplx,imgM,imgN);
    
    %[gCenter,dg_dnCenter,dg_dxCenter] = anisotropicImgFun(normals,midpoints,bwImg/255.0,imgM,imgN);
    %g = ones(size(g));
    %dg_dn = zeros(size(dg_dn)); %this is only to test mean curvature flow
    %dg_dx = zeros(size(dg_dx)); %this is only to test mean curvature flow
    %[g,dg_dn,dg_dx] = isotropicTestFun(normals,midpoints);
    %[g,dg_dn,dg_dx] = anisotropicTestFun5(normals,midpoints);
    %g = g+wCenter*gCenter;
    %dg_dn = dg_dn+wCenter*dg_dnCenter;
    %dg_dxCenter = dg_dx+wCenter*dg_dxCenter;
    dg_dx(end,:) = 0;
    %quiver(x(:,1),x(:,2),-dg_dx(:,1),-dg_dx(:,2),0);
    
    gamma_prime = dg_dn(:,1)-1i*dg_dn(:,2);
  
    G = g+1i*real(conj(tau_cmplx).*gamma_prime./q)+wReg;
    F = dot(dg_dx,normals,2);
    F_next = circshift(F,-1); F_next(end)=inf;
    
    xOrig = x;  
    
    GqNext = circshift(G./q,-1); GqNext(end)=inf;
    mainDiag = (1/(2*alpha))*(q+qNext)+G./q + GqNext - 0.5*1i*F + 0.5*1i*F_next;
    upperDiag = [0; -GqNext-0.5*1i*F_next]; upperDiag(end)=[];
    lowerDiag = -G./q+0.5*1i*F; lowerDiag(1) = []; lowerDiag = [lowerDiag; 0];
    
    rhs = (1/(2*alpha))*X.*(q+qNext);
    
    A = spdiags([lowerDiag, mainDiag, upperDiag], -1:1, n, n);
    
    A_cut = A(2:end-1,2:end-1);

    rhs_cut = rhs(2:end-1);
    rhs_cut(1) = rhs_cut(1)-A(2,1)*X(1);
    rhs_cut(end) = rhs_cut(end)-A(end-1,end)*X(end);
    
    X_new = A_cut\rhs_cut;
    X_new = [X(1); X_new; X(end)];
    if (any(isnan(X_new)))
        break;
    end
    %x = xnnn;
    x(:,1) = real(X_new);
    x(:,2) = -imag(X_new);
end

for i=1:numel(xHistory)
    %if mod(i,10)==0
       plot(xHistory{i}(:,2),-xHistory{i}(:,1),'Color',[double(i)/numel(xHistory) 0 1.0-double(i)/numel(xHistory)]);
    %end
    
    if i==numel(xHistory)
        plot(xHistory{i}(:,2),-xHistory{i}(:,1),'Color',[double(i)/numel(xHistory) 0 1.0-double(i)/numel(xHistory)],'LineWidth',5);
    end
end
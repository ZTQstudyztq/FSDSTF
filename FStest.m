function [result] = FStest(X,Y,lambda1,lambda2,lambda3,p,L,alpha,rho1,rho2,rho3,rho4)
epson = 0.1;
iter = 0;
max_iter = 30;

max_rho1 = 1e+10;
max_rho2 = 1e+10;
max_rho3 = 1e+10;
max_rho4 = 1e+10;


gt = Y;                          % the grount truth for testing
cls_num = length(unique(gt));    % the number of clusters
nV = length(X);                  % the number of views

%% variables
XX = NormalizeData(X{1});
for v = 2:nV
    X1 = NormalizeData(X{v});
    XX = cat(1,XX,X1);
end
M = size(XX,1);
N = size(XX,2);
sX = [N, N, 2];

Z = zeros(N,N);
O = Z;
Q3 = O;
C = Z;
G = cat(3,Z,C);
W = G; % W1 is Z, W2 is C

P = zeros(L,M);
J = P+eps;
Q1 = P;
D = zeros(size(P*XX));
Q2 = D;



%% =====================Optimization=====================

while iter < max_iter
    iter = iter + 1;
%     disp(iter);
    E_Z = Z - Z*C;

    %% Update Z
    Z = (2*D'*D + (rho2+rho4)*eye(size(D,2)))^(-1) * (2*D'*P*XX - W(:,:,1) + rho2*G(:,:,1) + rho4*O - Q3);
    
    %% Update O
    O = (rho4*Z + Q3)/2/lambda1 * (eye(size(C,2))*(rho4+2*lambda1)/2/lambda1 - C - C' - C*C')^(-1);

    %% Update C
    C = (2*lambda1*O'*O + rho2*eye(size(O,2)))^(-1)*(2*lambda1*O'*O + rho2*G(:,:,2) - W(:,:,2));
    
    %% Update P
    F1 = 2*D*Z*XX' + rho1*J - Q1 + rho3*D*XX' - Q2*XX';
    [U1,~,V1] = svd(F1, "econ");
    P = U1*V1';
    
    %% Update J
    F2 = P + Q1/rho1;
    JN = size(J,1);
    deta = zeros(JN);
    for i = 1:JN
        deta(i,i) = norm(J(i,:))^(p-2) * p/2;
    end
    J = (rho1*eye(size(deta))/2 + lambda3*deta)^(-1) * rho1*F2/2;
    
    %% Update D
    D = (2*P*XX*Z' + rho3*P*XX + Q2) * (2*Z*Z' + rho3*eye(size(Z,1)))^(-1);
    
    %% Update G
    H = cat(3,Z,C);
    g = zeros(N*N*2,1);    
    h = H(:);
    w = W(:);
    F3 = h+1/rho2*w;
    [g, ~] = wshrinkObj(F3, 2*lambda2/rho2, sX, 0, 3);
    G = reshape(g, sX);

    %% Update multipliers
    Q1 = Q1 + rho1*(P - J);
    Q2 = Q2 + rho3*(P*XX - D);
    Q3 = Q3 + rho4*(Z - O);
    W = W + rho2*(H - G);

    rho1 = min(alpha*rho1, max_rho1);
    rho2 = min(alpha*rho2, max_rho2);
    rho3 = min(alpha*rho3, max_rho3);
    rho4 = min(alpha*rho4, max_rho4);

     %% Check
    if (iter>15 && norm(Z - Z*C -E_Z, "inf")/norm(E_Z, "inf")<epson)
        break;
    end
    
end

CC = (abs(C) + abs(C'))/2;

% The required column is the number of samples
groups = SpectralClustering(CC, cls_num);
result = Clustering8Measure(Y, groups);
end
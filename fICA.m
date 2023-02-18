function [s,W] = fICA(y, N0, M) 
maxit = 200; tol   = 1e-8;
N = size(y,1);
leng = size(y,2);

% # Step 1 中心化
y = y - mean(y,2);
 
% # Step 2 白色化
Omega = y*y'/leng;% - N0 * eye(N);
[U,D,V] = svd(Omega);
V = diag(1./sqrt(diag(D)))*U';
z = V(1:M,1:N)*y;
 
% # Step 3 初期化
W_init = randn(M)+1i*randn(M);
 
% # Step 4 直交化
W = ica_parallel(z,M,maxit,tol,W_init);
     
% # Step 5 出力
s = W'*z;
W = W'*V(1:M,1:N);

function W_res = ica_parallel(z,M,maxit,tol,W_init)
    W_res = W_init;
    W_res = chol(inv(W_res*W_res')) * W_res;
    
    lim = 1000.0*ones(maxit+1,1);
    it = 1;
    while ((lim(it) > tol) && (it < maxit))
        W = zeros(size(W_res));
        for i = 1:M
            wz  = W_res(:,i)'*z;
            wz2 = abs(wz).^2;
            gwz = conj(wz).*wz2;
            gwz = gwz(ones(M,1),:);
            zgwz = z .* gwz;
            W(:,i) = mean(zgwz,2) - mean(2*wz2)*W_res(:,i);
        end
        W = chol(inv(W*W'))*W; %対称的直交化
         
        it = it + 1;
        lim(it) = 0;
        for i = 1:M
            lim(it) = lim(it) + abs( abs(sum(W(i,:) .* W_res(i,:))) - 1.0);
        end
        lim(it) = lim(it)/M;
        W_res = W;
    end
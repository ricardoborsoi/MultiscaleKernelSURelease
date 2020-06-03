function [err_norm_Xic] = unmix_coarse_scale_tmp(mu_0,  KM,MM,M,L,P,M1,numSuperpixels,Y_coarse_mtx,const_C0)


C1 = mu_0;



% Solve KHype for each superpixel ------------------------------
% Unmix each superpixel individually
K = 1*MM + KM +1/C1*eye(L);
H = [ K,    M,         -M1;
      M',   eye(P),    -ones(P,1);     
     -M1', -ones(P,1)', P];
H = (H+H')/2;
H = H+0.0000*eye(L+P+1);


beta_C = zeros(L, numSuperpixels);
a_est = zeros(P,numSuperpixels);
for n=1:numSuperpixels
	
    y = Y_coarse_mtx(n,:)';
	
    % max dual
 
    f = -[y;zeros(P,1);-1];
    A = -[zeros(P,L),eye(P),zeros(P,1)];
    b = zeros(P,1);
    
    z = qpas(H,f,A,b);%,Aeq,beq); 
    beta   = z(1:L);
    gam    = z(L+1:L+P);
    lambda = z(end);
    h = M'*beta+gam-lambda;
    
    a_est(:,n) = h;
    beta_C(:,n) = beta;
end




K = numSuperpixels;
err_norm_Xic = (1/mu_0^2)*norm(beta_C,'fro')^2 - K*const_C0;


% normalize to be around 1
err_norm_Xic = err_norm_Xic / K*const_C0;







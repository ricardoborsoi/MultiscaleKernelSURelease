function [err_norm_Xi, err_norm_dA_Xipsi] = unmix_orig_scale_tmp(vmu,  KM,MM,M,Mpinv,L,P,M1,N,Y,ad,psi_cM,const_C1,const_CY,const_CE)
% [err_norm_Xi, err_norm_dA_Xipsi] = unmix_orig_scale_tmp(mu_1,mu_2,  KM,MM,M,Mpinv,L,P,M1,N,Y,ad,psi_cM,const_C1,const_CY,const_CE)
% [a_est,beta,mu3] = unmix_orig_scale_tmp(mu_1,mu_2,  KM,MM,M,Mpinv,L,P,M1,N,Y,ad,psi_cM)

mu_1 = vmu(1);
mu_2 = vmu(2);


beta = zeros(L,N);
mu3  = zeros(P,N);
a_est = zeros(P,N);

% unmix back in the original scale ---------------------
K1 = KM + (1/(0+mu_2))*MM + (1/mu_1)*eye(L); 
K2 = (1/mu_2)*eye(P) + Mpinv * KM * Mpinv';
% K2 = zeros(P);
H = [K1,                -KM*Mpinv',     (1/(0+mu_2))*M,           -(1/(0+mu_2))*M1;
    -Mpinv*KM,           K2,            zeros(P),                  zeros(P,1);
     (1/(0+mu_2))*M',    zeros(P),      (1/(0+mu_2))*eye(P),      -(1/(0+mu_2))*ones(P,1);     
    -(1/(0+mu_2))*M1',   zeros(1,P),   -(1/(0+mu_2))*ones(P,1)',   (1/(0+mu_2))*P];
H = (H+H')/2;
H = H+0.0000*eye(L+2*P+1);

A = -[zeros(P,L),zeros(P),eye(P),zeros(P,1)];
b = zeros(P,1);

for n=1:N

    y = Y(:,n);

    % max dual

    f = -[-(mu_2/(0+mu_2))*M*ad(:,n)+y;   ...
              -Mpinv*psi_cM(:,n);   ... 
          -(mu_2/(0+mu_2))*ad(:,n);   ...
               (mu_2/(0+mu_2))*ad(:,n)'*ones(P,1)-1];

    z = qpas(H,f,A,b);%,Aeq,beq); 
    % z = quadprog(H,f,A,b);
    beta(:,n) = z(1:L);
    mu3(:,n)  = z(L+1:L+P);
    gam       = z(L+P+1:L+2*P);
    lambda    = z(end);

    h = (1/(0+mu_2)) * (mu_2 * ad(:,n) + M'*beta(:,n) + gam - lambda);
    a_est(:,n) = h;
end



err_norm_Xi = (1/mu_1^2)*norm(beta,'fro')^2 - N*const_C1;
err_norm_dA_Xipsi = norm(a_est-ad,'fro')^2+norm((1/mu_2)*mu3,'fro')^2 - (N*const_CY-N*const_CE);



% normalize to be around 1
err_norm_Xi = err_norm_Xi / N*const_C1;
err_norm_dA_Xipsi = err_norm_dA_Xipsi / (N*const_CY-N*const_CE);




    
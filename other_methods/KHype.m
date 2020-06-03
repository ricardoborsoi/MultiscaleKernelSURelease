function [a_est,beta_khype,rmse_r_KHYPE] = KHype(r,M,C) 
% r -- image
% M -- endmembers
% C -- regularization constant

P = size(M,2);
N = size(r,2);
L = size(r,1);


% Gaussian kernel bandwidth 
% par = 2;
% Regualrization parameter 
% C = 100; % \mu in the paper = 1/C

% Gaussian kernel calculation
% Q    = eye(P)/par^2;
% MQM  = M*Q*M';
% dMQM = diag(MQM);
% KM   = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));
% For using the polynomial proposed polynomial kernel, remove the comment symbol
KM = (1+1/P^2*(M-0.5)*(M-0.5)').^2;

M1    = M*ones(P,1);
a_est = zeros(P,N);
MM    = M*M';

K = 1*MM+KM;
K = [K+1/C*eye(L)];
H = [ K,    M,         -M1;
      M',   eye(P),    -ones(P,1);     
     -M1', -ones(P,1)', P];
H = (H+H')/2;
H = H + 0.0000*eye(L+P+1);


tic
beta_khype = zeros(L,N);
% Algorithm K-Hype 
for n = 1 : N
    y = r(:,n);
    % max dual

    f = -[y;zeros(P,1);-1];
    A = -[zeros(P,L),eye(P),zeros(P,1)];
    b = zeros(P,1);
    
    z = qpas(H,f,A,b);%,Aeq,beq); 
    beta   = z(1:L);
    gam    = z(L+1:L+P);
    lambda = z(end);
    h      = M'*beta+gam-lambda;
    a_est(:,n) = h;
    
    beta_khype(:,n) = beta;
    
end


rmse_r_KHYPE = sqrt(norm(r-M*a_est - KM*beta_khype,'fro')^2);



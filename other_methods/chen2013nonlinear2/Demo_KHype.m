% =========================================================
% K-Hype algorithm
% ---------------------------------------------------------
% Related paper: 
% 
% J. CHEN, C. RICHARD and P. HONEINE, "Nonlinear unmixing of hyperspectral
% data based on a linear-mixture/nonlinear-fluctuation model". IEEE Transactions on signal processing, 2012. 
% 
% ---------------------------------------------------------
% Utilization recommendation:
% 1. Parameters (Kernel bandwidth/regularization parameter) should be tuned according to scenes.
% 2. Here we use the function "qpas" for the optimization instead of the general optimization 
%    function of matlab: "quadprog", as the later is rather slow. One can
%    replace all "qpas" by "quadprog", or other optimization functions in case
%    there's a problem of running.
%    The full toolbox QuadprogC is also attached for respecting its completeness. 
%
%
% Version: 13-sept-2012
%
% Jie CHEN & Cedric Richard
% chen@unice.fr, cedric.richard@unice.fr
% Laboratoire Lagrange, Université de Nice Sophia-antipolis
% Nice, France


clear



% ============  Parameters to tune ======================== 
% Gaussian kernel bandwidth : 
par = 2;
% Regualrization parameter : 
% \mu in the paper = 1/C
C = 100;


% ============== Scene Parameters =========================
%  (see the function generate_image for detailes)
% endmember database (only 3 availabe in this demo due to endmembers.mat)
Material_database = 3;
% number of endmember (maximum to 8 in this demo due to the endmembers.mat)
R = 5; 
% mixture model (1: linear, 2: bilinear, 3:postnonlinear)
Mixture_model=3;
% number of pixels to generate
N = 10;
% snr
SNR = 30;

%  ============== Image generation  =======================
[r,M,a,noise_std]=generate_image(Material_database, R, Mixture_model, N, SNR);
[L, R] = size(M);


% ============= Gaussian kernel calculation =================
Q=eye(R)/par^2;
MQM=M*Q*M';
dMQM = diag(MQM);
KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));
% For using the polynomial proposed polynomial kernel, remove the
% comment symbol:
% KM = (1+1/R^2*(M-0.5)*(M-0.5)').^2;


M1 = M*ones(R,1);
a_est = zeros(R,N);
MM = M*M';

tic
% =================== Algorithm K-Hype =======================
for n = 1 : N
    y = r(:,n);
   % max dual
    K = 1*MM+KM;
    K = [K+1/C*eye(L)];
    H = [K,    M,    -M1;
         M',  eye(R), -ones(R,1);     
         -M1', -ones(R,1)',R];
    H=(H+H')/2;
    H=H+0.0000*eye(L+R+1);
 
    f = -[y;zeros(R,1);-1];
    A = -[zeros(R,L),eye(R),zeros(R,1)];
    b = zeros(R,1);
        
    
    z = qpas(H,f,A,b);%,Aeq,beq); 
    beta=z(1:L);
    gam=z(L+1:L+R);
    lambda=z(end);
    h = M'*beta+gam-lambda;
    a_est(:,n)=h;

end
toc

% =============== RMSE calculation ========================= 
[RMSE, std] = ErrComput(a,a_est)



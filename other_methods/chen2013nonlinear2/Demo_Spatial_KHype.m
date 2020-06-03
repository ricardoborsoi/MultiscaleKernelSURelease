% Adaptive spatial regularization for hyperspectral image unmixing
% - Nonlinear unmxing 
% - Grid image (DC1)
% - dynamic lambda_breg


clear;



%%%%%%%%paramters to tune%%%%%%%%
% Functional regularization
C = 200;
% Reg. parameter for spatial corr.
lambda_sp =0.5;
%%%%%%%%%%%%%%%%%   

%% Scene Parameters

% Generate the grid abundances
grid_image;
%image_a =image_a(1:20,46:65,:);
N = size(image_a,1);
N1=N;
N=N1*N1;

IHIV=IpHHT(N1);

disp('--inv--')

% Transform the square image to a matrix
a =reshape(image_a,N,1,R);
a = squeeze(a)';

% Materaial database
Material_database = 0;
% Mixture model
Mixture_model = 2;
% Number of endmember in this image is 5;
R = 5;
% SNR
SNR = 20;

%% Generate the image
[r,M,a,noise_std]=generate_image(Material_database, R, Mixture_model, N, SNR, a);
[L, R] = size(M);

%% Algorithm init.
% The size the image
w =N1; h =N1;
% Abundances init.
a_est=zeros(R,N);





lambda_breg = 10;
% Nonlinear kernel matrix

par = 2;
Q=eye(R)/par^2;
MQM=M*Q*M';
dMQM = diag(MQM);
KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));%-0.0352*M*M';
KM = (1+1/R^2*(M-0.5)*(M-0.5)').^2;
MM =M*M';
K = MM+KM;


                  
                         
                         
M1 = M*ones(R,1);

A = -[zeros(R,L),eye(R),zeros(R,1)];
b = zeros(R,1);
%% Algorithm run

V = zeros(R,N);
D1 = V;
D2 = zeros(R,4*N);
U=D2;

lambda_breg_c=10;
tic

for iter = 1 : 10
    iter
    if iter == 1 
      lambda_breg = 0; %iter;
    else
        lambda_breg=lambda_breg_c;
    end
%    lambda_breg = iter-1;
    
    
rho = 1/(1+lambda_breg);
Kn = KM + rho*MM;
Kn = [Kn+1/C*eye(L)];
Hn = rho* [Kn/rho,    M,    -M1;
          M',  eye(R), -ones(R,1);
           -M1', -ones(R,1)',R];         
Hn=(Hn+Hn')/2;
Hn=Hn+0.00001*eye(L+R+1);
  
    
    
    An = V + D1;
 %   An = convH(A);
    a_est0=a_est;
    
tic
for n = 1 : N
    if mod(n,100)==0, n, end
    
    rn = r(:,n);

    an = An(:,n);
    f = -[ rn-rho*lambda_breg*M*an;
        -rho*lambda_breg*an;
        rho*lambda_breg*an'*ones(R,1)-1];
    z = qpas(Hn,f,A,b);%,Aeq,beq);
    beta=z(1:L);
    gam=z(L+1:L+R);
    lambda=z(end);
    at =  rho*(M'*beta+gam-lambda+lambda_breg*an);
    a_est(:,n) = at;
end
toc

t2=U-D2;
%V = (a_est-D1+ConvHT(t2))*IHIV;
V_=V;
V=(a_est-D1+ConvHT(t2))/IHIV;
%V=max(V,0);
%V=V./(ones(R,1)*sum(V));
U0=ConvH(V);
U_=U;
U = sign(U0+D2).*max(abs(U0+D2)-lambda_sp/lambda_breg,0);

D1 = D1 + (V-a_est);
D2 = D2+(U0-U);


if norm(a_est-V,'fro')/R/N <1e-5 & norm(U-U0,'fro')/R/N/4 <1e-5
%if norm(U-ConvH(a_est))/R/N<1e-5
    break;
end



s1 = norm(lambda_breg_c*(V-V_),'fro');
s2 = norm(lambda_breg_c*ConvHT(U-U_),'fro');
s=s1+s2;
rx=norm(V-a_est,'fro')+norm(U-U0,'fro');


if rx > 10*s
      lambda_breg_c = 2*lambda_breg_c;
elseif s>10*rx
      lambda_breg_c = lambda_breg_c/2;
else
      lambda_breg_c = lambda_breg_c;
end




end

toc

[RMSE, std] = ErrComput(a, a_est)

b=reshape(a_est',w,h,[]);

figure,
for i = 1 : R
   subplot(1,R,i), imshow(b(:,:,i));
end
colormap jet





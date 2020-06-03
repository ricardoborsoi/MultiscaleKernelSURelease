function [a_est_spatialKhype,rmse_r_TV,mse_TV_pixelwise,beta_TV]=KhypeTV_wrapper2(r,N,P,nr,nc,M,L,lambda_sp,C)





% =========================================================================
% TV

% if false

% paramters to tune
% Functional regularization
% % % % % % % % % % % % % C = 200; 20; % 200;
% C = 200; 20; % 200;
% Reg. parameter for spatial corr.
% % % % % % % % % % % % lambda_sp = 0.5;
par = 2;

tic

% ASSUMING nr=nc ---------------------
% IHIV = IpHHT(nr);
IHIV = IpHHT2(N);
% IHIV = IpHHT3(nc,nr);%%%%%%%%%%%%%%%%%%%%%%

% The size the image
w = nc; h = nr;
% Abundances init.
a_est = zeros(P,N);

lambda_breg = 10;
% Nonlinear kernel matrix
% Q = eye(P)/par^2;
% MQM  = M*Q*M';
% dMQM = diag(MQM);
% KM = exp(-0.5*(dMQM*ones(L,1)'+ones(L,1)*dMQM'-2*M*Q*M'));%-0.0352*M*M';
% Polynomial kernel
KM = (1+1/P^2*(M-0.5)*(M-0.5)').^2;

MM = M*M';
K  = MM+KM;
M1 = M*ones(P,1);
A  = -[zeros(P,L),eye(P),zeros(P,1)];
b  = zeros(P,1);



% algorithm
V  = zeros(P,N);
D1 = V;
D2 = zeros(P,4*N);
U  = D2;

lambda_breg_c = 10;

beta_TV = zeros(L,N);
for iter = 1 : 10
    disp(iter)
    if iter == 1 
      lambda_breg = 0; %iter;
    else
        lambda_breg = lambda_breg_c;
    end
%    lambda_breg = iter-1;
    
    rho = 1/(1+lambda_breg);
    Kn = KM + rho*MM;
    Kn = [Kn+1/C*eye(L)];
    Hn = rho* [Kn/rho,    M,    -M1;
              M',  eye(P), -ones(P,1);
               -M1', -ones(P,1)',P];         
    Hn=(Hn+Hn')/2;
    Hn=Hn+0.00001*eye(L+P+1);



    An = V + D1;
 %   An = convH(A);
    a_est0=a_est;

%     tic
    for n = 1 : N
%         if mod(n,100)==0, disp(n), end

        rn = r(:,n);

        an = An(:,n);
        f = -[ rn-rho*lambda_breg*M*an;
            -rho*lambda_breg*an;
            rho*lambda_breg*an'*ones(P,1)-1];
        z = qpas(Hn,f,A,b);%,Aeq,beq);
        beta   = z(1:L);
        gam    = z(L+1:L+P);
        lambda = z(end);
        at =  rho*(M'*beta+gam-lambda+lambda_breg*an);
        a_est(:,n) = at;
        
        beta_TV(:,n) = beta;
    end
%     toc

    t2=U-D2;
    %V = (a_est-D1+ConvHT(t2))*IHIV;
    V_=V;
    V=(a_est-D1+ConvHT(t2))/IHIV;
    %V=max(V,0);
    %V=V./(ones(R,1)*sum(V));
    U0=ConvH(V);
    % V=(a_est-D1+ConvHT3(t2,nc,nr))/IHIV; %%%%%%%%%%%%%%%%%
    % U0=ConvH3(V,nc,nr); %%%%%%%%%%%%%%%%%
    
    U_=U;
    U = sign(U0+D2).*max(abs(U0+D2)-lambda_sp/lambda_breg,0);

    D1 = D1 + (V-a_est);
    D2 = D2+(U0-U);


    if norm(a_est-V,'fro')/P/N <1e-5 && norm(U-U0,'fro')/P/N/4 <1e-5
    %if norm(U-ConvH(a_est))/R/N<1e-5
        break;
    end

    s1 = norm(lambda_breg_c*(V-V_),'fro');
    s2 = norm(lambda_breg_c*ConvHT(U-U_),'fro');
    % s2 = norm(lambda_breg_c*ConvHT3(U-U_,nc,nr),'fro');%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = s1+s2;
    rx = norm(V-a_est,'fro')+norm(U-U0,'fro');
    
    
    if rx > 10*s
          lambda_breg_c = 2*lambda_breg_c;
    elseif s>10*rx
          lambda_breg_c = lambda_breg_c/2;
    else
          lambda_breg_c = lambda_breg_c;
    end

end
time_spatialKhype = toc;

a_est_spatialKhype = a_est;

% end


% time_spatialKhype = 10;
% 
% a_est_spatialKhype = a_est;

rmse_r_TV = norm(r-M*a_est_spatialKhype - KM*beta_TV,'fro');
mse_TV_pixelwise = sum((r-M*a_est_spatialKhype - KM*beta_TV).^2,1);




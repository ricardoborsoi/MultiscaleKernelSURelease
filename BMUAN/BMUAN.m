function [a_est_mscale,rmse_r_SPPX]=BMUAN(Y,M,nr,nc,sigma2_epsi)


% if the constants are not estimated from the data, please substitute "noise" by the true noise present in the image
noise = eps * ones(size(r)); % dummy placeholder (not used since the constants are estimated from the data)


L = size(Y,1);
P = size(M,2);
N = nr*nc;




run('vl_setup')

% Gaussian kernel bandwidth 
% par = 2;
% Regualrization parameter 
slic_reg   = 0.005;

%{
% Norm of modeling errors 
% sigma2_epsi = (1e-4) * norm(Y,'fro')^2/N; 
% sigma2_epsi = (1e-5) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
% sigma2_epsi = (1e-6) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
sigma2_epsi = (1e-8) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
% sigma2_epsi = (0e-8) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
%}

% estimate constants or use known (oracle) signals
flag_estimate_constants = true; 


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



% Y is the HSI ordered as a matrix (columnwise)
Y2 = reshape(Y', nr, nc, L);   
Y2a = Y2;


% ======================================================================
% Estimate superpixel size and noise covariance matrix -----------------
% slic_size = 5; 
[slic_size] = homogeneity_curvature_sppx_size(Y2,slic_reg,0.1); 
fprintf('\n Selected SLIC size: %d \n', slic_size)

% Estimate noise covariance matrix ---------------
[~,Sigma_e]=estimate_noise_covar(Y,P);

% try to "clean" the covariance matrix (assume a diagonal structure)
Sigma_et = diag(diag(Sigma_e));
% ASSUME CORRELATION BETWEEN NOISE AT ADJACENT BANDS (warning: off-diagonal estimates are very bad here)
% Sigma_et = diag(diag(Sigma_e)) + diag(diag(Sigma_e,1),1) + diag(diag(Sigma_e,-1),-1);
Sigma_e  = Sigma_et * norm(Sigma_e,'fro')/norm(Sigma_et,'fro'); % correct the lost noise power
% Sigma_e  = Sigma_et; % uncomment to override and use the original estimate



% ======================================================================
% reorder and rescale data into 2-D array -------------------------
[numRows,numCols,numSpectra] = size(Y2);
scfact = mean(reshape(sqrt(sum(Y2.^2,3)), numRows*numCols, 1));
Y2 = Y2./scfact;
imgVec = reshape(Y2, [numRows*numCols numSpectra]);

% compute superpixels
disp('Computing SLIC Superpixels...');
spSegs = vl_slic(single(Y2), slic_size, slic_reg);
numSuperpixels = double(max(spSegs(:)))+1; 


% Unmix the superpixels -------------------------
Y3 = zeros(size(Y2));
avg_superpx = zeros(1, numSuperpixels+1, L);
avg_sppxld_noise = zeros(1, numSuperpixels+1, L); % Compute the coarse representation of the noise
noiseim = reshape(noise', nr, nc, L);  
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    for j=1:length(rowi)
        % Averages all pixels inside each superpixel
        if j == 1
            avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
            avg_sppxld_noise(1,i+1,:) = (1/length(rowi)) * noiseim(rowi(j),coli(j),:);
        else
            avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
            avg_sppxld_noise(1,i+1,:) = avg_sppxld_noise(1,i+1,:) + (1/length(rowi)) * noiseim(rowi(j),coli(j),:);
        end
    end
    
    % This is optional (for visualization)
    for j=1:length(rowi)
        Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
    end
end

% Remove redundant dimensions
Y_coarse_mtx = squeeze(avg_superpx);



% ======================================================================
% Solve KHype for each superpixel ------------------------------------
% Perform bisection on the coarse scale
if flag_estimate_constants == false
    const_C0 = (1/(slic_size^2)) * (1/N) * norm(noise,'fro')^2 ...
               + sigma2_epsi; % true value cor C0
else
    const_C0 = (1/(slic_size^2)) * trace(Sigma_e) + sigma2_epsi;
end



ffunC = @(vmu)unmix_coarse_scale_tmp(vmu,  KM,MM,M,L,P,M1,numSuperpixels,Y_coarse_mtx,const_C0);


% for bisection method
mu_0_a = 0.001;
mu_0_b = 1000;

% evaluate at the first extremes
[err_norm_Xic_a] = ffunC(mu_0_a);
[err_norm_Xic_b] = ffunC(mu_0_b);


% check if a root exists in the interval
if err_norm_Xic_a*err_norm_Xic_b > 0
    disp('Warning: Coarse scale bisection: There is no root in the selected interval!!!')
    if err_norm_Xic_b <= err_norm_Xic_a
        mu_0 = mu_0_b;
    else
        mu_0 = mu_0_a;
    end
    
else % if a root exists, perform Bisection...
    fprintf('\n\n Bissecting ... \n')
    for i=1:10
        % compute the value at the middle and check corners
        % mu_0_c = (mu_0_a+mu_0_b)/2; % partition in linear space
        mu_0_c = 10^((log10(mu_0_a)+log10(mu_0_b))/2); % partition in log-space

        % lets check each half-space
        err_norm_Xic_c = ffunC(mu_0_c);

        if err_norm_Xic_c*err_norm_Xic_b < 0
            mu_0_a = mu_0_c;
            err_norm_Xic_a = err_norm_Xic_c;
            fprintf('\n (2)')
        else
            mu_0_b = mu_0_c;
            err_norm_Xic_b = err_norm_Xic_c;
            fprintf('\n (1)')
        end

        % check for termination criteria
        mu_0 = (mu_0_a+mu_0_b)/2;
        if (mu_0_b-mu_0_a)/2 < 0.2*mu_0
            break;
        end
    end
end



% Solve KHype for each superpixel --------
% Unmix each superpixel individually
K = 1*MM + KM +1/mu_0*eye(L);
H = [ K,    M,         -M1;
      M',   eye(P),    -ones(P,1);     
     -M1', -ones(P,1)', P];
H = (H+H')/2;
H = H+0.0000*eye(L+P+1);

beta_C = zeros(L, numSuperpixels);
for n=1:numSuperpixels
    y = Y_coarse_mtx(n,:)';
	
    % max dual
    f = -[y;zeros(P,1);-1];
    A = -[zeros(P,L),eye(P),zeros(P,1)];
    b = zeros(P,1);
    
    z = qpas(H,f,A,b);%,Aeq,beq); 
%     z = quadprog(H,f,A,b); %,Aeq,beq); 
    beta   = z(1:L);
    gam    = z(L+1:L+P);
    lambda = z(end);
    h = M'*beta+gam-lambda;
    
    a_est(:,n) = h;
    beta_C(:,n) = beta;
end

fprintf('\n norms: %f .... constants %f:', (1/mu_0^2)*norm(beta_C,'fro')^2, numSuperpixels*const_C0)
fprintf('\n')



% ======================================================================
% Re-attribute the abundances for the entire matrix ---------------
temp = zeros(size(Y2,1), size(Y2,2), P);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = a_est(:,i+1);
    end
end
ad = reshape(temp, [size(Y2,1)*size(Y2,2) P])'; % ad is good
a_est = zeros(P,N);


% Find coarse representation image in original domain
temp = zeros(size(Y2,1), size(Y2,2), L);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = Y_coarse_mtx(i+1,:)';
    end
end
Yd = reshape(temp, [size(Y2,1)*size(Y2,2) L])';


% Find coarse representation noise in original domain
temp = zeros(size(Y2,1), size(Y2,2), L);
noise_coarse_mtx = squeeze(avg_sppxld_noise);
for i=0:numSuperpixels
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = noise_coarse_mtx(i+1,:)';
    end
end
noise_d = reshape(temp, [size(Y2,1)*size(Y2,2) L])';


% Estimate nonlinear contribution psi_C(M) at the coarse scale
temp = zeros(size(Y2,1), size(Y2,2), L);
for i=0:numSuperpixels-1
    
    psi_cMi = KM * beta_C(:,i+1);
    [rowi, coli] = find(spSegs==i);
    
    % Attributes unmixing result to all pixels in a voxel
    for j=1:length(rowi)
        temp(rowi(j),coli(j),:) = psi_cMi;
    end
end
psi_cM = reshape(temp, [size(Y2,1)*size(Y2,2) L])';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the constants
if flag_estimate_constants == false
    Mpinv = pinv(M);
    const_C1 = (1/N)*norm(noise,'fro')^2;
    const_CY = (1/N)*norm(Mpinv*(Y-Yd),'fro')^2;
    const_CE = (1/N)*norm(Mpinv*(noise-noise_d),'fro')^2;
    const_C0 = (1/(slic_size^2)) * norm(noise,'fro')^2 + sigma2_epsi;
else
    Mpinv = pinv(M);
    const_C1 = trace(Sigma_e) + sigma2_epsi;
    const_CY = (1/N)*norm(Mpinv*(Y-Yd),'fro')^2;
    const_CE = trace(Mpinv*Sigma_e*Mpinv') * (slic_size^2-1)/(slic_size^2);
    const_C0 = (1/(slic_size^2)) * trace(Sigma_e) + sigma2_epsi;
end


% avoid an unfeasible constraint (i.e. like \|x\|<0) 
if const_CY < const_CE
    const_CE = const_CY - 1e-8;
end



% ======================================================================
% ------------------------------------
% Perform unmixing again in the fine scale

ffun = @(vmu)unmix_orig_scale_tmp(vmu,  KM,MM,M,Mpinv,L,P,M1,N,Y,ad,psi_cM,const_C1,const_CY,const_CE);

% for bisection method
mu_1_a = 0.001;
mu_1_b = 500;
mu_2_a = 0.001;
mu_2_b = 250;

% b _________________
%   |                |
%   |                |
%   |       c        |
%   |                |
% a |________________|
%   a'              b'

% evaluate at the first corners
[err_norm_Xi_aa, err_norm_dA_Xipsi_aa] = ffun([mu_1_a, mu_2_a]);
[err_norm_Xi_ab, err_norm_dA_Xipsi_ab] = ffun([mu_1_a, mu_2_b]);
[err_norm_Xi_ba, err_norm_dA_Xipsi_ba] = ffun([mu_1_b, mu_2_a]);
[err_norm_Xi_bb, err_norm_dA_Xipsi_bb] = ffun([mu_1_b, mu_2_b]);


fprintf('\n\n Bissecting ... \n')
for i=1:10
    % compute the value at the middle and check corners
    % partition in the linear space
%     mu_1_c = (mu_1_a+mu_1_b)/2;
%     mu_2_c = (mu_2_a+mu_2_b)/2;
    % partition in log-space:
	mu_1_c = 10^((log10(mu_1_a)+log10(mu_1_b))/2);
    mu_2_c = 10^((log10(mu_2_a)+log10(mu_2_b))/2);
    
    
    % lets check each half-space ---------------------
    % first mu_1
    [err_norm_Xi_ca, err_norm_dA_Xipsi_ca] = ffun([mu_1_c, mu_2_a]);
    [err_norm_Xi_cb, err_norm_dA_Xipsi_cb] = ffun([mu_1_c, mu_2_b]);
    
    % update mu_1
    if (abs(sign(err_norm_Xi_aa)+sign(err_norm_Xi_ab)+sign(err_norm_Xi_ca)+sign(err_norm_Xi_cb)) < 4) && ...
        (abs(sign(err_norm_dA_Xipsi_aa)+sign(err_norm_dA_Xipsi_ab)+sign(err_norm_dA_Xipsi_ca)+sign(err_norm_dA_Xipsi_cb)) < 4)
        % smaller rectangle
        mu_1_b = mu_1_c;
        err_norm_Xi_ba = err_norm_Xi_ca;
        err_norm_Xi_bb = err_norm_Xi_cb;
        err_norm_dA_Xipsi_ba = err_norm_dA_Xipsi_ca;
        err_norm_dA_Xipsi_bb = err_norm_dA_Xipsi_cb;
        
        fprintf('\n (1,')
    else % larger rectangle
        mu_1_a = mu_1_c;
        
        err_norm_Xi_aa = err_norm_Xi_ca;
        err_norm_Xi_ab = err_norm_Xi_cb;
        err_norm_dA_Xipsi_aa = err_norm_dA_Xipsi_ca;
        err_norm_dA_Xipsi_ab = err_norm_dA_Xipsi_cb;
        
        fprintf('\n (2,')
    end
    
    
    % now partition in mu_2 --------------------------
    [err_norm_Xi_ac, err_norm_dA_Xipsi_ac] = ffun([mu_1_a, mu_2_c]);
    [err_norm_Xi_bc, err_norm_dA_Xipsi_bc] = ffun([mu_1_b, mu_2_c]);
    
    % update mu_2
    if (abs(sign(err_norm_Xi_aa)+sign(err_norm_Xi_ba)+sign(err_norm_Xi_ac)+sign(err_norm_Xi_bc)) < 4) && ...
        (abs(sign(err_norm_dA_Xipsi_aa)+sign(err_norm_dA_Xipsi_ba)+sign(err_norm_dA_Xipsi_ac)+sign(err_norm_dA_Xipsi_bc)) < 4)
        % smaller rectangle
        mu_2_b = mu_2_c;

        err_norm_Xi_ab = err_norm_Xi_ac;
        err_norm_Xi_bb = err_norm_Xi_bc;
        err_norm_dA_Xipsi_ab = err_norm_dA_Xipsi_ac;
        err_norm_dA_Xipsi_bb = err_norm_dA_Xipsi_bc;
        
        fprintf('1)')
        
    elseif (abs(sign(err_norm_Xi_ac)+sign(err_norm_Xi_bc)+sign(err_norm_Xi_ab)+sign(err_norm_Xi_bb)) < 4) && ...
        (abs(sign(err_norm_dA_Xipsi_ac)+sign(err_norm_dA_Xipsi_bc)+sign(err_norm_dA_Xipsi_ab)+sign(err_norm_dA_Xipsi_bb)) < 4)
        % larger rectangle
        mu_2_a = mu_2_c;
        
        err_norm_Xi_aa = err_norm_Xi_ac;
        err_norm_Xi_ba = err_norm_Xi_bc;
        err_norm_dA_Xipsi_aa = err_norm_dA_Xipsi_ac;
        err_norm_dA_Xipsi_ba = err_norm_dA_Xipsi_bc;
        
        fprintf('2)')
        
    else % if there is not a root in the interval we have a problem -------

        % check if there is a root only for mu_2, otherwise take a good guess
        if (abs(sign(err_norm_dA_Xipsi_ac)+sign(err_norm_dA_Xipsi_bc)+sign(err_norm_dA_Xipsi_ab)+sign(err_norm_dA_Xipsi_bb)) < 4)
            % larger rectangle
            mu_2_a = mu_2_c;
            err_norm_Xi_aa = err_norm_Xi_ac;
            err_norm_Xi_ba = err_norm_Xi_bc;
            err_norm_dA_Xipsi_aa = err_norm_dA_Xipsi_ac;
            err_norm_dA_Xipsi_ba = err_norm_dA_Xipsi_bc;
            fprintf('2)')
        elseif (abs(sign(err_norm_dA_Xipsi_aa)+sign(err_norm_dA_Xipsi_ba)+sign(err_norm_dA_Xipsi_ac)+sign(err_norm_dA_Xipsi_bc)) < 4)
            % smaller rectangle
            mu_2_b = mu_2_c;
            err_norm_Xi_ab = err_norm_Xi_ac;
            err_norm_Xi_bb = err_norm_Xi_bc;
            err_norm_dA_Xipsi_ab = err_norm_dA_Xipsi_ac;
            err_norm_dA_Xipsi_bb = err_norm_dA_Xipsi_bc;
            fprintf('1)')
        else 
            % if there is no root at all try the best guess
            if min(err_norm_dA_Xipsi_aa,err_norm_dA_Xipsi_ba) <= min(err_norm_dA_Xipsi_ab,err_norm_dA_Xipsi_bb)
                % smaller error is for small mu_2 (smaller rectangle)
                mu_2_b = mu_2_c;
                err_norm_Xi_ab = err_norm_Xi_ac;
                err_norm_Xi_bb = err_norm_Xi_bc;
                err_norm_dA_Xipsi_ab = err_norm_dA_Xipsi_ac;
                err_norm_dA_Xipsi_bb = err_norm_dA_Xipsi_bc;
                fprintf('1)')
            else
                % smaller error is for large mu_2 (larger rectangle) 
                mu_2_a = mu_2_c;
                err_norm_Xi_aa = err_norm_Xi_ac;
                err_norm_Xi_ba = err_norm_Xi_bc;
                err_norm_dA_Xipsi_aa = err_norm_dA_Xipsi_ac;
                err_norm_dA_Xipsi_ba = err_norm_dA_Xipsi_bc;
                fprintf('2)')
            end
        end
    end
    
    
    % check for termination criteria
    mu_1 = (mu_1_a+mu_1_b)/2;
    mu_2 = (mu_2_a+mu_2_b)/2;
    if (mu_1_b-mu_1_a)/2 < 0.2*mu_1  &&  (mu_2_b-mu_2_a)/2 < 0.2*mu_2
        break;
    end
%     disp('iter...')
%     mu_1 = (mu_1_a+mu_1_b)/2
%     mu_2 = (mu_2_a+mu_2_b)/2
end


%%

% unmix back in the original scale ---------------------
K1 = KM + (1/(0+mu_2))*MM + (1/mu_1)*eye(L); 
K2 = (1/mu_2)*eye(P) + Mpinv * KM * Mpinv';
H = [K1,                -KM*Mpinv',     (1/(0+mu_2))*M,           -(1/(0+mu_2))*M1;
    -Mpinv*KM,           K2,            zeros(P),                  zeros(P,1);
     (1/(0+mu_2))*M',    zeros(P),      (1/(0+mu_2))*eye(P),      -(1/(0+mu_2))*ones(P,1);     
    -(1/(0+mu_2))*M1',   zeros(1,P),   -(1/(0+mu_2))*ones(P,1)',   (1/(0+mu_2))*P];
H = (H+H')/2;
H = H+0.0000*eye(L+2*P+1);

A = -[zeros(P,L),zeros(P),eye(P),zeros(P,1)];
b = zeros(P,1);

beta = zeros(L,N);
mu3  = zeros(P,N);
for n=1:N
    y = Y(:,n);

    % max dual
    f = -[-(mu_2/(0+mu_2))*M*ad(:,n)+y;   ...
          -Mpinv*psi_cM(:,n);   ... 
          -(mu_2/(0+mu_2))*ad(:,n);   ...
           (mu_2/(0+mu_2))*ad(:,n)'*ones(P,1)-1];

    z = qpas(H,f,A,b); %,Aeq,beq); 
    % z = quadprog(H,f,A,b);
    beta(:,n) = z(1:L);
    mu3(:,n)  = z(L+1:L+P);
    gam       = z(L+P+1:L+2*P);
    lambda    = z(end);

    h = (1/(0+mu_2)) * (mu_2 * ad(:,n) + M'*beta(:,n) + gam - lambda);
    a_est(:,n) = h;
end

fprintf('\n norms: %f .... constants %f:', (1/mu_1^2)*norm(beta,'fro')^2, N*const_C1)
fprintf('\n norms: %f .... constants %f:', norm(a_est-ad,'fro')^2+norm((1/mu_2)*mu3,'fro')^2, N*const_CY-N*const_CE)
fprintf('\n')




a_est_mscale = a_est;


% rmse_r_SPPX = norm(r-M*a_est_mscale-KM*beta,'fro');
rmse_r_SPPX = sqrt((1/mu_1^2)*norm(beta,'fro')^2);





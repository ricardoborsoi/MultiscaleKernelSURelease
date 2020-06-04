% =========================================================================
% 
% This code contains the example with the Urban image in the following paper:
% 
%    A Blind Multiscale Spatial Regularization Framework for Kernel-Based Spectral Unmixing 
%    R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
%    IEEE Transactions on Image Processing, 2020.
% =========================================================================



clear
close all
clc
warning off;


addpath(genpath('other_methods'))
addpath(genpath('DATA'))
addpath(genpath('utils'))
addpath(genpath('BMUAN'))


rng(10, 'twister') 




load('DATA/subimage_Urban_R162.mat')

[nr,nc,L] = size(im);
N = nr*nc;

r = reshape(im,nr*nc,L)';
Yim = reshape(r', nr, nc, L);





% figure;
% clear rgb
% v = [50 25 10 ];
% rgb(:,:,1) = imadjust(rescale(im(:,:,v(1)),0,1));
% rgb(:,:,2) = imadjust(rescale(im(:,:,v(2)),0,1));
% rgb(:,:,3) = imadjust(rescale(im(:,:,v(3)),0,1));
% imshow(1.5*rgb) % display RGB image
% set(gca,'ytick',[],'xtick',[])
% figure, imagesc(10*im(:,:,v))
% figure, imagesc(1.7*im(:,:,[100 60 40]))
% set(gca,'xtick',[]), set(gca,'xticklabel',[]), set(gca,'ytick',[]), set(gca,'yticklabel',[])




% Extract EMs from each class
P = 3;

% [M0, V, U, Y_bar, endm_proj, Y_proj] = find_endm(r,P,'vca');
load('DATA/M0_urban.mat')
M = M0;



% =========================================================================
%  FCLS

tic
A_FCLS = FCLSU(r,M);
time_fcls = toc;

rmse_r_FCLS = norm(r-M*A_FCLS','fro');



%%
% =========================================================================
%  K-Hype

% Regualrization parameter 
C = 200; % \mu in the paper = 1/C

tic
[a_est_khype,beta_khype,rmse_r_KHYPE] = KHype(r,M,C);
time_khype = toc;





%%
% =========================================================================
%  TV (process it by parts to circumvent a bug)

% parameters
lambda_sp = 0.1;
C = 1000;

tic 
r_temp_TV1 = im(1:nr, 1:nr, :);
r_temp_TV2 = im(1:nr, (nr+1):(2*nr), :);
r_temp_TV3 = im(1:nr, end-nr+1:end, :);

r_temp_TV1 = reshape(r_temp_TV1,nr*nr,[])';
r_temp_TV2 = reshape(r_temp_TV2,nr*nr,[])';
r_temp_TV3 = reshape(r_temp_TV3,nr*nr,[])';

[a_est_spatialKhype1,~,mse_TV_pixelwise1,beta_TV1]=KhypeTV_wrapper2(r_temp_TV1,nr*nr,P,nr,nc,M,L,lambda_sp,C);
[a_est_spatialKhype2,~,mse_TV_pixelwise2,beta_TV2]=KhypeTV_wrapper2(r_temp_TV2,nr*nr,P,nr,nc,M,L,lambda_sp,C);
[a_est_spatialKhype3,~,mse_TV_pixelwise3,beta_TV3]=KhypeTV_wrapper2(r_temp_TV3,nr*nr,P,nr,nc,M,L,lambda_sp,C);

a_est_spatialKhype1 = reshape(a_est_spatialKhype1', nr, nr, P);
a_est_spatialKhype2 = reshape(a_est_spatialKhype2', nr, nr, P);
a_est_spatialKhype3 = reshape(a_est_spatialKhype3', nr, nr, P);
a_est_spatialKhype_im = zeros(nr,nc,P);
a_est_spatialKhype_im(1:nr, 1:nr, :)          = a_est_spatialKhype1;
a_est_spatialKhype_im(1:nr, (nr+1):(2*nr), :) = a_est_spatialKhype2;
a_est_spatialKhype_im(1:nr, end-nr+1:end, :)  = a_est_spatialKhype3;
a_est_spatialKhype = reshape(a_est_spatialKhype_im,nr*nc,P)';


% compute RMSE
mse_TV_pixelwise1 = reshape(mse_TV_pixelwise1', nr, nr, []);
mse_TV_pixelwise2 = reshape(mse_TV_pixelwise2', nr, nr, []);
mse_TV_pixelwise3 = reshape(mse_TV_pixelwise3', nr, nr, []);

mse_TV_pixelwise_im = zeros(nr,nc);
mse_TV_pixelwise_im(1:nr, 1:nr, :)          = mse_TV_pixelwise1;
mse_TV_pixelwise_im(1:nr, (nr+1):(2*nr), :) = mse_TV_pixelwise2;
mse_TV_pixelwise_im(1:nr, end-nr+1:end, :)  = mse_TV_pixelwise3;
rmse_r_TV = reshape(mse_TV_pixelwise_im,nr*nc,1)';
rmse_r_TV = sqrt(sum(rmse_r_TV));


beta_TV1 = reshape(beta_TV1', nr, nr, L);
beta_TV2 = reshape(beta_TV2', nr, nr, L);
beta_TV3 = reshape(beta_TV3', nr, nr, L);
beta_TV_im = zeros(nr,nc,L);
beta_TV_im(1:nr, 1:nr, :)          = beta_TV1;
beta_TV_im(1:nr, (nr+1):(2*nr), :) = beta_TV2;
beta_TV_im(1:nr, end-nr+1:end, :)  = beta_TV3;
beta_TV = reshape(beta_TV_im,nr*nc,L)';


time_spatialKhype = toc;









%%
% =========================================================================
% BMUAN


% Norm of modeling errors
% sigma2_epsi = (1e-4) * norm(r,'fro')^2/N; 
sigma2_epsi = (1e-5) * norm(r,'fro')^2/N; % 40 db of Signal-to-modeling errors ratio
% sigma2_epsi = (1e-3) * norm(r,'fro')^2/N; 
% sigma2_epsi = (1e-2) * norm(r,'fro')^2/N; 
% sigma2_epsi = (1e-1) * norm(r,'fro')^2/N; 

tic
[a_est_mscale,rmse_r_SPPX]=BMUAN(r,M,nr,nc,sigma2_epsi);
timeMscale = toc;

% fprintf('\n\n mu0: %f, mu1: %f, mu2: %f \n\n',mu_0, mu_1, mu_2)





%% ========================================================================
% =========================================================================
% =========================================================================
% OTHER METHODS

%% Halimis TIP16 method

flag_sto_halimi = -1;
tic
% [a_NL] = Unmix_CDA_NL_TIP_v1(M,Yim,0,-1,0);
% [a_NL] = Unmix_CDA_NL_TIP_v1(M,Yim,0,-1,-1);
% [a_NL, sigma2t, gam_NLt_halimi, s, w, c_Ev_halimi, info, Cost_fun] = Unmix_CDA_NL_TIP_v1(M,Yim,0,-1,-1);
% [a_NL,~,gam_NLt_halimi,~,~,c_Ev_halimi,~,~] = Unmix_CDA_NL_TIP_v1(M,Yim,0,-1,-1);
[a_NL,~,gam_NLt_halimi,~,~,c_Ev_halimi,~,~] = Unmix_CDA_NL_TIP_v1(M,Yim,0,-1,flag_sto_halimi);
a_NL = a_NL(:,:,end);
time_NL_halimi = toc;
a_NL_im = reshape(a_NL',nr,nc,P);  

M_NL_halimi = [];
for i=1:P, for j=i+1:P, M_NL_halimi = [M_NL_halimi  M(:,i).*M(:,j)]; end, end
M_NL_halimi = [sqrt(2)*M_NL_halimi M.^2]; % LxD 
Y_rec_halimi  = (M*a_NL) .*(c_Ev_halimi(:,end)*ones(1,L))' ...
    + (M_NL_halimi*gam_NLt_halimi(:,:,end)) .*(c_Ev_halimi(:,end).^2*ones(1,L))';

rmse_r_NL_halimi = sqrt(sum(sum((r - Y_rec_halimi).^2)));



%% Rita's TIP17 method (OK)
tic
lambdaVVal = 5;
% strVVal = 'Transformable'; % 'full';
strVVal = 'Separable';
[a_NDU,F_NDU_rita] = nonlinearU_vectValKernels(Yim,M,lambdaVVal,strVVal);
a_NDU_im = reshape(a_NDU',nr,nc,P);  %figure, imagesc(a_NDU_im(:,:,2)), colormap jet
a_NDU    = reshape(a_NDU_im,nr*nc,P)';
time_vectValKernel = toc;

rmse_r_vecVal_rita = sqrt(sum(sum((r - (M*a_NDU+F_NDU_rita)).^2)));




%% ========================================================================
% =========================================================================
% =========================================================================

fprintf('\n\n TIME \n')
fprintf('FCLS.............: %f \n',time_fcls)
fprintf('KHype............: %f \n',time_khype)
fprintf('K-Hype-TV........: %f \n',time_spatialKhype)
fprintf('Halimi NL TIP16..: %f \n',time_NL_halimi)
fprintf('VectorValKernel..: %f \n',time_vectValKernel)
fprintf('BMUA-N...........: %f \n',timeMscale)


fprintf('\n\n Reconstruction errors: \n') 
fprintf('FCLS............: RMSE_R = %f \n', rmse_r_FCLS        /sqrt(L*N) )
fprintf('K-Hype..........: RMSE_R = %f \n', rmse_r_KHYPE       /sqrt(L*N) )
fprintf('K-Hype-TV.......: RMSE_R = %f \n', rmse_r_TV          /sqrt(L*N) )
fprintf('Halimi TIP16....: RMSE_R = %f \n', rmse_r_NL_halimi   /sqrt(L*N) )
fprintf('Rita TIP17......: RMSE_R = %f \n', rmse_r_vecVal_rita /sqrt(L*N) )
fprintf('BMUA-N..........: RMSE_R = %f \n', rmse_r_SPPX        /sqrt(L*N) )






%%

A_FCLS_im         = reshape(A_FCLS,nr,nc,[]);
a_est_khype_im    = reshape(a_est_khype',nr,nc,[]);
a_est_khypeTV_im  = reshape(a_est_spatialKhype',nr,nc,[]);
a_est_mscale_im   = reshape(a_est_mscale',nr,nc,[]);

idx_ems2 = 1:P;
P_reduced = length(idx_ems2); % 4
fh = figure; 
[ha, pos] = tight_subplot(5, P_reduced, 0.01, 0.1, 0.1);

maxval = 1;
for i=1:P_reduced
    kk = 0;
    axes(ha(i+kk*P_reduced)); kk = kk+1;
    imagesc(a_est_khype_im(:,:,idx_ems2(i)),[0 maxval])
    set(gca,'ytick',[],'xtick',[])
    axes(ha(i+kk*P_reduced)); kk = kk+1;
    imagesc(a_est_khypeTV_im(:,:,idx_ems2(i)),[0 maxval])
    set(gca,'ytick',[],'xtick',[])
    
    axes(ha(i+kk*P_reduced)); kk = kk+1;
    imagesc(a_NL_im(:,:,idx_ems2(i)),[0 maxval])
    set(gca,'ytick',[],'xtick',[])
    
    axes(ha(i+kk*P_reduced)); kk = kk+1;
    imagesc(a_NDU_im(:,:,idx_ems2(i)),[0 maxval])
    set(gca,'ytick',[],'xtick',[])
    
    axes(ha(i+kk*P_reduced)); kk = kk+1;
    imagesc(a_est_mscale_im(:,:,idx_ems2(i)),[0 maxval])
    set(gca,'ytick',[],'xtick',[])
    
end

% set(fh, 'Position', [0 0 650 700])
axes(ha(1));
title('Asphalt','interpreter','latex')
axes(ha(2));
title('Tree','interpreter','latex')
axes(ha(3));
title('Ground','interpreter','latex')

fontSizeNum = 12;

kk = 0;
axes(ha(kk*P_reduced+1)); kk = kk+1;
ylabel('K-Hype','interpreter','latex')
axes(ha(kk*P_reduced+1)); kk = kk+1;
ylabel('K-Hype TV','interpreter','latex')
axes(ha(kk*P_reduced+1)); kk = kk+1;
ylabel('CDA-NL','interpreter','latex')
axes(ha(kk*P_reduced+1)); kk = kk+1;
ylabel('NDU','interpreter','latex')
axes(ha(kk*P_reduced+1)); kk = kk+1;
ylabel('BMUA-N','interpreter','latex')

colormap(jet)


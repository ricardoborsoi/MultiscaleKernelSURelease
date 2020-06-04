% =========================================================================
% 
% This code contains the example with the Jasper Ridge image in the following paper:
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

load('DATA/jasperRidge2_R198.mat')

im = reshape(Y', 100, 100, 198 )/6000;
[nr,nc,L] = size(im);
N = nr*nc;
r = reshape(im,nr*nc,L)';
Yim = im;



% figure;
% clear rgb
% v = [25,15,5];
% rgb(:,:,1) = imadjust(rescale(im(:,:,v(1)),0,1));
% rgb(:,:,2) = imadjust(rescale(im(:,:,v(2)),0,1));
% rgb(:,:,3) = imadjust(rescale(im(:,:,v(3)),0,1));
% imshow(rgb) % display RGB image
% set(gca,'ytick',[],'xtick',[])



% Extract EMs 
P = 4;
% [M0, V, U, Y_bar, endm_proj, Y_proj] = find_endm(r,P,'vca');
load('DATA/M0_jasperRidge.mat')
M = M0;



% =========================================================================
%  FCLS

tic
A_FCLS = FCLSU(r,M);
time_fcls = toc;

rmse_r_FCLS = sqrt(norm(r-M*A_FCLS','fro')^2);



%%
% =========================================================================
%  K-Hype

% Regualrization parameter 
C = 500; % \mu in the paper = 1/C

tic
[a_est_khype,beta_khype,rmse_r_KHYPE] = KHype(r,M,C);
time_khype = toc;



%%
% =========================================================================
% TV

% Functional regularization
C = 500;
% Reg. parameter for spatial corr.
lambda_sp = 0.1;

tic
[a_est_spatialKhype,rmse_r_TV,~,~]=KhypeTV_wrapper2(r,N,P,nr,nc,M,L,lambda_sp,C);
time_spatialKhype = toc;







%%
% =========================================================================
% BMUAN


% Norm of modeling errors
% sigma2_epsi = (1e-4) * norm(r,'fro')^2/N; 
% sigma2_epsi = (1e-5) * norm(r,'fro')^2/N; % 40 db of Signal-to-modeling errors ratio
% sigma2_epsi = (1e-2) * norm(r,'fro')^2/N; 
% sigma2_epsi = (1e-3) * norm(r,'fro')^2/N; 
sigma2_epsi = (1e-4) * norm(r,'fro')^2/N; 

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



%%


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





%%


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
title('Vegetation','interpreter','latex')
axes(ha(2));
title('Road','interpreter','latex')
axes(ha(3));
title('Ground','interpreter','latex')
axes(ha(4));
title('Water','interpreter','latex')

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







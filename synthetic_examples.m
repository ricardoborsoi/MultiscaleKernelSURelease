% =========================================================================
% 
% This code contains the synthetic examples in the following paper:
% 
%    A Blind Multiscale Spatial Regularization Framework for Kernel-Based Spectral Unmixing 
%    R.A. Borsoi, T. Imbiriba, J.C.M. Bermudez, C. Richard.
%    IEEE Transactions on Image Processing, 2020.
% =========================================================================

clear
close all
clc
warning off;


disp('Using polynomial kernel...')


addpath(genpath('other_methods'))
addpath(genpath('DATA'))
addpath(genpath('utils'))
addpath(genpath('BMUAN'))


rng(10, 'twister') % reproducible



% --------------------- OPTIONS ---------------------
% Noise SNR (20 or 30 dB) 
SNR = 20; 30;


% abundance maps ----
ab_maps = 'im1_squares';
% ab_maps = 'im2_spatial'; 

% nonlinearity model (2 or 3) ----
% nlmix_fcn = 2; nlmm_str = 'blmm'; % bilinear mixing model
nlmix_fcn = 3; nlmm_str = 'pnmm'; % post-nonlinear mixing model


% estimate EMs ----
flag_estimate_EM = false; eeaStr = 'trueEMs'; % use known EMs
% flag_estimate_EM = true; eeaStr = 'VCA';  % use estimated EMs (VCA)


% Gaussian kernel bandwidth (not used) ----
% par = 2;

% =========================================================================
% generate synthetic image

% load abundance maps
if strcmp(ab_maps,'im1_squares')
    load('DATA/abundance_squares_3em.mat')
elseif strcmp(ab_maps,'im2_spatial')
    load('DATA/abundance_spatial_3em.mat')
else
    error('Unknown abundance maps selected!')
end

nr = size(A_cube,1);
nc = size(A_cube,2);
P  = size(A_cube,3);
N  = nr*nc;

% THIS ORDERING IS OK ---------------
% convert matrix to image
conv2im  = @(A)  reshape(A',nr,nc,P);
% convert image to matrix
conv2mat = @(A)  reshape(A,nr*nc,P)';
% -----------------------------------

% convert to matrix
Ath = conv2mat(A_cube);


% load endmembers ---------------------------------------------------------
[M,namesM] = load_endmembers(P);
L = size(M,1);

% generate noiseless image according to the nonlinearity model
% [Y, ~] = hypermix(M, nr*nc, 2, Ath); % Bilinear mixing model
% [Y, ~] = hypermix(M, nr*nc, 3, Ath); % Post nonlinear mixing model
[Y, ~] = hypermix(M, nr*nc, nlmix_fcn, Ath);

% add noise
pw_signal = norm(Y,'fro')^2/(nr*nc*L);
pw_noise  = pw_signal/(10^(SNR/10));
std_noise = sqrt(pw_noise);
noise     = std_noise*randn(L,nr*nc);
r = Y + noise;

% write the image over other variables too
Y = r;
Yim = reshape(Y', nr, nc, L);




% --------------------------------------------------------------
% Endmember initialization 

% save true endmember matrix
Mth = M;

% load pre-saved estimated endmembers if available, or extract it using VCA
EM_str = ['M0_' ab_maps '_SNR' num2str(SNR) '_' nlmm_str];
if exist([EM_str '.mat'],'file') == 2
    load(EM_str)
else
    [M0, V, U, Y_bar, endm_proj, Y_proj] = find_endm(r,P,'vca');
    save(EM_str,'M0')
end


% Sort M0 with respect to real/desired EM signatures to ease the comparison 
% of estimated abundance maps
id = zeros(P,1);
for k = 1:P
    for l = 1:P
        s(l) = 180*acos( (Mth(:,k).')*M0(:,l) /(norm(Mth(:,k))*norm(M0(:,l))) )/pi; 
    end
    [~, id(k)] = min(s);
end
M0 = M0(:,id);


% -------------------------------
if flag_estimate_EM
    M = M0;
    disp('Using VCA...')
else
    M = Mth;
    disp('Using TRUE EMs...')
end







% =========================================================================
%  FCLS
tic
A_FCLS = FCLSU(r,M);
time_fcls = toc;

rmse_r_FCLS = sqrt(norm(r-M*A_FCLS','fro')^2);


%%
% =========================================================================
%  K-Hype

% Select parameters
if flag_estimate_EM
    % VCA
    if strcmp(ab_maps,'im1_squares')
        if SNR == 20
            if nlmix_fcn == 2
                C = 1;
            elseif nlmix_fcn == 3
                C = 10;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 200;
            elseif nlmix_fcn == 3
                C = 200; 
            end
        else
            error('No parameters for K-Hype!')
        end

    elseif strcmp(ab_maps,'im2_spatial')
        if SNR == 20
            if nlmix_fcn == 2
                C = 200;
            elseif nlmix_fcn == 3
                C = 50;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 500;
            elseif nlmix_fcn == 3
                C = 500;
            end
        else
            error('No parameters for K-Hype!')
        end
    end
    
    
else % TRUE EMs
    if strcmp(ab_maps,'im1_squares')
        if SNR == 20
            if nlmix_fcn == 2
                C = 10;
            elseif nlmix_fcn == 3
                C = 200;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 50;
            elseif nlmix_fcn == 3
                C = 1000; 
            end
        else
            error('No parameters for K-Hype!')
        end

    elseif strcmp(ab_maps,'im2_spatial')
        if SNR == 20
            if nlmix_fcn == 2
                C = 50;
            elseif nlmix_fcn == 3
                C = 1000;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 500;
            elseif nlmix_fcn == 3
                C = 1000;
            end
        else
            error('No parameters for K-Hype!')
        end
    end
end


tic
[a_est_khype,beta_khype,rmse_r_KHYPE] = KHype(r,M,C);
time_khype = toc;



%%
% =========================================================================
% TV

% Select parameters
if flag_estimate_EM
    % VCA
    if strcmp(ab_maps,'im1_squares')
        if SNR == 20
            if nlmix_fcn == 2
                C = 1;
                lambda_sp = 0;
            elseif nlmix_fcn == 3
                C = 100;
                lambda_sp = 0.1;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 500;
                lambda_sp = 0;
            elseif nlmix_fcn == 3
                C = 100;
                lambda_sp = 0;
            end
        else
            error('No parameters for TV!')
        end

    elseif strcmp(ab_maps,'im2_spatial')
        if SNR == 20
            if nlmix_fcn == 2
                C = 100;
                lambda_sp = 0;
            elseif nlmix_fcn == 3
                C = 500;
                lambda_sp = 0.1;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 500;
                lambda_sp = 0;
            elseif nlmix_fcn == 3
                C = 500;
                lambda_sp = 0;
            end
        else
            error('No parameters for TV!')
        end
    end

    
else % TRUE EMs
    if strcmp(ab_maps,'im1_squares')
        if SNR == 20
            if nlmix_fcn == 2
                C = 50;
                lambda_sp = 0.1;
            elseif nlmix_fcn == 3
                C = 1000;
                lambda_sp = 0.1;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 50;
                lambda_sp = 0.01;
            elseif nlmix_fcn == 3
                C = 1000;
                lambda_sp = 0.01;
            end
        else
            error('No parameters for TV!')
        end

    elseif strcmp(ab_maps,'im2_spatial')
        if SNR == 20
            if nlmix_fcn == 2
                C = 500;
                lambda_sp = 0.1;
            elseif nlmix_fcn == 3
                C = 1000;
                lambda_sp = 0.1;
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                C = 500;
                lambda_sp = 0.01;
            elseif nlmix_fcn == 3
                C = 1000;
                lambda_sp = 0.01;
            end
        else
            error('No parameters for TV!')
        end
    end
end



tic
[a_est_spatialKhype,rmse_r_TV,~,~]=KhypeTV_wrapper2(r,N,P,nr,nc,M,L,lambda_sp,C);
time_spatialKhype = toc;











%%
% =========================================================================
% BMUAN

% Norm of modeling errors 
% sigma2_epsi = (1e-4) * norm(Y,'fro')^2/N; 
% sigma2_epsi = (1e-5) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
% sigma2_epsi = (1e-6) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
sigma2_epsi = (1e-8) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio
% sigma2_epsi = (0e-8) * norm(Y,'fro')^2/N; % Many db of Signal-to-modeling errors ratio

tic
[a_est_mscale,rmse_r_SPPX]=BMUAN(r,M,nr,nc,sigma2_epsi);
timeMscale = toc;





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
ErrComput(Ath, a_NL);
time_NL_halimi = toc;

M_NL_halimi = [];
for i=1:P, for j=i+1:P, M_NL_halimi = [M_NL_halimi  M(:,i).*M(:,j)]; end, end
M_NL_halimi = [sqrt(2)*M_NL_halimi M.^2]; % LxD 
Y_rec_halimi  = (M*a_NL) .*(c_Ev_halimi(:,end)*ones(1,L))' ...
    + (M_NL_halimi*gam_NLt_halimi(:,:,end)) .*(c_Ev_halimi(:,end).^2*ones(1,L))';

rmse_r_NL_halimi = sqrt(sum(sum((Y - Y_rec_halimi).^2)));



%% Rita's TIP17 method (OK)


% Select parameters
if flag_estimate_EM
    error('No parameters for VCA here!')
    
else % TRUE EMs
    if strcmp(ab_maps,'im1_squares')
        if SNR == 20
            if nlmix_fcn == 2
                lambdaVVal = 5; 
                strVVal = 'Transformable'; % 'full';
            elseif nlmix_fcn == 3
                lambdaVVal = 5000;
                strVVal = 'Transformable'; % 'full';
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                lambdaVVal = 5; 
                strVVal = 'Transformable'; % 'full';
            elseif nlmix_fcn == 3
                lambdaVVal = 5000;
                strVVal = 'Transformable'; % 'full';
            end
        else
            error('No parameters here!')
        end

    elseif strcmp(ab_maps,'im2_spatial')
        if SNR == 20
            if nlmix_fcn == 2
                lambdaVVal = 5;
                strVVal = 'Transformable'; 
            elseif nlmix_fcn == 3
                lambdaVVal = 5000;
                strVVal = 'Transformable'; % 'full';
            end
        elseif SNR == 30
            if nlmix_fcn == 2
                lambdaVVal = 5;
                strVVal = 'Transformable';
            elseif nlmix_fcn == 3
                lambdaVVal = 5000;
                strVVal = 'Transformable'; % 'full';
            end
        else
            error('No parameters here!')
        end
    end
end

tic
[a_NDU,F_NDU_rita] = nonlinearU_vectValKernels(Yim,M,lambdaVVal,strVVal);
a_NDU_im = reshape(a_NDU',nr,nc,P); 
a_NDU    = reshape(a_NDU_im,nr*nc,P)';
ErrComput(Ath, a_NDU);
time_vectValKernel = toc;

rmse_r_vecVal_rita = sqrt(sum(sum((Y - (M*a_NDU+F_NDU_rita)).^2)));





%% ========================================================================
% =========================================================================
% =========================================================================
% Results


clc

% RMSE calculation 
[RMSE0a, std0a] = ErrComput(Ath, A_FCLS');
[RMSE1, std1] = ErrComput(Ath, a_est_khype);
[RMSE2, std2] = ErrComput(Ath, a_est_spatialKhype);
[RMSE3, std3] = ErrComput(Ath, a_est_mscale);
RMSE5 = ErrComput(Ath, a_NL);
RMSE6 = ErrComput(Ath, a_NDU);


fprintf('\n\n RMSE \n')
fprintf('FCLS.............: %f \n',RMSE0a)
fprintf('KHype............: %f \n',RMSE1)
fprintf('K-Hype-TV........: %f \n',RMSE2)
fprintf('Halimi TIP16.....: %f \n',RMSE5)
fprintf('Rita TIP17.......: %f \n',RMSE6)
fprintf('BMUA-N...........: %f \n',RMSE3)

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

% ====================================================
A_FCLSU_im  = reshape(A_FCLS,nr,nc,P);
A_KHYPE_im  = reshape(a_est_khype',nr,nc,P);
A_TV_im     = reshape(a_est_spatialKhype',nr,nc,P);
A_NL_im     = reshape(a_NL',nr,nc,P);
A_NDU_im    = reshape(a_NDU',nr,nc,P);
A_BMUAN_im  = reshape(a_est_mscale',nr,nc,P);





fh = figure;
[ha, pos] = tight_subplot(6, P, 0.01, 0.1, 0.1);
for i=1:P
    kk = 0;
    
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_cube(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_KHYPE_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_TV_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
    
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_NL_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_NDU_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
    
    axes(ha(i+kk*P)); kk = kk+1;
    imagesc(A_BMUAN_im(:,:,i),[0 1])
    set(gca,'ytick',[],'xtick',[])%, axis image
end
% set(fh, 'Position', [0 0 650 700])
axes(ha(1));
title('EM #1','interpreter','latex')
axes(ha(2));
title('EM #2','interpreter','latex')
axes(ha(3));
title('EM #3','interpreter','latex')

fontSizeNum = 10;

axes(ha(1));
ylabel('True','interpreter','latex','fontsize',fontSizeNum)
axes(ha(4));
ylabel('K-Hype','interpreter','latex','fontsize',fontSizeNum)
axes(ha(7));
ylabel('K-Hype TV','interpreter','latex','fontsize',fontSizeNum)

axes(ha(10));
ylabel('CDA-NL','interpreter','latex','fontsize',fontSizeNum)
axes(ha(13));
ylabel('NDU','interpreter','latex','fontsize',fontSizeNum)
axes(ha(16));
ylabel('BMUA-N','interpreter','latex','fontsize',fontSizeNum)
colormap(jet)



%%
% 
% % print('example2/abundances_houston','-dpng')
% fnameStr = ['examples/abundances_polyK_' eeaStr '_' ab_maps '_' nlmm_str '_SNR' num2str(SNR) '.pdf'];
% print(fnameStr,'-dpdf')
% system(['pdfcrop ' fnameStr ' ' fnameStr]) 


%fprintf('\n\n mu0: %f, mu1: %f, mu2: %f \n\n',mu_0, mu_1, mu_2)









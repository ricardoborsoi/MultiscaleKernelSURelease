function [alphat sigma2t Outliers s w C_EV info Cost_fun] = Unmix_CDA_ME_TIP_v1(MPlus,y,Bands,L0,Init,Illum);

%% Unmix_CDA_ME_TIP(MPlus,y,Bands,L0,Init);
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  CDA-ME
%  Version (April, 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% MPlus  : endmembers of size (L x R)
%%% y      : pixels of size (row x col x N)
%%% Bands  : Index of the considered bands (used to relax the
%%%          spectral correlation from the removed water absorption bands)
%%% L0     : Number of spectral bands before removing bands >= L
%%% Init   : if =0  Initialze randomly
%%%          else   Initialization with Sunsal for the abundances
%%%                 and Hysime for the noise variance
%%% Illum  :if =0  Impose abundance sum-to-one (c=1)
%%%        :else   Estimate the illumation coefficient c
%%%
%% --------------- Outputs --------------------------------------------
%%% alphat   : Chain of abundance values  (R x N x Iter/10)
%%% sigma2t  : Chain of the noise variance  (L x Iter/10)
%%% Outliers : Residual components or ME (L x N) with N= row x col
%%% s        : Variances or energies of the residual components  (N x 1)
%%% w        : Auxiliary  variable of size (N x 1)
%%% C_EV     : Chain of the illumination coefficient (N x Iter)
%%% info     : Function ended because
%%%            info(1)=1: the criterion associated with (a.w.) the cost function is
%%%            satisfied
%%%            info(2)=1: the criterion a.w. the abundances was satisfied
%%%            info(3)=1: the criterion a.w. the outliers was satisfied
%%%            info(4)=1: the maximum number of iterations was attained
%% --------------------------------------------------------------------
%% -------------------------------------------------------------------------
%
% Copyright (April, 2016):        Abderrahim Halimi (a.halimi@hw.ac.uk)
%
% CDA-ME is distributed under the terms of
% the GNU General Public License 3.0.
%
% Permission to use, copy, modify, and distribute this software for
% any purpose without fee is hereby granted, provided that this entire
% notice is included in all copies of any software which is or includes
% a copy or modification of this software and in all copies of the
% supporting documentation for such software.
% This software is being provided "as is", without any express or
% implied warranty.  In particular, the authors do not make any
% representation or warranty of any kind concerning the merchantability
% of this software or its fitness for any particular purpose."
% ---------------------------------------------------------------------



%% Constants
[row col L] = size(y);
N           = row*col;
N1          = (row+2)*(col+2);
Y_bloc      = reshape(y,N,L)'; % L x N
R           = size(MPlus,2);
a           = 0.52;   % Coupling parameter controling the degree of spatial smoothness
StpKmax     = 500;   % Maximum number of iterations
Stpfunc     = 10^(-5);  % Threshold of the cost function criterion
StpX        = 10^(-6);  % Threshold of the abundance criterion
StpXG       = 10^(-11); % Threshold of the outlier criterion

if(L~=size(MPlus,1))
    error('mixing matrix MPlus and data set y are inconsistent')
end
if(L~=length(Bands))
    error('Bands is  inconsistent with MPlus or y')
end
if(L0<L)
    error('L0 is inconsistent with MPlus or y')
end

% Squred exponential covariance function
Hdir        = toeplitz(exp(-[1:L0].^2/100^2));Hdir=Hdir(Bands,Bands);
% Hdir        = toeplitz(exp(-[1:L0].^2/(length(Bands)/2)^2));Hdir=Hdir(Bands,Bands);

Kinv        = inv(Hdir);
Hdir        = chol(Kinv);
Hinv        = inv(Hdir);
[VecK,ValK] = svd(Kinv);
ValK        = repmat(diag(ValK),1,N); % L x N


%% Initialization
Outliers  = zeros(L,N);
m_Outl    = Outliers;
[w Rn]    = estNoise(Y_bloc,'additive','off');
sigma2    = diag(Rn);                   % Noise variances
bHy       = (100* sigma2+ 1) .*sigma2 ; % Hyperparameter a.w. Sigma2
aHy       = (100* sigma2+ 2);           % Hyperparameter a.w. Sigma2
s         = 10^(-7)*ones(N1,1);
w         = ones(N,1);
sbi_j     = reshape(s,row+2,col+2);
s         = reshape(sbi_j(2:end-1,2:end-1),N,1);
wbi_j     = reshape(w,row,col);
wbi_j1    = circshift(wbi_j, [0 -1]);    %LxN  gauche i j+1
wbi1_j1   = circshift(wbi_j, [-1 -1]);   %LxN  haut gauche i+1 j+1
wbi1_j    = circshift(wbi_j, [-1 0]);    %LxN  haut i+1 j 
alpha     = sunsal(MPlus,Y_bloc ,'lambda',0,'ADDONE','no','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
C_EV      = sum(alpha,1)';  % Illumination coefficients
 
% sum(C_EV>2)
% (sum(C_EV>2) > N/10) 
% pause
if(sum(C_EV>2) > N/10) 
    sigCEV      = 0.0000001; % Hyperparameter a.w. C_EV
%     sigCEV      = 0.0000000001; %For Moffett
else
    sigCEV      = 0.01;    % Hyperparameter a.w. C_EV
end
 
alpha    = sunsal(MPlus,Y_bloc ,'lambda',0,'ADDONE','yes','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
if(Init ==0)
    alpha   = 1/R*ones(R,N);
    sigma2  = 10^(-2)*ones(L,1);
    C_EV = ones(N,1);
end

%% iterative algorithm
info      = zeros(4,1);
condition = 0;compt=1;
while condition == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Abundance (alpha of size R x N)
    M2        = MPlus./repmat(sqrt(sigma2),1,R);
    Y_Lin2    = (Y_bloc-m_Outl)./repmat(sqrt(sigma2),1,N);
    alphaIni  = alpha;
    lambda    = 0; 
    alpha     = sunsal(M2,Y_Lin2./repmat(C_EV(:,compt)',L,1),'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
        'AL_iters',200,'TOL', 1e-4, 'verbose','no','X0',abs(alphaIni)./repmat(sum(abs(alphaIni),1),R,1));
    m_Lin     = MPlus*alpha.*(C_EV(:,compt)*ones(1,L))'; % LxN
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Outliers  (Outliers of size L x N)
    OutliersIni = Outliers;
    Y_EV        = Y_bloc-m_Lin; % Lx N
    [Vec,Val]   = svd(Hinv' * diag(1./sigma2) * Hinv);
    Vec         = Hinv*Vec;
    Val         = diag(Val); % Lx1
    vali        = 1./(Val*ones(1,N)+repmat(1./s(:,compt)',L,1)); % LxN
    mu_EV22     = Y_EV./repmat(sigma2,1,N); % L x N
    mu_EV22     = Vec' *mu_EV22;
    Outliers    = Vec * (mu_EV22.*vali)  ;     % LxN
    m_Outl      = Outliers ;  % LxN
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Noise variance (sigma of size Lx1)
    sigma2 = max((sum( (Y_bloc - m_Outl - m_Lin).^2,2) + 2*bHy)./(N+2 +2*aHy),eps);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Outlier energies (s of size Nx1)
    bw_IG  = wbi_j + wbi_j1 + wbi1_j1 + wbi1_j; %rowxcol
    bw_IG  = a*bw_IG; %rowxcol
    gmat   = reshape(2* sum((m_Outl ).^2,1)',row,col); % rowxcol
    sbi_j(2:end-1,2:end-1)    = (bw_IG + gmat   ) / (4*a+L/2+1);  % row-2 xcol-2
    
    sbi_j1  = circshift(sbi_j, [0 1]);   %LxN   i j+1
    sbi1_j1 = circshift(sbi_j, [1 1]);   %LxN    i+1 j+1
    sbi1_j  = circshift(sbi_j, [1 0]);   %LxN   i+1 j
    s(:,compt+1)       = reshape(sbi_j(2:end-1,2:end-1),N,1);
    s(:,compt+1)       = max(s(:,compt+1) ,10^(-100));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Auxilliary variables (w of size Nx1)
    bs_IG   = 1./sbi_j + 1./sbi_j1 + 1./sbi1_j1 + 1./sbi1_j; % row x col
    wbi_j   = (4*a-1)./(a* bs_IG(2:end-1,2:end-1)) ;   %row-1 x col-1
    wbi_j1  = circshift(wbi_j, [0 -1]);  %LxN    i j+1
    wbi1_j1 = circshift(wbi_j, [-1 -1]);   %LxN    i+1 j+1
    wbi1_j  = circshift(wbi_j, [-1 0]);  %LxN    i+1 j
    w(:,compt+1)       = reshape(wbi_j,N,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%  Illumination coefficients (C_EV of size Nx1)
    if(Illum==0)
        C_EV(:,compt+1)  = ones(N,1);
    else
        m_Ling           = MPlus*alpha; % LxN
        m_Lin2           = m_Ling./repmat(sqrt(sigma2),1,N);
        Y_bloc2          = (Y_bloc - m_Outl)./repmat(sqrt(sigma2),1,N);
        x1               = sum(m_Lin2.^2,1)'; % N x 1
        x2               = sum(m_Lin2.*Y_bloc2,1)'; % N x 1
        sigc             = 1./(x1+1/sigCEV );
        C_EV(:,compt+1)  = max(0.2, sigc.* ( x2+1/sigCEV )); % N x 1
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Test of convergence
    m_Lin  = MPlus*alpha.*(C_EV(:,compt+1)*ones(1,L))'; % LxN
    bw_IG  = wbi_j + wbi_j1 + wbi1_j1 + wbi1_j; %r x c
    bw_IG  = a*bw_IG; %r x c
    
    if mod(compt,10) == 1
        indCov             = (compt-1)/10+1;
        alphat(:,:,indCov) = alpha;
        sigma2t(:,indCov)  = sigma2;
        Cost_fun(indCov)   =  N/2*log(sigCEV)  -4*a*sum(log(bw_IG(:))) +  L/2*log(2*pi) ...
            + sum((N/2+1+aHy).*(log(sigma2)-log(eps))) + sum((0.5*sum((Y_bloc - m_Outl - m_Lin).^2,2)+bHy)./sigma2) ...
            + sum((C_EV(:,compt+1)-1).^2)/2/sigCEV +  sum((4*a+1+L/2).* (log(s(:,compt+1))+100) ) ...
            + sum( (0.5*sum((VecK'*Outliers).^2.*ValK,1)' +  bw_IG(:))./s(:,compt+1) ) ;
        
        cf1(indCov)   = N/2*log(sigCEV)  -4*a*sum(log(bw_IG(:))) +  L/2*log(2*pi) ...
            + sum((N/2+1+aHy).*(log(sigma2)-log(eps)));
        cf2(indCov)   = sum((4*a+1+L/2).* (log(s(:,compt+1))+100) );
         
        nh(indCov)         = norm(alpha-alphaIni,'fro');
        maxDiff(indCov)    = max(abs(alpha(:)-alphaIni(:)));
        nx(indCov)         = StpX + norm(alphaIni,'fro');
        nhG(indCov)        = norm(Outliers-OutliersIni,'fro');
        maxDiffG(indCov)   = max(abs(Outliers(:)-OutliersIni(:)));
        nxG(indCov)        = StpX + norm(OutliersIni,'fro');
        
        if(compt>10)
            ng = abs(Cost_fun(indCov)-Cost_fun(indCov-1));
            if(ng <= Stpfunc*abs(Cost_fun(indCov-1)))  info(1)=1; end
        end
        if(nh(indCov) <= StpX*nx(indCov))    info(2)=1; end
        if(nhG(indCov) <= StpXG*nxG(indCov))  info(3)=1; end
        if(compt>StpKmax)    info(4)=1; end
        
        if(sum(info)>0)
            condition =1;
        end
         
        
        fprintf(' iter = %f \n',compt)
    end
    compt = compt+1;
    
end

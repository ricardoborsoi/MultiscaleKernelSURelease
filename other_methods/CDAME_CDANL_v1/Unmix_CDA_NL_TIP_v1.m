function [alphat sigma2t gam_NLt s w C_EV info Cost_fun] = Unmix_CDA_NL_TIP_v1(MPlus,y,Posi,Init,Illum);

%% Unmix_CDA_NL_TIP(MPlus,y,Posi,Init);
%  Paper :  A. Halimi, P. Honeine, J. Bioucas-Dias, "Hyperspectral Unmixing
%          in Presence of Endmember Variability, Nonlinearity or Mismodelling
%          Effects", IEEE Trans. Image Process., 2016.
%  Code  :  CDA-NL
%  Version (April, 2016) by Abderrahim Halimi (a.halimi@hw.ac.uk)
%  For any comments contact the author
%
%% --------------- Inputs --------------------------------------------
%%% MPlus  : endmembers of size (L x R)
%%% y      : pixels of size (row x col x N)
%%% Posi   : Positivity constraint on Gamma
%%%          Posi=0: without positivity constraint (gam \in R)
%%%          else  : with positivity constraint (gam >0),
%%% Init   : if =0  Initialze randomly
%%%          else   Initialization with Sunsal for the abundances
%%%                 and Hysime for the noise variance
%%% Illum  :if =0  Impose abundance sum-to-one (c=1)
%%%        :else   Estimate the illumation coefficient c
%%%
%% --------------- Outputs --------------------------------------------
%%% alphat   : Chain of abundance values  (R x N x Iter/10)
%%% sigma2t  : Chain of the noise variance  (L x Iter/10)
%%% gam_NLt  : Chain of nonlinear coefficients (D x N x Iter/10) with N= row x col
%%% s        : Variances or energies of the nonlinear coefficients (N x 1)
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
% CDA-NL is distributed under the terms of
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
D           = R*(R+1)/2;
a           = 0.52;    % Coupling parameter controling the degree of spatial smoothness
StpKmax     = 500;
Stpfunc     = 10^(-5); % Threshold of the cost function criterion
StpX        = 10^(-6); % Threshold of the abundance criterion
StpXG       = 10^(-6); % Threshold of the nonlinear coefficients criterion

if(L~=size(MPlus,1))
    error('mixing matrix MPlus and data set y are inconsistent')
end

% matrix of endmember interactions
M_NL    = [];
for i=1:R
    for j=i+1:R
        M_NL   = [M_NL  MPlus(:,i).*MPlus(:,j)];
    end
end
M_NL      = [sqrt(2)*M_NL MPlus.^2]; % LxD
small     = 10^(-6);
[EVe,EVa] = svd(M_NL'*M_NL +small*eye(D));
EVa       = diag(EVa);

%%% Initialization
gam_NL  =  10^(-2)*rand(D,N);
[w Rn]  = estNoise(Y_bloc,'additive','off');
sigma2  = diag(Rn);                   % Noise variances
bHy     = (100* sigma2+ 1) .*sigma2 ; % Hyperparameter a.w. Sigma2
aHy     = (100* sigma2+ 2);           % Hyperparameter a.w. Sigma2
s       = var(10^(-2)*rand(D,N1))';
w       = ones(N,1);
sbi_j   = reshape(s,row+2,col+2);
s       = reshape(sbi_j(2:end-1,2:end-1),N,1);
wbi_j   = reshape(w,row,col);
wbi_j1  = circshift(wbi_j, [0 -1]);    %LxN  gauche i j+1
wbi1_j1 = circshift(wbi_j, [-1 -1]);   %LxN  haut gauche i+1 j+1
wbi1_j  = circshift(wbi_j, [-1 0]);    %LxN  haut i+1 j
sigCEV  = 0.01; %10000000;% %0.04;
C_EV    = ones(N,1);  % Illumination coefficients
alpha   = sunsal(MPlus,Y_bloc ,'lambda',0,'ADDONE','no','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
if(sum( sum(alpha,1)>2) > L/20)
    sigCEV      = 0.0000001; % Hyperparameter a.w. C_EV
else
    sigCEV      = 0.01;    % Hyperparameter a.w. C_EV
end

alpha   = sunsal(MPlus,Y_bloc ,'lambda',0,'ADDONE','yes','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
m_NL    = (M_NL*gam_NL).*(C_EV(:,1).^2*ones(1,L))';  % LxN
if(Init ==0)
    alpha   = 1/R*ones(R,N);
    sigma2  = 10^(-2)*ones(L,1);
end

%% iterative algorithm
info      = zeros(4,1);
condition = 0;compt=1;
Cost_fun  = 0;
while condition == 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Abundance (alpha of size R x N)
    M2        = MPlus./repmat(sqrt(sigma2),1,R);
    Y_Lin2    = (Y_bloc-m_NL)./repmat(sqrt(sigma2),1,N);
    alphaIni  = alpha;
    lambda    = 0;
    alpha     = sunsal(M2,Y_Lin2./repmat(C_EV(:,compt)',L,1),'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
        'AL_iters',400,'TOL', 1e-6, 'verbose','no');
    m_Lin     = MPlus*alpha.*(C_EV(:,compt)*ones(1,L))'; % LxN
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Nonlinear coefficients (D x N)
    M_NL2     = M_NL./repmat(sqrt(sigma2),1,D);
    Y_NL2     = (Y_bloc-m_Lin)./repmat(sqrt(sigma2),1,N);
    [EVe,EVa] = svd(M_NL2'*M_NL2+small*eye(D));
    EVa       = diag(EVa);
    gam_NLIni = gam_NL;
    gam_NL    = sunsal_gam_c2(Y_NL2,M_NL2,EVe,EVa,s(:,compt),R,Posi,C_EV(:,compt));
    m_NL      = (M_NL*gam_NL).*(C_EV(:,compt).^2*ones(1,L))';  % LxN
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Noise variance (sigma of size Lx1)
    sigma2 = (sum( (Y_bloc - m_NL - m_Lin).^2,2) + 2*bHy)./(N+2 +2*aHy);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Nonlinear energies (s of size Nx1)
    bw_IG  = wbi_j + wbi_j1 + wbi1_j1 + wbi1_j; % row x col
    bw_IG  = a*bw_IG; % row x col
    gmat   = reshape(0.5* sum((gam_NL ).^2,1)',row,col); % row x col
    sbi_j(2:end-1,2:end-1) = (bw_IG + gmat    ) / (4*a+D/2+1);  % row x col
    
    sbi_j1  = circshift(sbi_j, [0 1]);    %LxN1  i j+1
    sbi1_j1 = circshift(sbi_j, [1 1]);    %LxN1  i+1 j+1
    sbi1_j  = circshift(sbi_j, [1 0]);    %LxN1  i+1 j
    s(:,compt+1)    = reshape(sbi_j(2:end-1,2:end-1),N,1);
    %     s(:,compt+1)    = max(s(:,compt+1) ,10^(-100));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Auxilliary variables (w of size Nx1)
    bs_IG        = 1./sbi_j + 1./sbi_j1 + 1./sbi1_j1 + 1./sbi1_j; % row+2 x col+2
    wbi_j        = (4*a-1)./(a* bs_IG(2:end-1,2:end-1)) ;   % row x col
    wbi_j1       = circshift(wbi_j, [0 -1]);    %LxN    i j+1
    wbi1_j1      = circshift(wbi_j, [-1 -1]);   %LxN    i+1 j+1
    wbi1_j       = circshift(wbi_j, [-1 0]);    %LxN    i+1 j
    w(:,compt+1) = reshape(wbi_j,N,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Illumination coefficients (C_EV of size Nx1)
    
    if(Illum==0)
        C_EV(:,compt+1)  = ones(N,1);
    else
        m_Ling   = MPlus*alpha; % LxN
        m_NLg    = (M_NL*gam_NL);  % LxN
        m_NL2    = m_NLg./repmat(sqrt(2*sigma2),1,N);
        m_Lin2   = m_Ling./repmat(sqrt(2*sigma2),1,N);
        Y_bloc2  = Y_bloc./repmat(sqrt(2*sigma2),1,N);
        x1       = sum(m_NL2.^2,1)'; % N x 1
        x2       = 2*sum(m_Lin2.*m_NL2 ,1)'; % N x 1
        x3       = sum(m_Lin2.^2 - 2*m_NL2.*Y_bloc2,1)' +0.5/sigCEV; % N x 1
        x4       = -2*sum(m_Lin2.*Y_bloc2,1)' - 1/sigCEV; % N x 1
        x5       = sum(Y_bloc2.*Y_bloc2,1)' + 0.5/sigCEV; % N x 1
        vectP3   = [4*x1,3*x2,2*x3,x4]; % N x 4
        for n=1:N
            P       = roots(vectP3(n,:)); % 3xN
            Der2P   = 12*x1(n)*P.^2 + 6*x2(n)*P+2*x3(n);
            Pol2P   = x1(n)*P.^4 + x2(n)*P.^3+ x3(n)*P.^2+ x4(n)*P+x5(n);
            ind     = find(Der2P>0 & real(P)>0.5 & real(P)<3);
            if(length(ind)==1)
                C_EV(n,compt+1)  =  min(real(P(ind)),3);
                C_EV(n,compt+1)  =  max(C_EV(n,compt+1),0.5);
            elseif(length(ind)==2)
                [v ind2]         = min(real(Pol2P(ind)));
                C_EV(n,compt+1)  = real(P(ind(ind2)));
                C_EV(n,compt+1)  = max(C_EV(n,compt+1),0.5);
            else
                C_EV(n,compt+1)  = C_EV(n,compt);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Test of convergence
    m_Lin  = MPlus*alpha.*(C_EV(:,compt+1)*ones(1,L))'; % LxN
    bw_IG  = wbi_j + wbi_j1 + wbi1_j1 + wbi1_j; %r x c
    bw_IG  = a*bw_IG; %r x c
    
    if mod(compt,10) == 1
        indCov              = (compt-1)/10+1;
        alphat(:,:,indCov)  = alpha;
        sigma2t(:,indCov)   = sigma2;
        gam_NLt(:,:,indCov) = gam_NL;
        Cost_fun(indCov)    =  N/2*log(sigCEV)  -4*a*sum(log(bw_IG(:))) +  L/2*log(2*pi)+ sum((N/2+1+aHy).*(log(sigma2)-log(eps))) + sum((0.5*sum((Y_bloc - m_NL - m_Lin).^2,2)+bHy)./sigma2) ...
            + sum((C_EV(:,compt+1)-1).^2)/2/sigCEV +  sum((4*a+1+D/2).* (log(s(:,compt+1))+100) ) ...
            + sum( ( sum(0.5*gam_NL.^2 ,1)' +  bw_IG(:))./s(:,compt+1) ) ;
        nh(indCov)          = norm(alpha-alphaIni,'fro');
        maxDiff(indCov)     = max(abs(alpha(:)-alphaIni(:)));
        nx(indCov)          = StpX + norm(alphaIni,'fro');
        nhG(indCov)         = norm(gam_NL-gam_NLIni,'fro');
        maxDiffG(indCov)    = max(abs(gam_NL(:)-gam_NLIni(:)));
        nxG(indCov)         = StpX + norm(gam_NLIni,'fro');
        if(compt>10)
            ng  = abs(Cost_fun(indCov)-Cost_fun(indCov-1));
            if(ng <= Stpfunc*abs(Cost_fun(indCov-1)))  info(1)=1; end
        end
        if(nh(indCov)  <= StpX*nx(indCov))     info(2)=1; end
        if(nhG(indCov) <= StpXG*nxG(indCov))   info(3)=1; end
        if(compt>StpKmax)    info(4)=1; end
        
        if(sum(info)>0)
            condition =1;
        end
        
        fprintf(' iter = %f \n',compt)
    end
    compt = compt+1;
    
end

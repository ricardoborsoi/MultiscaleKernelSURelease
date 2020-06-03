clc
clear all
global betag

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%        Constants       %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
Synth_Im  = 1;          %  1: LMM  2: LMM+NL   3: ME+EV 
noise     = 0;          %  For synthetic data
                        %  0: iid noise     1: Gauss. shape of the variance
eta       = 50;         %  Gauss. Width of the variance
                        %  eta > L/6,  L being the number of bands
SNR       = 25;         %  Noise level
R         = 3;          %  Number of endmembers (R in {3,6}) 
N         = 50;         %  Number of pixels for rows and columns
                        %  available: 50x50  or 100x100 pixels
Trunc     = 0.9999;     %  Maximum value for the abundances
D         = R*(R+1)/2;  %  Number of bilinear terms
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Synthetic/Real Data   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Generation_Synth_Data

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%  Unmixing algorithms   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%        Sunsal      %%%%%%%%%%%%%%%
disp('Sunsal')
algo = 1;
tic
lambda = 0;
alpha_sunsal  = sunsal(MPlus,Y_bloc,'lambda',lambda,'ADDONE','yes','POSITIVITY','yes', ...
    'AL_iters',200,'TOL', 1e-4, 'verbose','no');
time(algo) = toc;
[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:) REn(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_sunsal,alpha0,ones(row*col,1),c_Ev0,zeros(L,row*col),RComp,K,Label,1);


% % %%%%%%%%%%%%%      CDA + NL       %%%%%%%%%%%%%%%
disp('CDA-NL')
algo        = 2;
if(Synth_Im==1) 
      Illum       = 1;   %  0: Fixed c=1,  1: Variable c
else  Illum       = 0;   %  0: Fixed c=1,  1: Variable c
end 
Posi        = 1;   %  1: positive gamma, 0: unconstrained gamma
Init_NL     = 1;   %  0: Random initialization    1: good initialization
tic
[alpha_NLt sig2NLt gam_NLt sNLt wNLt C_NLt infoNL]  = Unmix_CDA_NL_TIP_v1(MPlus,Y,Posi,Init_NL,Illum);
time(algo) = toc;
Iter(algo) = size(alpha_NLt,3)-1;
alpha_NL   = alpha_NLt(:,:,Iter(algo));  sig2NL   = sig2NLt(:,Iter(algo)); gam_NL  = gam_NLt(:,:,Iter(algo)); sNL   = sNLt(:,Iter(algo));
wNL   = wNLt(:,Iter(algo));  C_NL   = C_NLt(:,Iter(algo));

[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_NL,alpha0,C_NL,c_Ev0,gam_NL,RComp,K,Label,2); %%%


% % %%%%%%%%%%%%%      CDA + ME      %%%%%%%%%%%%%%%
disp('CDA-ME')
algo        = 3;
if(Synth_Im==1) 
      Illum       = 1;   %  0: Fixed c=1,  1: Variable c
else  Illum       = 0;   %  0: Fixed c=1,  1: Variable c
end
Init_ME     = 1;    %  0: Random initialization    1: good initialization
tic
[alpha_MEt sig2MEt Out_ME sMEt wMEt C_MEt infoME]  = Unmix_CDA_ME_TIP_v1(MPlus,Y,Bands,L0,Init_ME,Illum);
time(algo) = toc;

Iter(algo) = size(alpha_MEt,3)-1;
alpha_ME   = alpha_MEt(:,:,Iter(algo));  sig2ME   = sig2MEt(:,Iter(algo)); Out_ME  = Out_ME; sME   = sMEt(:,Iter(algo));
wME   = wMEt(:,Iter(algo));  C_ME   = C_MEt(:,Iter(algo));

[RMSE(algo) RE(algo) SAM(algo) SAMn(algo,:) RCest(:,:,algo) RE_RC(algo) SAM_RC(algo) SAMn_RC(algo,:)   ...
    RMSE_Class(algo,:) RE_Class(algo,:) SAM_Class(algo,:) RE_Class_RC(algo,:) SAM_Class_RC(algo,:)]  = ...
    exploitation_TIP(Y_bloc,MPlus,alpha_ME,alpha0,C_ME,c_Ev0,Out_ME,RComp,K,Label,3); %%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%        Display results      %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CodeName = ['SUNSAL'; 'CDA-NL'; 'CDA-ME';];

disp(' ')
disp('**************************')
disp('Processing time (seconds)')
disp('**************************')
for i=1:3
disp([num2str(CodeName(i,:)) ':  '  num2str(time(i))]);
end
disp('**************************')
disp(' ')
disp('****************')
disp('RMSE (x 10^(-2))')
disp('****************')
for i=1:3
disp([num2str(CodeName(i,:)) ':  '  num2str(RMSE(i)*100)]);
end
disp('****************')
disp(' ')
disp('****************')
disp('RE (x 10^(-2))')
disp('****************')
for i=1:3
disp([num2str(CodeName(i,:)) ':  '  num2str(RE(i)*100)]);
end
disp('****************')
disp(' ')
disp('****************')
disp('SAM (x 10^(-2))')
disp('****************')
for i=1:3
disp([num2str(CodeName(i,:)) ':  '  num2str(SAM(i)*100)]);
end
disp('****************')



 
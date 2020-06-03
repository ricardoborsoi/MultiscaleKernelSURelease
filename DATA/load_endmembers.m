function [M,namesM] = load_endmembers(P)
% =========================================================================
% Load endmembers from the USGS library
% 
% P - number of endmembers to load
% M - endmember library
% =========================================================================


load USGS_1995_Library.mat
%  order bands by increasing wavelength
[dummy index] = sort(datalib(:,1));
A =  datalib(index,4:end);
names = names(4:end,:);

% prune the library 
% min angle (in degres) between any two signatures 
% the larger min_angle the easier is the sparse regression problem
min_angle = 4.44;       
[A, index] = prune_library2(A,min_angle); % 240  signature 
names = names(index',:);

% order  the columns of A by decreasing angles 
[A, index, angles] = sort_library_by_angle(A);
names = names(index',:);

% Names of the first 10 ordered materials, with 4.44 deg. prunning:
% 1 - Jarosite GDS99 K,Sy 200C
% 2 - Jarosite GDS101 Na,Sy 200
% 3 - Anorthite HS349.3B 
% 4 - Calcite WS272 
% 5 - Alunite GDS83 Na63 
% 6 - Howlite GDS155
% 7 - Corrensite CorWa-1
% 8 - Fassaite HS118.3B  
% 9 - Adularia GDS57 Orthoclase  
% 10 - Andradite NMNH113829 

% select p endmembers  from A
% angles (a_1,a_j) \sisizemeq min_angle)
supp = 2:(P+1);

% Remove an endmember
supp = supp(find(supp~=4));
supp = [supp P+2];


% % Sample endmembers at random
% supp = randsample(size(A,2), p);

M = A(:,supp);
[L,p] = size(M);  % L = number of bands; p = number of material

% get the names of the materials
namesM = names(supp);


% Reorder M with decreasing energy columns
[~,ii] = sort(diag(M'*M));
M = M(:,ii);


function [sppx_variance2,rec_error_sppx2,sppx_sizes2,eig_spread2,eig_ratio2] = compute_sppx_variance(Y2,slic_reg,size_max)
% =========================================================================
% 
% Computes the superpixel error and variance for different superpixel sizes
% 
% INPUTS: Y2 - an  nr * nc * L  hyperspectral image
%         slic_reg - superpixels regularity (usually 0.01, from 0.005 to 0.05)
% 
% OUTPUTS: rec_error_sppx - error as a function of superpixels size
% 
% =========================================================================



if nargin < 3
    size_max = 20;
end


% slic_reg 
% slic_size

L = size(Y2,3);
nr = size(Y2,1);
nc = size(Y2,2);
Y2a = Y2;

if nargin < 4
    gamma_rel = 1;
end

% % Estimate noise power
% [numRows,numCols,numSpectra] = size(Y2);
% imgVec = reshape(Y2, [numRows*numCols numSpectra])';
% SNR_est = estimate_SNR_im(imgVec,P);
% noise_energy = norm(Y2a(:),2)^2 * 10^(-SNR_est/10);



flag_selected = false;
idxo = 0;
rec_error_sppx = [];
sppx_variance  = [];
eig_spread     = [];
eig_ratio      = [];
KK             = [];

% Model_complexity = 1:0.5:20;
% Model_complexity = 2:0.5:20;
% Model_complexity = 1:1:12;

% Model_complexity = 2:1:12;
% Model_complexity = [2:1:20];
% Model_complexity = 2:1:15;
% Model_complexity = [3:1:25];
% Model_complexity = 1.5:1:40;


Model_complexity = [2:1:size_max];

for slic_size = Model_complexity
    idxo = idxo + 1;
    
    
    % reorder and rescale data into 2-D array
    [numRows,numCols,numSpectra] = size(Y2a);
    scfact = mean(reshape(sqrt(sum(Y2a.^2,3)), numRows*numCols, 1));
    Y2 = Y2a./scfact;
%     imgVec = reshape(Y2, [numRows*numCols numSpectra]);

    % compute superpixels
    disp('Computing SLIC Superpixels...');
    spSegs = vl_slic(single(Y2), slic_size, slic_reg);
    numSuperpixels = double(max(spSegs(:)))+1; 


    % Unmix the superpixels ------------------------- 
    Y3 = zeros(size(Y2));
    avg_superpx = zeros(1, numSuperpixels+1, L);
    all_superpx = cell(1,numSuperpixels);
    for i=0:numSuperpixels
        [rowi, coli] = find(spSegs==i);

        for j=1:length(rowi)
            % Averages all pixels inside each superpixel
            if j == 1
                avg_superpx(1,i+1,:) = (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
                all_superpx{i+1} = squeeze(Y2a(rowi(j),coli(j),:));
            else
                avg_superpx(1,i+1,:) = avg_superpx(1,i+1,:) + (1/length(rowi)) * Y2a(rowi(j),coli(j),:);
                all_superpx{i+1} = [all_superpx{i+1} squeeze(Y2a(rowi(j),coli(j),:))];
            end
        end

        % This is optional (for visualization)
        for j=1:length(rowi)
            Y3(rowi(j),coli(j),:) = avg_superpx(1,i+1,:);
        end
    end

    
    % compute the error ----------------------
%     rec_error_sppx(idxo) = norm(Y2(:) - Y3(:))^2;
    rec_error_sppx(idxo) = norm(1 * (Y2a(:) - Y3(:)))^2;

    
    % compute variance of each superpixel
    sppx_variance_tmp = 0;
    count_sppx = 0;
    for i=1:numSuperpixels
        if numel(all_superpx{i}) > 0
            count_sppx = count_sppx + 1;
            sppx_variance_tmp = sppx_variance_tmp + mean(var(all_superpx{i}'));
        end
    end
    sppx_variance(idxo) = sppx_variance_tmp / count_sppx;
    KK(idxo) = numSuperpixels;
    
    
    
    % compute singular value spread in each superpixel
    eig_spread_tmp = 0;
    eig_ratio_tmp  = 0;
    count_sppx = 0;
    for i=1:numSuperpixels
        if size(all_superpx{i},2) > 1 % numel(all_superpx{i}) > 0
            count_sppx = count_sppx + 1;
            Sv = svd(all_superpx{i});
            eig_spread_tmp = eig_spread_tmp + abs(Sv(1))/sum(abs(Sv));
            eig_ratio_tmp  = eig_ratio_tmp  + abs(Sv(1))/sum(abs(Sv(2)));
        end
    end
    eig_spread(idxo) = eig_spread_tmp / count_sppx;
    eig_ratio(idxo)  = eig_ratio_tmp  / count_sppx;
    KK(idxo) = numSuperpixels;
    
end

% sppx_sizes = Model_complexity.^2;
sppx_sizes = (nr*nc)./KK;

% % add superpixel size of 1:
% sppx_sizes     = [1 sppx_sizes];
% sppx_variance  = [0 sppx_variance];
% rec_error_sppx = [0 rec_error_sppx];


% clear vector of repeated values
sppx_sizes2     = unique(sppx_sizes);
sppx_variance2  = zeros(size(sppx_sizes2));
rec_error_sppx2 = zeros(size(sppx_sizes2));
eig_spread2     = zeros(size(sppx_sizes2));
eig_ratio2      = zeros(size(sppx_sizes2));

for i=1:length(sppx_sizes2)
    iidx = find(sppx_sizes == sppx_sizes2(i));
    sppx_variance2(i)  = mean(sppx_variance(iidx));
    rec_error_sppx2(i) = mean(rec_error_sppx(iidx));
    eig_spread2(i)     = mean(eig_spread(iidx));
    eig_ratio2(i)      = mean(eig_ratio(iidx));
end

% sppx_sizes2     = sppx_sizes;
% sppx_variance2  = sppx_variance;
% rec_error_sppx2 = rec_error_sppx;
% eig_spread2     = eig_spread;
% eig_ratio2      = eig_ratio;
end




function [SNR] = estimate_SNR_im(R,P)
% Estimate the SNR of a hyperspectral image
% R - L * N HSI
% P - number of endmembers

[L, N]=size(R);

verbose = 'on';

r_m = mean(R,2);      
R_m = repmat(r_m,[1 N]); % mean of each band
R_o = R - R_m;           % data with zero-mean 
[Ud,Sd,Vd] = svds(R_o*R_o'/N,P);  % computes the p-projection matrix 
x_p =  Ud' * R_o;                 % project the zero-mean data onto p-subspace

SNR = estimate_snr(R,r_m,x_p);

if strcmp (verbose, 'on'), fprintf(1,'SNR estimated = %g[dB]\n',SNR); end
end

function snr_est = estimate_snr(R,r_m,x)

 [L N]=size(R);           % L number of bands (channels)
                          % N number of pixels (Lines x Columns) 
 [p N]=size(x);           % p number of endmembers (reduced dimension)

 P_y = sum(R(:).^2)/N;
 P_x = sum(x(:).^2)/N + r_m'*r_m;
 snr_est = 10*log10( (P_x - p/L*P_y)/(P_y- P_x) );
end






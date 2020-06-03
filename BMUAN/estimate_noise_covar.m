function [SNR,Sigma_e]=estimate_noise_covar(R,P)
% Estimates the SNR and the noise covariance matrix given an HSI
%
% IN: R - HSI, with dimension L * N
%     P - Number of endmembers
%
% OUT: SNR     - Estimated signal to noise ratio
%      Sigma_e - Estimated noise covariance matrix
% 
% The noise covariance matrix is estimated using standard regression 
% residual analysis. It is described in the papers:
% 
%     R E Roger. Principal components transform with simple, automatic noise adjustment.
%     International journal of remote sensing, 1996.
% and
% 
%     Modified Residual Method for the Estimation of Noise in Hyperspectral Images
%     Asad Mahmood, Amandine Robin, and Michael Sears
% 
% Author: Ricardo Borsoi, 2018
% -----------------------------------

[L, N]=size(R);

verbose = 'on';

r_m = mean(R,2);      
R_m = repmat(r_m,[1 N]); % mean of each band
R_o = R - R_m;           % data with zero-mean 
[Ud,Sd,Vd] = svds(R_o*R_o'/N,P);  % computes the p-projection matrix 
x_p =  Ud' * R_o;                 % project the zero-mean data onto p-subspace

SNR = estimate_snr(R,r_m,x_p);


% Lets compute the residuals for a least-squares fit of the data
X = R_o';
alphai = zeros(L-1,1);
resd = zeros(N,L);
for i=1:L
    % solve Xi = X^{\{i}} alphai
    X_i = X(:,i);
    X_ni = X(:,[1:(i-1) (i+1):L]);
    alphai = X_ni \ X_i;
    % compute residual
    resd(:,i) = X_i - X_ni * alphai;
end

Sigma_e = (1/N) * (resd' * resd);



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


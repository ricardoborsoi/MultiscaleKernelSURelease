
function [sel_slic_size] = homogeneity_curvature_sppx_size(Y2,slic_reg,der_thr)
% =========================================================================
% 
% Computes the optimal superpixel size
% 
% INPUTS: Y2       : an  nr * nc * L  hyperspectral image
%         slic_reg : superpixels regularity (usually 0.01, from 0.005 to 0.05)
%         kappa    : parameter weighting model complexity vs.
%                    reconstruction error (default = 10)
% 
% 
% OUTPUTS: sel_slic_size : selected superpixels equivalent side (i.e. 
%                         the number of pixels is sel_slic_size^2)
% 
% =========================================================================

size_max = 15; 
[sppx_variance,rec_error_sppx,sppx_sizes,eig_spread2,eig_ratio2] = compute_sppx_variance(Y2,slic_reg,size_max);


% discard first sample?

sppx_sizes  = sppx_sizes(2:end);
eig_spread2 = eig_spread2(2:end);
eig_ratio2  = eig_ratio2(2:end);


% sppx_sizes
% eig_spread2
% eig_ratio2 


% ------------------------------
% first value within 15% of the maximum ratio:

ydata = eig_ratio2;
xdata = sppx_sizes;
[pp, ~] = csaps( xdata, ydata , 1.2);

xdata = [4:1:13].^2; 
% xdata = [3:1:13].^2; 
% xdata = [2:1:13].^2; 
yy = fnval(pp, xdata);

maxratio = max(yy);
%maxratio = sort(yy,2,'descend');
%maxratio = mean(maxratio(1:4))
posi = min(find( yy >= maxratio*(1-der_thr) ));

% % % % if posi == 1
% % % %     posi = min(find( yy <= maxratio*(1-der_thr) ));
% % % % end

sel_slic_size = sqrt(xdata(posi));
sel_slic_size = round(sel_slic_size);


% if false
% % ------------------------------
% % first value within 15% of the minimum spread:
% 
% ydata = eig_spread2;
% xdata = sppx_sizes;
% [pp, ~] = csaps( xdata, ydata , 1.2);
% 
% xdata = [3:1:20].^2; 
% yy = fnval(pp, xdata);
% 
% minratio = min(yy);
% posi = min(find( yy <= minratio*(1+der_thr) ));
% sel_slic_size = sqrt(xdata(posi));
% sel_slic_size = round(sel_slic_size);
% end



end


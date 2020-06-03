function IHHT=IpHHT3(w,h)
% HI = sparse(eye(N));
N = w*h;
HI = speye(N);
% H=ConvH(HI);
H=ConvH3(HI,w,h);

% HHT=ConvHT(H);
HHT=ConvHT3(H,w,h);
IHHT=HI+HHT;
%INV_IHHT=inv(IHHT);
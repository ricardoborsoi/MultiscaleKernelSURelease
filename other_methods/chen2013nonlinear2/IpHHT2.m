function IHHT=IpHHT2(N);
% HI = sparse(eye(N));
HI = speye(N);
H=ConvH(HI);

HHT=ConvHT(H);
IHHT=HI+HHT;
%INV_IHHT=inv(IHHT);
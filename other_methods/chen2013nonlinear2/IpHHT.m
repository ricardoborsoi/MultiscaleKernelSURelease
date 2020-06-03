function IHHT=IpHHT(N);
HI = sparse(eye(N*N));
H=ConvH(HI);

HHT=ConvHT(H);
IHHT=HI+HHT;
%INV_IHHT=inv(IHHT);
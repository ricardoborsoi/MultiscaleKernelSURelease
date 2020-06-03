function U0 = ConvH(V,w,h)

if nargin == 1
     % h=sqrt(size(V,2));
     h=floor(sqrt(size(V,2)));
end

N = size(V,2);


 U1= V(:,1:end-1)-V(:,2:end);
 U1 = [U1,V(:,end)];
 pos = h:h:N;
 U1(:,pos)=0;
 
 
 U2= V(:,2:end)-V(:,1:end-1);
 U2 = [V(:,1),U2];
 pos = 1:h:N;
 U2(:,pos)=0;
 
 U3= V(:,1:end-h)-V(:,h+1:end);
 U3 = [U3,V(:,end-h+1:end)];
 U3(:,end-h+1:end)=0;
 
 U4= V(:,h+1:end)-V(:,1:end-h);
 U4 = [V(:,1:h), U4];
 U4(:,1:h) = 0;
 
 U0 = [U1,U2,U3,U4];
 
 
 
 
 
 
%  N=sqrt(size(V,2));
% 
% 
% 
%  U1= V(:,1:end-1)-V(:,2:end);
%  U1 = [U1,V(:,end)];
%  
%  U2= V(:,2:end)-V(:,1:end-1);
%  U2 = [V(:,1),U2];
%  
%  U3= V(:,1:end-N)-V(:,N+1:end);
%  U3 = [U3,V(:,end-N+1:end)];
%  
%  U4= V(:,N+1:end)-V(:,1:end-N);
%  U4 = [V(:,1:N), U4];
%  
%  U0 = [U1,U2,U3,U4];
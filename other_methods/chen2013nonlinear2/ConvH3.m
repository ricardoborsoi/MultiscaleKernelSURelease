function U0 = ConvH3(V,w,h)

V = full(V);

Vim = reshape(V', h,w,[]);
P = size(Vim,3);

ff1 = [0 0 0; 0 1 -1; 0 0 0];
ff2 = [0 0 0; -1 1 0; 0 0 0];
ff3 = [0 0 0; 0 1 -1; 0 0 0]';
ff4 = [0 0 0; -1 1 0; 0 0 0]';

temp1 = zeros(size(Vim));
temp2 = zeros(size(Vim));
temp3 = zeros(size(Vim));
temp4 = zeros(size(Vim));

for j=1:P
    temp1(:,:,j) = imfilter(Vim(:,:,j), ff1, 'circular', 'same', 'conv');
    temp2(:,:,j) = imfilter(Vim(:,:,j), ff2, 'circular', 'same', 'conv');
    temp3(:,:,j) = imfilter(Vim(:,:,j), ff3, 'circular', 'same', 'conv');
    temp4(:,:,j) = imfilter(Vim(:,:,j), ff4, 'circular', 'same', 'conv');
end

U1 = reshape(temp1,h*w,P)';
U2 = reshape(temp2,h*w,P)';
U3 = reshape(temp3,h*w,P)';
U4 = reshape(temp4,h*w,P)';

U0 = [U1,U2,U3,U4];




% reshape(V',)

% 
% if nargin == 1
%      % h=sqrt(size(V,2));
%      h=floor(sqrt(size(V,2)));
% end
% 
% N = size(V,2);
% 
% 
%  U1= V(:,1:end-1)-V(:,2:end);
%  U1 = [U1,V(:,end)];
%  pos = h:h:N;
%  U1(:,pos)=0;
%  
%  
%  U2= V(:,2:end)-V(:,1:end-1);
%  U2 = [V(:,1),U2];
%  pos = 1:h:N;
%  U2(:,pos)=0;
%  
%  U3= V(:,1:end-h)-V(:,h+1:end);
%  U3 = [U3,V(:,end-h+1:end)];
%  U3(:,end-h+1:end)=0;
%  
%  U4= V(:,h+1:end)-V(:,1:end-h);
%  U4 = [V(:,1:h), U4];
%  U4(:,1:h) = 0;
%  
%  U0 = [U1,U2,U3,U4];
 
 
 
 
 
 
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
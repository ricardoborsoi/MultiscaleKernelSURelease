function T = ConvHT3(X,w,h)

X = full(X);

P = size(X,1);
N = size(X,2)/4;

X1=X(:,1:N);
X2=X(:,N+1:2*N);
X3=X(:,2*N+1:3*N);
X4=X(:,3*N+1:end);

X1im = reshape(X1', h,w,[]);
X2im = reshape(X2', h,w,[]);
X3im = reshape(X3', h,w,[]);
X4im = reshape(X4', h,w,[]);



ff1 = [0 0 0; 0 1 -1; 0 0 0];
ff2 = [0 0 0; -1 1 0; 0 0 0];
ff3 = [0 0 0; 0 1 -1; 0 0 0]';
ff4 = [0 0 0; -1 1 0; 0 0 0]';

temp1 = zeros(size(X1im));
temp2 = zeros(size(X2im));
temp3 = zeros(size(X3im));
temp4 = zeros(size(X4im));

for j=1:P
    temp1(:,:,j) = imfilter(X1im(end:-1:1,end:-1:1,j), ff1, 'circular', 'same', 'conv');
    temp2(:,:,j) = imfilter(X2im(end:-1:1,end:-1:1,j), ff2, 'circular', 'same', 'conv');
    temp3(:,:,j) = imfilter(X3im(end:-1:1,end:-1:1,j), ff3, 'circular', 'same', 'conv');
    temp4(:,:,j) = imfilter(X4im(end:-1:1,end:-1:1,j), ff4, 'circular', 'same', 'conv');
end
temp1 = temp1(end:-1:1,end:-1:1,:);
temp2 = temp2(end:-1:1,end:-1:1,:);
temp3 = temp3(end:-1:1,end:-1:1,:);
temp4 = temp4(end:-1:1,end:-1:1,:);


U1 = reshape(temp1,h*w,P)';
U2 = reshape(temp2,h*w,P)';
U3 = reshape(temp3,h*w,P)';
U4 = reshape(temp4,h*w,P)';

T = U1 + U2 + U3 + U4;




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
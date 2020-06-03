function T = ConvHT(X,w,h)

N = size(X,2)/4;

if nargin == 1
     % h=sqrt(N);
     h=floor(sqrt(N));
end

X1=X(:,1:N);
X2=X(:,N+1:2*N);
X3=X(:,2*N+1:3*N);
X4=X(:,3*N+1:end);



pos = h:h:N;
XT = X1;
XT(:, pos)=0;

U1= XT(:,2:end)-XT(:,1:end-1);
U1 = [X1(:,1),U1];



pos = 1:h:N;
XT = X2;
XT(:, pos)=0;

U2= XT(:,1:end-1)-XT(:,2:end);
U2 = [U2,X2(:,end)];


XT = X3;
XT(:,end-h+1:end)=0;
U3= XT(:,h+1:end)-XT(:,1:end-h);
U3 = [X3(:,1:h),U3];

XT = X4;
XT(:,1:h)=0;
U4= XT(:,1:end-h)-XT(:,h+1:end);
U4 = [U4,X4(:,end-h+1:end)];



T = U1+U2+U3+U4;




% N = size(X,2)/4;
% N1 = sqrt(N);
% 
% X1=X(:,1:N);
% X2=X(:,N+1:2*N);
% X3=X(:,2*N+1:3*N);
% X4=X(:,3*N+1:end);
% 
% U1= X2(:,1:end-1)-X2(:,2:end);
% U1 = [U1,X2(:,end)];
% 
% U2= X1(:,2:end)-X1(:,1:end-1);
% U2 = [X1(:,1),U2];
% 
% 
%  U3= X4(:,1:end-N1)-X4(:,N1+1:end);
%  U3 = [U3,X4(:,end-N1+1:end)];
%  
%  U4= X3(:,N1+1:end)-X3(:,1:end-N1);
%  U4 = [X3(:,1:N1), U4];
% 
% 
% T = U1+U2+U3+U4;

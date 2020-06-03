function [X_NDU,F_NDU] = nonlinearU_vectValKernels(Yim,R,lambda,str,FCLS_init)
% R % endmembers


[nl,nc,L] = size(Yim);
nt = nc*nl; N = nt;
% S = reshape(Yim,nt,L)';

P = size(R,2);

% FCLS_init_im = reshape(FCLS_init,nl,nc,P);

% some parameters
% % lambda = 10; % nonlinearity regularization (can tune)
mu = 0.0001; % abundance norm regularization (can keep small)

% str = 'Transformable'; % 'full'; %
% % str = 'Separable';
par = 0; % Polynomial if 0, Gaussian otherwise (can keep Polynomial, worked well in general)
% par = 0;

if nargin < 3
    lambda = 10;
    str = 'Transformable';
end

%% =============== NDU (P) =====================

X3d = zeros(nl,nc,P); 
F3d = zeros(nl,nc,L);
t = 0;

% h = 10;
% w = 10;
% % patch 1x100
h = 1;
w = 1;

totalIters = nl*nc;
countPercent = 1;
fprintf('\n Vector valued kernels... \n')

for i=1:floor(nl/h)
    for j=1:floor(nc/w)
        
        currentIter = (i-1)*nc + j;
        if 100*currentIter/totalIters > 10*countPercent
            fprintf('%d percent...',10*countPercent)
            countPercent = countPercent + 1;
        end
        
        % select block of h*w
        xi = (i-1)*h+1:i*h; 
        yj = (j-1)*w+1:j*w;
        S3d = Yim(xi,yj,:);
        Ni = h*w; 
        Si = reshape(S3d,Ni,L).';
        
        % Setting the neighbours of each pixel 
        V = cell(1,Ni);
        [S1, S2, S3, S4] = FourNeighbours(Yim,h,w,i,j); % returns 4 neighbours of the (i,j)-th patch
        for ii=1:Ni
            V{1,ii} = [Si(:,ii) S1(:,ii) S2(:,ii) S3(:,ii) S4(:,ii)];
        end
        
        [X_FCLSi, t_FCLS] = FCLS(Si,R);
% % % %         X_FCLSi = FCLS_init_im(xi,yj,:);
        
        
        abd_map2_i = reshape(X_FCLSi,Ni,P)';
        abd_map2_i = reshape(abd_map2_i,Ni,P)';
        
        % run the algorithm
        [Xi, Fi, t_NDU, K_NDU] = NDU(Si,R,V,lambda,mu,par,str);
                                    
        X3di = reshape(Xi',h,w,P);
        F3di = reshape(Fi',h,w,L);
        X3d(xi,yj,:)=X3di;
        F3d(xi,yj,:)=F3di;
        t = t+ t_NDU; 
%         disp(strcat(num2str(i),'-',num2str(j)))

    end
end

fprintf('\n Done! \n')


X_NDU = reshape(X3d,N,P)';
F_NDU = reshape(F3d,N,L)';

% RE_NDU = ErrComput(S,R*X_NDU+F_NDU);
% Angle_NDU = Compute_avg_angle(S,R*X_NDU+F_NDU);
end






% =========================================================================
function [S1, S2, S3, S4] = FourNeighbours(M3d,h,w,i,j)
[nl,nc,L] = size(M3d);

S1 = zeros(L,h*w); 
S2 = zeros(L,h*w); 
S3 = zeros(L,h*w); 
S4 = zeros(L,h*w); 

xi = (i-1)*h+1:i*h; 
yj = (j-1)*w+1:j*w;

cnt = 0;
for yp = yj
    for xp = xi
        cnt = cnt + 1;
        
        % left neighbour
        xl = xp;
        yl = max(yp-1,1);
        % right neighbour
        xr = xp;
        yr = min(yp+1,nc);
        % upper
        xu = max(xp-1,1);
        yu = yp;
        % lower
        xlo = min(xp+1,nl);
        ylo = yp;
        
        S1(:,cnt) = squeeze(M3d(xl,yl,:));
        S2(:,cnt) = squeeze(M3d(xr,yr,:));
        S3(:,cnt) = squeeze(M3d(xu,yu,:));
        S4(:,cnt) = squeeze(M3d(xlo,ylo,:));
    end
end
end





% =========================================================================
function [X_FCLS, t_FCLS] = FCLS(S,R)
H = R'*R; 
f = -R'*S;
P = size(R,2);
N = size(S,2);
l = zeros(P,1); 
A = ones(1,P);
b = 1;
X_FCLS = zeros(P,N);
t = clock;
for i=1:N
    % [x,err,lm] = qp(H,f,L,k,A,b,l,u,display);
    X_FCLS(:,i) = qpas(H,f(:,i),[],[],A,b,l);
end
t_FCLS = etime(clock,t);
end














% =========================================================================










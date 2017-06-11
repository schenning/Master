function w = rmusic(x,n,m,fs)
% Root-MUSIC
% Written by Henning Schei 


%      x  ->  the data vector
%      n  ->  the model order
%      m  ->  the order of the covariance matrix 
%      w  <-  the frequency estimates

% variable check
    narginchk(3,4)
    if nargin == 3
        fs = 1;
    end

    
    R = compute_autocovariance(x,m);
% singular value decomposition 
    [U,S,V] = svd(R);
    G = U(:,n+1:end);
    P = G * conj(G)';

% Forward backward approach?

% Calculate sum of each diagonal in P
    [r ,c]=size(P);
    idx=bsxfun(@plus,(r:-1:1)',0:c-1);
    Q=flipud(accumarray(idx(:),P(:)));

% Find the roots inside the unit circle and imag-val != 0    
    rts = roots(Q);
    rts = rts(abs(rts) < 1);
    rts = rts(imag(rts) ~=0);
% Find the n roots closest to the unit circle
    dfc = abs(abs(rts)-1);
    idx_sort= argsort(dfc);
    component_roots = rts(idx_sort(1:n));
    

    ang = -angle(component_roots);
    w = fs*ang;
%%
function R = compute_autocovariance(x,M)
    % Henning Schei
    % Function to compute the corvarinace matrix of a signal. 
    
    
    % x --> array of input signal of size N
    % M --> int, opitional. Size of signal block
    % covmtx <-- NxN auto-covariance matrix
    
    N = length(x); 
    
    x_vect = x(:); % Force x to be a column vector
        
    % Initial covariance matrix
    yn = flipud(x_vect(1:M));
    R = yn * conj(yn)';
    for i = 2:N-M
        yn = x_vect(M-1+i:-1:i);
        R = R + yn * conj(yn)';
    end
    
    R = R/N;

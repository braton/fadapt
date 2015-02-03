function [w,y] = fcDHTBpNLMS( u,d,m,b,B,N )
% fcDHTBpNLMS - Constrained DHT Block Power Normalized Least Mean Square
%               with error estimate by convolution in fullband

% [w,y] = fcDHTBpNLMS( u,d,m,b,B,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% b         - smoothing factor;
% B         - block length;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Mar-2006$
% $Revision: 03-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 648
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
K = 3*B;
uLen = length(u);
uFb = zeros(1,N);
y = zeros(1,uLen);
eB = zeros(B,1);
uB = zeros(2*B,1);
uBH = zeros(K,ceil(N/B));
p = zeros(K,1);
lkc = zeros(K,ceil(N/B));
Ip = eye(B);
Ip = Ip( :,B:-1:1 );
% DHT Matrix
H=zeros(K);
for i=1:K
    for j=1:K
        H(i,j)=(cos(2*pi*(i-1)*(j-1)/K)-sin(2*pi*(i-1)*(j-1)/K))/sqrt(K);
    end;
end;
% Filtering ---------------------------------------------------------------
for i=1:ceil(uLen/B)
    for j=1:B
      uFb = [ (((i-1)*B+j)<=uLen)*u(min((i-1)*B+j,uLen)),uFb(1:N-1) ];    
      y(min((i-1)*B+j,uLen)) = uFb*w;
      eB(B-j+1) = d(min((i-1)*B+j,uLen))-y(min((i-1)*B+j,uLen)); 
    end;
    uB = [ [zeros(1,i*B-uLen),u(min(i*B,uLen):-1:(i-1)*B+1)].';uB(1:B) ];
    uBH = [ H*([0;uB(1:2*B-1);zeros(B,1)]),uBH(:,1:length(uBH(1,:))-1) ];
    ek = H*([eB;zeros(K-B,1)]);
    p = b*p + (1-b)*abs(uBH(:,1)).^2;
    for j=1:K
        lkc(j,:) = lkc(j,:) + m*(ek(j)/p(j))*conj(uBH(j,:));
    end;
    % Subband/fullband mapping
    g = H*(lkc);
    g = [zeros(B,1),eye(B,B),zeros(B,B-1),Ip]*g;
    for j = 1:ceil(N/B)
        w((j-1)*B+1:j*B) = g(:,j); 
    end ;
    w = (1/2)*w;
end;
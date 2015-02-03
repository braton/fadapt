function [w,y] = fcDCTBpNLMS( u,d,m,b,B,N )
% fcDCTBpNLMS - Constrained DCT Block Power Normalized Least Mean Square
%               with error estimate by convolution in fullband
%
% [w,y] = fcDCTBpNLMS( u,d,m,b,B,N )
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
% p. 640
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
if mod(B,2)
    K = 3.5*B-1.5;
    a = 1.5*B-0.5;
    n = 0.5*B+0.5;
else
    K = 3.5*B-2;
    a = 1.5*B-1;
    n = 0.5*B;
end;
uLen = length(u);
uFb = zeros(1,N);
y = zeros(1,uLen);
eB = zeros(B,1);
uB = zeros(2*B,1);
uBC = zeros(K,ceil(N/B));
p = zeros(K,1);
lkc = zeros(K,ceil(N/B));
% Filtering ---------------------------------------------------------------
for i=1:ceil(uLen/B)
    for j=1:B
      uFb = [ (((i-1)*B+j)<=uLen)*u(min((i-1)*B+j,uLen)),uFb(1:N-1) ];    
      y(min((i-1)*B+j,uLen)) = uFb*w;
      eB(B-j+1) = d(min((i-1)*B+j,uLen))-y(min((i-1)*B+j,uLen)); 
    end;
    uB = [ [zeros(1,i*B-uLen),u(min(i*B,uLen):-1:(i-1)*B+1)].';uB(1:B) ];
    uBC = [ idct([zeros(a,1);uB(1:2*B-1)]),uBC(:,1:length(uBC(1,:))-1) ];
    ek = idct([zeros(n,1);eB;zeros(2*B-2,1)]);
    p = b*p + (1-b)*abs(uBC(:,1)).^2;
    for j=1:K
        lkc(j,:) = lkc(j,:) + m*(ek(j)/p(j))*conj(uBC(j,:));
    end;
    % Subband/fullband mapping
    g = dct(lkc);
    g = g(B:2*B-1,:);
    for j=1:ceil(N/B)
        w((j-1)*B+1:j*B) = (1/sqrt(2*K))*g(:,j); 
    end;
end;
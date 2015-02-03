function [w,y] = fcDFTBpNLMS( u,d,m,b,B,N )
% fcDFTBpNLMS - Constrained DFT Block Power Normalized Least Mean Square
%               with error estimate by convolution in fullband
%
% [w,y] = fcDFTBpNLMS( u,d,m,b,B,N )
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
% p. 628
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
uFb = zeros(1,N);
y = zeros(1,uLen);
eB = zeros(B,1);
uB = zeros(2*B,1);
uBF = zeros(2*B,ceil(N/B));
p = zeros(2*B,1);
lkc = zeros(2*B,ceil(N/B));
fn = 1/sqrt(2*B);
% Filtering ---------------------------------------------------------------
for n=1:ceil(uLen/B)
    for i=1:B
      uFb = [ (((n-1)*B+i)<=uLen)*u(min((n-1)*B+i,uLen)),uFb(1:N-1) ];    
      y(min((n-1)*B+i,uLen)) = uFb*w;
      eB(B-i+1) = d(min((n-1)*B+i,uLen))-y(min((n-1)*B+i,uLen)); 
    end;
    uB = [ [zeros(1,n*B-uLen),u(min(n*B,uLen):-1:(n-1)*B+1)].';uB(1:B) ];
    uBF = [ fn*fft(uB),uBF(:,1:length(uBF(1,:))-1) ];
    ek = fn*fft( [eB;zeros(B,1)] );
    p = b*p + (1-b)*abs(uBF(:,1)).^2;
    for i=1:2*B
        lkc(i,:) = lkc(i,:) + m*(ek(i)/p(i))*conj(uBF(i,:));
    end;
    % Subband/fullband mapping
    g = fn*fft(lkc);
    g = g(1:B,:);
    for i=1:ceil(N/B)
        w((i-1)*B+1:i*B) = g(:,i); 
    end;
end;
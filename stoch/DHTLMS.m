function [w,y] = DHTLMS( u,d,m,b,N )
% DHTLMS - Discrete Hartley Transform Least Mean Square
%
% [w,y] = DHTLMS( u,d,m,b,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% b         - smoothing factor;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Mar-2006$
% $Revision: 02-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 573
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
p = zeros(1,N);
% DHT Matrix
H=zeros(N);
for i=1:N
    for j=1:N
        H(i,j)=(cos(2*pi*(i-1)*(j-1)/N)-sin(2*pi*(i-1)*(j-1)/N))/sqrt(N);
    end;
end;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    uiH = H*ui.';
    p = b*p + (1-b)*abs(uiH.').^2;
    D = diag(p);
    y(i) = uiH.'*w;
    w = w + m*inv(D)*conj(uiH)*(d(i)-y(i));
end;
w = H*w;
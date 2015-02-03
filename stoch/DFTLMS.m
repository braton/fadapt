function [w,y] = DFTLMS( u,d,m,b,N )
% DFTLMS - Discrete Fourier Transform Least Mean Square
%
% [w,y] = DFTLMS( u,d,m,b,N )
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
% p. 581
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
% DFT Matrix
F = zeros(N);
for k=1:N
    for n=1:N
        F(k,n) = exp( -j*2*pi*((k-1)/N)*(n-1) )/sqrt(N);
    end
end
S = diag(F(2,:)*sqrt(N));
p = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = ui*S + ((u(i)-u(max(i-N,1))*(i>N))/sqrt(N))*ones(1,N);
    p = b*p + (1-b)*abs(ui).^2;
    D = diag(p);
    y(i) = ui*w;
    w = w + m*inv(D)*ui'*(d(i)-y(i));
end;
w = F*w;
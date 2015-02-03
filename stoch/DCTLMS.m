function [w,y] = DCTLMS( u,d,m,b,N )
% DCTLMS - Discrete Cosine Transform Least Mean Square
%
% [w,y] = DCTLMS( u,d,m,b,N )
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
% p. 582
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
alfa = (((0:N-1)==0)/sqrt(N)+((0:N-1)~=0)*sqrt(2/N));
% DCT Matrix
C = zeros(N);
for i=1:N
    for j=1:N
        C(i,j) = cos( (i-1)*pi*(2*j-1)/(2*N) );
    end
end
C(2:N,:) = C(2:N,:)*sqrt(2/N);
C(1,:) = C(1,:)*(1/sqrt(2));
S = diag( 2*cos((0:N-1)*pi/N) );
p = zeros(1,N);
uip = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    a = (u(i)-u(max(i-1,1))*(i>1))*cos((0:N-1)*pi/(2*N));
    z = (u(max(i-N,1))*(i>N)-u(max(i-N-1,1))*(i>N+1))*((-1).^(0:N-1)).* ...
    cos((0:N-1)*pi/(2*N));
    h = (a-z).*alfa-uip;
    uip = ui;
    ui = ui*S + h;
    p = b*p + (1-b)*abs(ui).^2;
    D = diag(p);
    y(i) = ui*w;
    w = w + m*inv(D)*ui'*(d(i)-y(i));
end;
w = C.'*w;
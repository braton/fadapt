function [w,y] = epNLMS( u,d,m,e,b,N )
% epNLMS - e-Power Normalized Least Mean Square
%
% [w,y] = epNLMS( u,d,m,e,b,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% e         - regularization factor;
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
% p. 227
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen=length(u);
y=zeros(1,uLen);p=0;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i:-1:max(i-N+1,1)),zeros(1,N-i) ];
    p = b*p + (1-b)*abs(ui(1)).^2;
    y(i) = ui*w;
    w = w+diag((m/N)./((e/N)+p))*ui'*(d(i)-y(i));
end;
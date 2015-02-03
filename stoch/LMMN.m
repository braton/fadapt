function [w,y] = LMMN( u,d,m,b,N )
% LMMN - Least Mean Mixed Norm
%
% [w,y] = LMMN( u,d,m,b,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% b         - mixing factor;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 07-Mar-2006$
% $Revision: 02-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 233
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i:-1:max(i-N+1,1)),zeros(1,N-i) ];
    y(i) = ui*w;
    w = w + m*ui'*(d(i)-y(i))*(b+(1-b)*abs(d(i)-y(i)).^2);
end;
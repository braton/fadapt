function [w,y] = LMS( u,d,m,N )
% LMS - Least Mean Square
%
% [w,y] = LMS( u,d,m,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
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
% p. 214
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w = zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    w = w+m*ui'*(d(i)-y(i));
end;
function [w,y] = RLS( u,d,a,e,N )
% RLS - Recursive Least Squares
%
% [w,y] = RLS( u,d,a,e,N )
%
% u         - Input signal;
% d         - desired signal;
% a         - forgetting factor;
% e         - regularization factor;
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
% p. 733
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
P = (1/e)*eye(N);
y = zeros(1,uLen);
ui = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    gamma = 1/(1+(1/a)*ui*P*ui');
    g = (1/a)*P*ui'*gamma;
    w = w + g*(d(i)-y(i));
    P = (P/a) - (g*g')/gamma;
end
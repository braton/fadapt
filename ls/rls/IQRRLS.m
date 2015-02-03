function [w,y] = IQRRLS( u,d,a,e,N )
% IQRRLS - Inverse QR-Recursive Least Squares
%
% [w,y] = IQRRLS( u,d,a,e,N )
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
% p. 785
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
P = sqrt(1/e)*eye(N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    [Q,R]= qr( ([1,sqrt(1/a)*ui*P;zeros(N,1),sqrt(1/a)*P])' );
    w = w + (R(1,2:N+1)'/R(1,1))*(d(i)-y(i));
    P = R(2:N+1,2:N+1)';
end
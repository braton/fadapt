function [w,y] = eAPA( u,d,m,e,p,N )
% eAPA - e-Affine Projection Algorithm
%
% [w,y] = eAPA( u,d,m,e,p,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% e         - regularization factor;
% p         - projection order;
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
% p. 238
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
Ui = zeros(p,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    Ui(2:p,:) = Ui(1:p-1,:);
    Ui(1,:) = [ u(i:-1:max(i-N+1,1)),zeros(1,N-i) ];
    di = transpose([ d(i:-1:max(i-p+1,1)),zeros(1,p-i) ]);
    Yi = Ui*w;
    y(i) = Yi(1);
    w = w + m*Ui'*inv(e*eye(p)+Ui*Ui')*(di-Yi);
end;
function [w,y] = ePRA( u,d,m,e,p,N )
% ePRA - e-Partial Rank Algorithm
%
% [w,y] = ePRA( u,d,m,e,p,N )
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
% p. 242
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
    Ui(1+p*sign(mod(i,p))-mod(i,p),:) = [ u(i:-1:max(i-N+1,1)),zeros(1,N-i) ];    
    if ~mod(i,p)
        di = transpose(d(i:-1:i-p+1));
        Yi = Ui*w;
        y(i:-1:i-p+1) = transpose(Yi);
        w = w + m*Ui'*inv(e*eye(p)+Ui*Ui')*(di-Yi);
    end;
end;
y(p*floor(uLen/p)+1:uLen) = d(p*floor(uLen/p)+1:uLen);
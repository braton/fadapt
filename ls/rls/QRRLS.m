function [w,y] = QRRLS( u,d,a,e,N )
% QRRLS - QR-Recursive Least Squares
%
% [w,y] = QRRLS( u,d,a,e,N )
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
% p. 788
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
winit=zeros(N,1);
w=winit;
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
O = sqrt(e)*eye(N);
q = zeros(N,1);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    [Q,R] = qr( ([sqrt(a)*O,ui';sqrt(a)*q',(d(i)-ui*winit)'])' );
    O = R(1:N,1:N)';
    q = R(1:N,N+1);
    w = inv(O)'*q + O'*winit;
    %{
    % Fastest for real data and taps initialized to zero
    [Q,R]=qr([sqrt(a)*O,sqrt(a)*q;ui,d(i)]);
    O=R(1:N,1:N);q=R(1:N,N+1);
    w=O\q;
    %}
end
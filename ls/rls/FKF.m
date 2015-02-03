function [w,y] = FKF( u,d,a,e,N )
% FKF - Fast Kalman Filter
%
% [w,y] = FKF( u,d,a,e,N )
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
% p. 838
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
uN = zeros(1,N);
gNe = zeros(N+1,1);
wf = zeros(N,1);
wb = zeros(N,1);
rf = e*a^-2;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    fpri = u(i)-uN*wf;
    wf = wf+fpri*gNe(1:N);
    fpos = u(i)-uN*wf;
    rf = a*rf+fpri'*fpos;
    gNe = [0;gNe(1:N)]+(fpos'/rf)*[1;-wf];
    uiN = uN(N);uN = [u(i),uN(1:N-1)];
    bpri = uiN-uN*wb;
    gNe(1:N) = (gNe(1:N)+gNe(N+1)*wb)/(1-gNe(N+1)*bpri);
    wb = wb+bpri*gNe(1:N);
    y(i) = uN*w;
    w = w+gNe(1:N)*(d(i)-y(i));
end
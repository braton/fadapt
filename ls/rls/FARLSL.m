function [y] = FARLSL( u,d,a,e,N )
% FARLSL - Fast Array Recursive Least Squares Lattice
%          based on QR decomposition
%
% [y] = FARLSL( u,d,a,e,N )
%
% u         - Input signal;
% d         - desired signal;
% a         - forgetting factor;
% e         - regularization factor;
% N         - number of taps.
%
% y         - Output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Jun-2006$
% $Revision: 04-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 920
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u=u.'; end
s = size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
uLen = length(u);
Pf=sqrt(e*a^(-2))*ones(1,N);
Pb=sqrt(e*a.^(-2:-1:-N-1));Pbp=Pb;
q=zeros(1,N);qf=zeros(1,N);qb=zeros(1,N);
y = zeros(1,uLen);
f=zeros(1,N+1);bp=zeros(1,N+1);b=zeros(1,N+1);
r=zeros(1,N+1);gamma=ones(1,N+1);
cb=zeros(1,N);sb=zeros(1,N);
cbp=zeros(1,N);sbp=zeros(1,N);
cf=zeros(1,N);sf=zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gamma(1)=1;b(1)=u(i);f(1)=u(i);r(1)=d(i);
    for m=1:N
        % Backward prediction     
        Pbpp=Pbp(m);Pbp(m)=sqrt(a*Pbp(m)^2+abs(bp(m))^2);
        cbp(m)=sqrt(a)*Pbpp/Pbp(m);sbp(m)=bp(m)'/Pbp(m);
        f(m+1)=f(m)*cbp(m)-sqrt(a)*qf(m)*sbp(m)';
        qf(m)=sqrt(a)*cbp(m)*qf(m)+sbp(m)*f(m);
        % Joint estimation
        Pb(m)=sqrt(a*Pbp(m)^2+abs(b(m))^2);
        cb(m)=sqrt(a)*Pbp(m)/Pb(m);sb(m)=b(m)'/Pb(m);
        r(m+1)=r(m)*cb(m)-sqrt(a)*q(m)*sb(m)';
        q(m)=sqrt(a)*cb(m)*q(m)+sb(m)*r(m);
        gamma(m+1)=gamma(m)*cb(m);
        % Forward prediction
        Pfp=Pf(m);Pf(m)=sqrt(a*Pf(m)^2+abs(f(m))^2);
        cf(m)=sqrt(a)*Pfp/Pf(m);sf(m)=f(m)'/Pf(m);
        b(m+1)=cf(m)*bp(m)-sqrt(a)*sf(m)'*qb(m);
        qb(m)=sqrt(a)*cf(m)*qb(m)+sf(m)*bp(m);
    end
    bp=b;
    y(i)=d(i)-r(N+1)'/gamma(N+1);
end
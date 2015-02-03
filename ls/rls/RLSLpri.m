function [y] = RLSLpri( u,d,a,e,N )
% RLSLpri - A Priori Based Recursive Least Squares Lattice
%
% [y] = RLSLpri( u,d,a,e,N )
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
% $Revision: 03-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 895
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
uLen = length(u);
Pf=e*a^(-2)*ones(1,N);
Pb=e*a.^(-2:-1:-N-1);
gconvp=ones(1,N+1);bprip=zeros(1,N+1);
ddot=zeros(1,N);pdot=zeros(1,N);
y = zeros(1,uLen);
gconv=ones(1,N+1);bpri=zeros(1,N+1);
fpri=zeros(1,N+1);epri=zeros(1,N+1);
k=zeros(1,N);kf=zeros(1,N);kb=zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gconv(1)=1;bpri(1)=u(i);fpri(1)=u(i);epri(1)=d(i);
    for m=1:N
        Pf(m)=a*Pf(m)+abs(fpri(m))^2*gconvp(m);
        Pbp=Pb(m);Pb(m)=a*Pb(m)+abs(bpri(m))^2*gconv(m);
        ddot(m)=a*ddot(m)+fpri(m)'*bprip(m)*gconvp(m);
        pdot(m)=a*pdot(m)+epri(m)'*bpri(m)*gconv(m);
        bpri(m+1)=bprip(m)-kb(m)*fpri(m);
        fpri(m+1)=fpri(m)-kf(m)*bprip(m);
        epri(m+1)=epri(m)-k(m)*bpri(m);
        gconv(m+1)=gconv(m)-abs(gconv(m)*bpri(m))^2/Pb(m);
        kb(m)=ddot(m)/Pf(m);
        kf(m)=ddot(m)'/Pbp;
        k(m)=pdot(m)'/Pb(m);
    end
    gconvp=gconv;bprip=bpri;
    y(i)=d(i)-epri(N+1);
end
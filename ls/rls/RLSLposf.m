function [y] = RLSLposf( u,d,a,e,N )
% RLSLposf - A Posteriori Based Error Feedback 
%            Recursive Least Squares Lattice
%
% [y] = RLSLposf( u,d,a,e,N )
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
% p. 902
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
uLen = length(u);
Pf=e*a^(-2)*ones(1,N);
Pb=e*a.^(-2:-1:-N-1);
Pbp=Pb;Pbpp=Pbp;
gconvp=ones(1,N+1);bposp=zeros(1,N+1);
y = zeros(1,uLen);
gconv=ones(1,N+1);bpos=zeros(1,N+1);
fpos=zeros(1,N+1);rpos=zeros(1,N+1);
k=zeros(1,N);kf=zeros(1,N);kb=zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gconv(1)=1;bpos(1)=u(i);fpos(1)=u(i);rpos(1)=d(i);
    for m=1:N
        Pfp=Pf(m);Pf(m)=a*Pf(m)+abs(fpos(m))^2/gconvp(m);
        Pb(m)=a*Pb(m)+abs(bpos(m))^2/gconv(m);
        gconv(m+1)=gconv(m)-abs(bpos(m))^2/Pb(m);
        k(m)=(gconv(m+1)/gconv(m))*(k(m)+bpos(m)'*rpos(m)/(a*gconv(m)*Pbp(m)));
        kf(m)=(gconvp(m+1)/gconvp(m))*(kf(m)+bposp(m)'*fpos(m)/(a*gconvp(m)*Pbpp(m)));
        kb(m)=(gconv(m+1)/gconvp(m))*(kb(m)+fpos(m)'*bposp(m)/(a*gconvp(m)*Pfp));
        rpos(m+1)=rpos(m)-k(m)*bpos(m);
        bpos(m+1)=bposp(m)-kb(m)*fpos(m);
        fpos(m+1)=fpos(m)-kf(m)*bposp(m);        
    end
    gconvp=gconv;bposp=bpos;Pbpp=Pbp;Pbp=Pb;
    y(i)=d(i)-rpos(N+1)/gconv(N+1);
end
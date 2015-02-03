function [y] = NRLSL( u,d,a,e,N )
% NRLSL - Normalized Recursive Least Squares Lattice
%
% [y] = NRLSL( u,d,a,e,N )
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
% p. 904
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
uLen = length(u);
Pb=e*a^(-2);
P=e*a^(-2);
bn=zeros(1,N+1);
y = zeros(1,uLen);
fn=zeros(1,N+1);rn=zeros(1,N+1);bnp=zeros(1,N+1);
sgm=zeros(1,N+1);gconv=zeros(1,N+1);
rpos=zeros(1,N+1);epri=zeros(1,N+1);
p=zeros(1,N);pb=zeros(1,N);pf=zeros(1,N);pbp=ones(1,N);
ka=zeros(1,N);kc=zeros(1,N);pa=zeros(1,N);pc=zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    Pb=a*Pb+abs(u(i))^2;P=a*P+abs(d(i))^2;
    bn(1)=u(i)/sqrt(Pb);fn(1)=bn(1);rn(1)=d(i)/sqrt(P);
    sgm(1)=sqrt(P);gconv(1)=1;
    for m=1:N
        pb(m)=sqrt(1-abs(bn(m))^2);
        pf(m)=sqrt(1-abs(fn(m))^2);
        p(m)=sqrt(1-abs(rn(m))^2);
        ka(m)=ka(m)*pbp(m)*pf(m)+fn(m)*bnp(m)';
        kc(m)=kc(m)*pb(m)*p(m)+bn(m)'*rn(m);
        pa(m)=sqrt(1-abs(ka(m))^2);
        pc(m)=sqrt(1-abs(kc(m))^2);
        rn(m+1)=(pb(m)*pc(m))^(-1)*(rn(m)-kc(m)*bn(m));
        bn(m+1)=(pf(m)*pa(m))^(-1)*(bnp(m)-ka(m)'*fn(m));
        fn(m+1)=(pbp(m)*pa(m))^(-1)*(fn(m)-ka(m)*bnp(m));
        sgm(m+1)=sgm(m)*pc(m)*pb(m);
        gconv(m+1)=gconv(m)*(pb(m)^2);
        rpos(m+1)=sgm(m+1)*rn(m+1);
        epri(m+1)=rpos(m+1)/gconv(m+1);
    end
    pbp=pb;bnp=bn;
    y(i)=d(i)-epri(N+1);
end
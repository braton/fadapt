function [y] = SWRLSLposf( u,d,e,L,N )
% SWRLSLposf - Sliding Window A Posteriori Error Feedback 
%              Recursive Least Squares Lattice
%
% [y] = SWRLSLposf( u,d,e,L,N )
%
% u         - Input signal;
% d         - desired signal;
% e         - regularization factor;
% L         - memory length;
% N         - number of taps.
%
% y         - Output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Jun-2006$
% $Revision: 03-Nov-2006$
%
% References:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% ch. 11-15 
%
% K. ZHAO, F. LING, H. LEV-ARI, J.G. PROAKIS, "Sliding Window
% Order-Recursive Least-Squares Algorithms", IEEE Transactions
% On Signal Processing, vol. 42, no. 8, August 1994.
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
uLen=length(u);
Pfu=e*ones(1,N);Pbd=Pfu;Pbu=Pfu;Pbdp=Pfu;Pfd=Pfu;
gconvdp=ones(1,N+1);gconvup=gconvdp;
bposdp=zeros(1,N+1);bposup=bposdp;
kd=zeros(1,N);kbd=kd;kfd=kd;
y=zeros(1,uLen);
Pfdp=zeros(1,N);Pbdpp=Pfdp;Pbup=Pfdp;
gconvd=zeros(1,N+1);gconvu=gconvd;
bposd=zeros(1,N+1);bposu=bposd;
fposd=bposd;fposu=bposd;rposd=bposd;rposu=bposd;
ku=zeros(1,N);kfu=ku;kbu=ku;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gconvd(1)=1;gconvu(1)=1;
    bposd(1)=(i-L+1>0)*u(max(i-L+1,1));bposu(1)=u(i);
    fposd(1)=(i-L+1>0)*u(max(i-L+1,1));fposu(1)=u(i);
    rposd(1)=(i-L+1>0)*d(max(i-L+1,1));rposu(1)=d(i);
    for m=1:N
        Pbup(m)=Pbu(m);Pbu(m)=Pbd(m)+abs(bposu(m))^2/gconvu(m);
        Pbdpp(m)=Pbdp(m);Pbdp(m)=Pbd(m);Pbd(m)=Pbu(m)-abs(bposd(m))^2/gconvd(m);
        Pfu(m)=Pfd(m)+abs(fposu(m))^2/gconvup(m);  
        Pfdp(m)=Pfd(m);Pfd(m)=Pfu(m)-abs(fposd(m))^2/gconvdp(m);
        gconvu(m+1)=gconvu(m)-abs(bposu(m))^2/Pbu(m);
        gconvd(m+1)=gconvd(m)+abs(bposd(m))^2/Pbd(m);
        ku(m)=(gconvu(m+1)/gconvu(m))*(kd(m)+bposu(m)'*rposu(m)/(gconvu(m)*Pbdp(m)));
        kd(m)=(gconvd(m+1)/gconvd(m))*(ku(m)-bposd(m)'*rposd(m)/(gconvd(m)*Pbu(m)));
        kbu(m)=(gconvu(m+1)/gconvup(m))*(kbd(m)+fposu(m)'*bposup(m)/(gconvup(m)*Pfdp(m)));
        kbd(m)=(gconvd(m+1)/gconvdp(m))*(kbu(m)-fposd(m)'*bposdp(m)/(gconvdp(m)*Pfu(m)));
        kfu(m)=(gconvup(m+1)/gconvup(m))*(kfd(m)+bposup(m)'*fposu(m)/(gconvup(m)*Pbdpp(m)));
        kfd(m)=(gconvdp(m+1)/gconvdp(m))*(kfu(m)-bposdp(m)'*fposd(m)/(gconvdp(m)*Pbup(m)));
        rposu(m+1)=rposu(m)-ku(m)*bposu(m);
        rposd(m+1)=rposd(m)-kd(m)*bposd(m);
        bposu(m+1)=bposup(m)-kbu(m)*fposu(m);
        bposd(m+1)=bposdp(m)-kbd(m)*fposd(m);
        fposu(m+1)=fposu(m)-kfu(m)*bposup(m);
        fposd(m+1)=fposd(m)-kfd(m)*bposdp(m);
    end
    gconvdp=gconvd;gconvup=gconvu;
    bposdp=bposd;bposup=bposu;
    y(i)=d(i)-rposu(N+1)/gconvu(N+1);
end
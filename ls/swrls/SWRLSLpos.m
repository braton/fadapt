function [y] = SWRLSLpos( u,d,e,L,N )
% SWRLSLpos - Sliding Window A Posteriori Based
%             Recursive Least Squares Lattice
%
% [y] = SWRLSLpos( u,d,e,L,N )
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
uLen = length(u);
Pfu=e*ones(1,N);Pbd=Pfu;Pbu=Pfu;
gconvpd=ones(1,N+1);gconvpu=gconvpd;
pdotu=zeros(1,N);ddotu=pdotu;
bpospd=zeros(1,N+1);bpospu=bpospd;
y=zeros(1,uLen);
Pfd=zeros(1,N);Pbdp=Pfd;Pbup=Pfd;
pdotd=zeros(1,N);ddotd=pdotd;
gconvd=zeros(1,N+1);gconvu=gconvd;
bposd=zeros(1,N+1);bposu=bposd;
fposd=bposd;fposu=bposd;rposd=bposd;rposu=bposd;
kd=zeros(1,N);kfd=kd;kbd=kd;ku=kd;kfu=kd;kbu=kd;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gconvd(1)=1;gconvu(1)=1;
    bposd(1)=(i-L>0)*u(max(i-L,1));bposu(1)=u(i);
    fposd(1)=(i-L>0)*u(max(i-L,1));fposu(1)=u(i);
    rposd(1)=(i-L>0)*d(max(i-L,1));rposu(1)=d(i);
    for m=1:N
        % Time Update
        Pfd(m)=Pfu(m)-abs(fposd(m))^2/gconvpd(m);
        Pfu(m)=Pfd(m)+abs(fposu(m))^2/gconvpu(m);
        Pbdp(m)=Pbd(m);Pbd(m)=Pbu(m)-abs(bposd(m))^2/gconvd(m);
        Pbup(m)=Pbu(m);Pbu(m)=Pbd(m)+abs(bposu(m))^2/gconvu(m);
        pdotd(m)=pdotu(m)-rposd(m)'*bposd(m)/gconvd(m);
        pdotu(m)=pdotd(m)+rposu(m)'*bposu(m)/gconvu(m);        
        ddotd(m)=ddotu(m)-fposd(m)'*bpospd(m)/gconvpd(m);        
        ddotu(m)=ddotd(m)+fposu(m)'*bpospu(m)/gconvpu(m);        
        % Reflection Coefficients
        kd(m)=pdotd(m)'/Pbd(m);ku(m)=pdotu(m)'/Pbu(m);
        kbd(m)=ddotd(m)/Pfd(m);kbu(m)=ddotu(m)/Pfu(m);
        kfd(m)=ddotd(m)'/Pbdp(m);kfu(m)=ddotu(m)'/Pbup(m);
        % Order Update
        gconvd(m+1)=gconvd(m)+abs(bposd(m))^2/Pbd(m);
        gconvu(m+1)=gconvu(m)-abs(bposu(m))^2/Pbu(m);
        bposd(m+1)=bpospd(m)-kbd(m)*fposd(m);
        bposu(m+1)=bpospu(m)-kbu(m)*fposu(m);
        fposd(m+1)=fposd(m)-kfd(m)*bpospd(m);
        fposu(m+1)=fposu(m)-kfu(m)*bpospu(m);
        rposd(m+1)=rposd(m)-kd(m)*bposd(m);
        rposu(m+1)=rposu(m)-ku(m)*bposu(m);
    end
    gconvpd=gconvd;gconvpu=gconvu;
    bpospd=bposd;bpospu=bposu;
    y(i)=d(i)-rposu(N+1)/gconvu(N+1);
end
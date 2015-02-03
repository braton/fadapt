function [y] = SWFARLSL( u,d,e,L,N )
% SWFARLSL - Sliding Window Fast Array Recursive Least Squares Lattice
%	     Based on QR Decomposition
%
% [y] = SWFARLSL( u,d,e,L,N )
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
% $Date: Aug-2006$
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
Pfd=sqrt(e)*ones(1,N);Pbdp=Pfd;
qd=zeros(1,N);qbd=qd;qfd=qd;
bposdp=zeros(1,N+1);bposup=bposdp;
gammaup=ones(1,N+1);
y=zeros(1,uLen);
Pfu=zeros(1,N);Pbd=Pfu;Pbu=Pfu;Pbup=Pbu;
qu=zeros(1,N);qbu=qu;qfu=qu;
bposd=zeros(1,N+1);bposu=bposd;
fposd=bposd;fposu=bposd;rposd=bposd;rposu=bposd;
gammad=zeros(1,N+1);gammau=gammad;
% Filtering ---------------------------------------------------------------
for i=1:uLen
    gammad(1)=1;gammau(1)=1;
    bposd(1)=(i-L+1>0)*u(max(i-L+1,1));bposu(1)=u(i);
    fposd(1)=(i-L+1>0)*u(max(i-L+1,1));fposu(1)=u(i);
    rposd(1)=(i-L+1>0)*d(max(i-L+1,1));rposu(1)=d(i);
    for m=1:N
        %-------------------------------------------%
        % Backward Update
        Pbup(m)=sqrt(Pbdp(m)^2+abs(bposup(m))^2);
        cNb=Pbdp(m)/Pbup(m);sNb=bposup(m)'/Pbup(m);
        qfu(m)=qfd(m)*cNb+fposu(m)'*sNb';
        fposu(m+1)=fposu(m)*cNb-qfd(m)'*sNb';
        % Backward Downdate
        Pbdp(m)=sqrt(Pbup(m)^2-abs(bposdp(m))^2);
        chNb=Pbup(m)/Pbdp(m);shNb=bposdp(m)'/Pbdp(m);
        qfd(m)=qfu(m)*chNb-fposd(m)'*shNb';
        fposd(m+1)=fposd(m)*chNb-qfu(m)'*shNb';
        %-------------------------------------------%
        % Joint Update
        Pbu(m)=sqrt(Pbdp(m)^2+abs(bposu(m))^2);
        cNr=Pbdp(m)/Pbu(m);sNr=bposu(m)'/Pbu(m);
        qu(m)=qd(m)*cNr+rposu(m)'*sNr';
        rposu(m+1)=rposu(m)*cNr-qd(m)'*sNr';
        % Joint Downdate
        Pbd(m)=sqrt(Pbu(m)^2-abs(bposd(m))^2);
        chNr=Pbu(m)/Pbd(m);shNr=bposd(m)'/Pbd(m);
        qd(m)=qu(m)*chNr-rposd(m)'*shNr';
        rposd(m+1)=rposd(m)*chNr-qu(m)'*shNr';
        %-------------------------------------------%
        % Forward Update
        Pfu(m)=sqrt(Pfd(m)^2+abs(fposu(m))^2);
        cNf=Pfd(m)/Pfu(m);sNf=fposu(m)'/Pfu(m);
        qbu(m)=qbd(m)*cNf+bposup(m)'*sNf';
        bposu(m+1)=bposup(m)*cNf-qbd(m)'*sNf';
        gammau(m+1)=gammaup(m)*cNf;
        % Forward Downdate
        Pfd(m)=sqrt(Pfu(m)^2-abs(fposd(m))^2);
        chNf=Pfu(m)/Pfd(m);shNf=fposd(m)'/Pfd(m);
        qbd(m)=qbu(m)*chNf-bposdp(m)'*shNf';
        bposd(m+1)=bposdp(m)*chNf-qbu(m)'*shNf';
    end
    bposdp=bposd;bposup=bposu;gammaup=gammau;
    y(i)=d(i)-rposu(N+1)'/gammau(N+1);
end
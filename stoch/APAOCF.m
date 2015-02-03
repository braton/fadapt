function [w,y] = APAOCF( u,d,m,e,p,N )
% APAOCF - Affine Projection Algorithm with Orthogonal Correction Factors
%
% [w,y] = APAOCF( u,d,m,e,p,N )
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
% $Date: 07-Mar-2006$
% $Revision: 02-Nov-2006$
%
% Reference:
% S.G. SANKARAN, A.A.(LOUIS) BEEX, "Normalized LMS algorithm with
% orthogonal correction factors", Proc. Asilomar Conf. on Signals,
% Systems, and Computers, Pacific Grove, CA, pp. 1670-1673, 1997.
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen=length(u);
y=zeros(1,uLen);
Ui=e*ones(p,N);
Uig=e*ones(p,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    Ui(2:p,:)=Ui(1:p-1,:);
    Ui(1,:)=[ u(i:-1:max(i-N+1,1)),zeros(1,N-i) ];
    di=[ d(i:-1:max(i-p+1,1)),zeros(1,p-i) ].';
    Yi=Ui*w;
    y(i)=Yi(1);
    wi=w;
    S=zeros(N,N);
    for k=0:p-1
        Uig(k+1,:)=Ui(k+1,:)*(eye(N)-S);
        err=di(k+1)-Ui(k+1,:)*wi;
        Un=(Uig(k+1,:)*Uig(k+1,:)');
        if Un>e, Si=Uig(k+1,:)'/Un; else Si=zeros(N,1); end;
        wi=wi+Si*err;
        S=S+Si*Uig(k+1,:);
    end;
    w=(1-m)*w + m*wi;
end;
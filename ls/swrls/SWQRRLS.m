function [w,y] = SWQRRLS( u,d,e,L,N )
% SWQRRLS - Sliding Window QR Recursive Least Squares
%
% [w,y] = SWQRRLS( u,d,e,L,N )
%
% u         - Input signal;
% d         - desired signal;
% e         - regularization factor;
% L         - memory length;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Jun-2006$
% $Revision: 03-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% ch. 13,14
%

function [c,s] = CircGivens(a,b)
if b==0, c=1;s=0; elseif a==0, c=0;s=1; else
    p = b/a; c=1/sqrt(1+abs(p)^2); s=p*c;
end
end

function [ch,sh] = HypGivens(a,b)
if b==0, ch=1;sh=0; else
    if abs(b)<abs(a), p = b/a; else p = (a/b)'; end
    ch=1/sqrt(1-abs(p)^2); sh=p*ch;
end
end

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
O = sqrt(e)*eye(N);
qu = zeros(N,1);
ui = zeros(1,N);
uiL = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    uiL = [ (i-L>0)*u(max(i-L,1)),uiL(1:N-1) ];
    ui = [ u(i),ui(1:N-1) ];
    % Downdate
    M=[O,qu;uiL,(i-L>0)*d(max(i-L,1))];
    for j=N:-1:1,
        [ch,sh]=HypGivens(M(j,j),M(N+1,j));
        Mj=M(j,:)*ch-M(N+1,:)*conj(sh);
        M(N+1,:)=M(N+1,:)*ch-M(j,:)*sh;M(j,:)=Mj;
    end
    O=M(1:N,1:N);qd=M(1:N,N+1);
    w = O\qd;
    % Update
    y(i)=ui*w;
    M=[O,qd;ui,d(i)];
    for j=N:-1:1
        [c,s]=CircGivens(M(j,j),M(N+1,j));
        Mj=M(j,:)*c+M(N+1,:)*conj(s);
        M(N+1,:)=M(N+1,:)*c-M(j,:)*s;M(j,:)=Mj;
    end
    O=M(1:N,1:N);qu=M(1:N,N+1);
    w=O\qu;
end
end
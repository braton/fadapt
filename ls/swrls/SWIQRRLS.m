function [w,y] = SWIQRRLS( u,d,e,L,N )
% SWIQRRLS - Sliding Window Inverse QR Recursive Least Squares
%
% [w,y] = SWIQRRLS( u,d,e,L,N )
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
% p. 851
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
P = sqrt(1/e)*eye(N);
y = zeros(1,uLen);
ui = zeros(1,N);
uiL = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    uiL = [ (i-L>0)*u(max(i-L,1)),uiL(1:N-1) ];
    ui = [ u(i),ui(1:N-1) ];
    % Downdate
    eL = (i-L>0)*d(max(i-L,1))-uiL*w;
    M=[1,uiL*P;zeros(N,1),P];
    for j=N+1:-1:2,
        [ch,sh]=HypGivens(M(1,1),M(1,j));
        M1=M(:,1)*ch-M(:,j)*conj(sh);
        M(:,j)=M(:,j)*ch-M(:,1)*sh;M(:,1)=M1;
    end
    P=M(2:N+1,2:N+1);
    w = w - (M(2:N+1,1)/M(1,1))*eL;
    % Update
    y(i)=ui*w;
    e = d(i)-y(i);
    M=[1,ui*P;zeros(N,1),P];
    for j=N+1:-1:2
        [c,s]=CircGivens(M(1,1),M(1,j));
        M1=M(:,1)*c+M(:,j)*conj(s);
        M(:,j)=M(:,j)*c-M(:,1)*s;M(:,1)=M1;
    end
    P=M(2:N+1,2:N+1);
    w = w + (M(2:N+1,1)/M(1,1))*e;
end
end
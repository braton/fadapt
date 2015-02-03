function [w,y] = SWFARLS( u,d,e,L,N )
% SWFARLS - Sliding Window Fast Array Recursive Least Squares
%
% [w,y] = SWFARLS( u,d,e,L,N )
%
% u         - Input signal;
% d         - desired signal;
% e         - regularization factor;
% L         - block length;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Jun-2006$
% $Revision: 04-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% ch. 11-15 
%
% K. ZHAO, F. LING, H. LEV-ARI, J.G. PROAKIS, "Sliding Window
% Order-Recursive Least-Squares Algorithms", IEEE Transactions
% On Signal Processing, vol. 42, no. 8, August 1994.
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
uLen=length(u);
y=zeros(1,uLen);
ui=zeros(1,N);
ul=zeros(1,N);
gammau=1; ggammau=zeros(N,1);
gammad=1; ggammad=zeros(N,1);
M=zeros(N+1,2);
M(1,1)=sqrt(1/e);M(N+1,2)=sqrt(1/e);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    % Update
    E = [ [gammau;0;ggammau],[[u(i),ui]*M;M] ];
    [c,s] = CircGivens(E(1,1),E(1,2));
    E = E*[c,s,0;conj(s),-c,0;0,0,1];
    [ch,sh] = HypGivens(E(1,1),E(1,3));
    E = E*[ch,0,-sh;0,1,0;-conj(sh),0,ch];
    gammau = E(1,1); ggammau = E(2:N+1,1);
    M = E(2:N+2,2:3);
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    w = w + (ggammau/gammau)*(d(i)-y(i));
    % Downdate
    E = [ [gammad;0;ggammad],[[ (i-L>0)*u(max(i-L,1)),ul ]*M;-M] ];
    [ch,sh] = HypGivens(E(1,1),E(1,2));
    E = E*[ch,sh,0;-conj(sh),-ch,0;0,0,1];
    [c,s] = CircGivens(E(1,1),E(1,3));
    E = E*[c,0,-s;0,1,0;conj(s),0,c];
    gammad = E(1,1); ggammad = E(2:N+1,1);
    M = -E(2:N+2,2:3);
    ul=[ (i-L>0)*u(max(i-L,1)),ul(1:N-1) ];
    w = w - (ggammad/gammad)*((i-L>0)*d(max(i-L,1))-ul*w);
end
end
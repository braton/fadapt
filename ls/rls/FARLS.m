function [w,y] = FARLS( u,d,a,e,N )
% FARLS - Fast Array Recursive Least Squares
%
% [w,y] = FARLS( u,d,a,e,N )
%
% u         - Input signal;
% d         - desired signal;
% a         - forgetting factor;
% e         - regularization factor;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: Mar-2006$
% $Revision: 03-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 
%

function [c,s] = GivensRot(a,b)
if b==0, c=1;s=0; elseif a==0, c=0;s=1; else
    p = b/a; c=1/sqrt(1+abs(p)^2); s=p*c;
end
end

function [ch,sh] = HyperbolicRot(a,b)
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
ui = zeros(1,N);
gamma = 1; ggamma = zeros(N,1);
L = zeros(N+1,2);
L(1,1) = sqrt(a/e); L(N+1,2) = sqrt(a/e)*a^(N/2);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    E = [ [gamma;0;ggamma],[[u(i),ui]*L;L] ];
    [c,s] = GivensRot(E(1,1),E(1,2));
    E = E*[c,s,0;conj(s),-c,0;0,0,1];
    [ch,sh] = HyperbolicRot(E(1,1),E(1,3));
    E = E*[ch,0,-sh;0,1,0;-conj(sh),0,ch];
    gamma = E(1,1); ggamma = E(2:N+1,1);
    L = E(2:N+2,2:3)/sqrt(a);
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    w = w + (ggamma/gamma)*(d(i)-y(i));
end
end
function [w,y] = uDHTBpNLMS( u,d,m,b,B,N )
% uDHTBpNLMS - Unconstrained DHT Block Power Normalized Least Mean Square
%
% [w,y] = uDHTBpNLMS( u,d,m,b,B,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% b         - smoothing factor;
% B         - block length;
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
% p. 648
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
A = ceil(N/B);
AB = ceil(B/2);
K = 3*B-1 + 2*AB-1;
uLen = length(u);
uB = zeros(2*B,1);
uBH = zeros(K,A);
p = zeros(K,1);
lkc = zeros(K,A);
y = zeros(1,uLen);
% DHT Matrix
H=zeros(K);
for i=1:K
    for j=1:K
        H(i,j)=(cos(2*pi*(i-1)*(j-1)/K)-sin(2*pi*(i-1)*(j-1)/K))/sqrt(K);
    end;
end;
% Filtering ---------------------------------------------------------------
for i=1:ceil(uLen/B)
    uB = [ [zeros(1,i*B-uLen),u(min(i*B,uLen):-1:(i-1)*B+1)].';uB(1:B) ];
    dB = [zeros(1,i*B-uLen),d(min(i*B,uLen):-1:(i-1)*B+1)].';
    uBH = [ H*([zeros(AB,1);uB(1:2*B-1);zeros(B+AB-1,1)]),...
        uBH(:,1:length(uBH(1,:))-1) ];
    p = b*p + (1-b)*abs(uBH(:,1)).^2;
    yB = H*(sum(uBH.*lkc,2));
    y(min(i*B,uLen):-1:(i-1)*B+1) = yB(1+(i*B-uLen)*(i*B-uLen>=0):B);    
    ek = H*([ (dB-yB(1:B));zeros(K-B,1) ]);
    for j=1:K
        lkc(j,:) = lkc(j,:) + m*(ek(j)/p(j))*conj(uBH(j,:));
    end;
end;
% Coefficients transformation ---------------------------------------------
G = H*lkc/sqrt(K);
P = (G(AB+1:AB+B,:)+G(AB+3*B-1:-1:AB+2*B,:))/2;
for j = 1:ceil(N/B)
    w((j-1)*B+1:j*B) = P(:,j); 
end;
function [y] = BDHTconv( g,s,B )
% Block convolution via DHT
%
% [y] = BDHTconv( g,s,B )
%
% g         - Long filter taps;
% s         - input signal;
% B         - block length.
%
% y         - Output signal.
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
v = size(g); if v(1)>v(2), g = g.'; end;
v = size(s); if v(1)>v(2), s = s.'; end;
% Initialization ----------------------------------------------------------
N = length(g);
A = ceil(N/B);
sLen = length(s);
K = 3*B-1 + 2*ceil(B/2)-1;
% DHT Matrix
H=zeros(K);
for i=1:K
    for j=1:K
        H(i,j)=(cos(2*pi*(i-1)*(j-1)/K)-sin(2*pi*(i-1)*(j-1)/K))/sqrt(K);
    end;
end;
y = zeros( 1,sLen );
P = zeros( B,A );
uB = zeros(2*B,1);
uBH = zeros(K,A);
for i=1:A
    P(:,i) = [ g( (i-1)*B+1:min(i*B,N) ),zeros(1,i*B-N) ];
end;
L = H*([zeros(ceil(B/2),A);P;zeros(B-1,A);P(B:-1:1,:);zeros(ceil(B/2)-1,A)])*sqrt(K);
% Filtering ---------------------------------------------------------------
for i=1:ceil(sLen/B)
    uB = [ [zeros(1,i*B-sLen),s(min(i*B,sLen):-1:(i-1)*B+1)].';uB(1:B) ];
    uBH = [ H*([zeros(ceil(B/2),1);uB(1:2*B-1);zeros(B+ceil(B/2)-1,1)]),...
        uBH(:,1:length(uBH(1,:))-1) ];
    yB = H*(sum(uBH.*L,2));
    y(min(i*B,sLen):-1:(i-1)*B+1) = yB(1+(i*B-sLen)*(i*B-sLen>=0):B);
end;
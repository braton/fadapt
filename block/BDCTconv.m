function [y] = BDCTconv( g,s,B )
% Block convolution via DCT
%
% [y] = BDCTconv( g,s,B )
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
% p. 640
%

% Ensure row vectors ------------------------------------------------------
v = size(g); if v(1)>v(2), g = g.'; end;
v = size(s); if v(1)>v(2), s = s.'; end;
% Initialization ----------------------------------------------------------
N = length(g);
sLen = length(s);
if mod(B,2)
    K = 3.5*B-1.5;
    a = 1.5*B-0.5;
    n = 0.5*B+0.5;
else
    K = 3.5*B-2;
    a = 1.5*B-1;
    n = 0.5*B;
end;
y = zeros( 1,sLen );
P = zeros( B,ceil(N/B) );
uB = zeros(2*B,1);
uBC = zeros(K,ceil(N/B));
for i=1:ceil(N/B)
    P(:,i) = [ g( (i-1)*B+1:min(i*B,N) ),zeros(1,i*B-N) ];
end;
L = idct(sqrt(2*K)*[ zeros(B-1,ceil(N/B));P;zeros(a,ceil(N/B)) ]);
% Filtering ---------------------------------------------------------------
for i=1:ceil(sLen/B)
    uB = [ [zeros(1,i*B-sLen),s(min(i*B,sLen):-1:(i-1)*B+1)].';uB(1:B) ];
    uBC = [ idct([zeros(a,1);uB(1:2*B-1)]),uBC(:,1:length(uBC(1,:))-1) ];
    yB = dct(sum(uBC.*L,2));
    y(min(i*B,sLen):-1:(i-1)*B+1) = yB(n+1+(i*B-sLen)*(i*B-sLen>=0):n+B);
end;
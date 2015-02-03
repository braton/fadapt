function [y] = Bconv( g,s,B )
% Block convolution via DFT
%
% [y] = Bconv( g,s,B )
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
% p. 589
%

% Ensure row vectors ------------------------------------------------------
v = size(g); if v(1)>v(2), g = g.'; end;
v = size(s); if v(1)>v(2), s = s.'; end;
% Initialization ----------------------------------------------------------
N = length(g);
sLen = length(s);
y = zeros( 1,sLen );
P = zeros( 2*B,ceil(N/B) );
uB = zeros(2*B,1);
uBF = zeros(2*B,ceil(N/B));
for n=1:ceil(N/B)
    P(:,n) = [ g( (n-1)*B+1:min(n*B,N) ),zeros(1,n*B-N),zeros(1,B) ];
end;
L = ifft(P);
% Filtering ---------------------------------------------------------------
for n=1:ceil(sLen/B)
    uB = [ [zeros(1,n*B-sLen),s(min(n*B,sLen):-1:(n-1)*B+1)].';uB(1:B) ];
    uBF = [ fft(uB),uBF(:,1:length(uBF(1,:))-1) ];
    yB = ifft( sum(uBF.*L,2) )*2*B;
    y(min(n*B,sLen):-1:(n-1)*B+1) = yB(max(n*B-sLen+1,1):B);
end;
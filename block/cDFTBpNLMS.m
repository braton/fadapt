function [w,y] = cDFTBpNLMS( u,d,m,b,B,N )
% cDFTBpNLMS - Constrained DFT Block Power Normalized Least Mean Square
%              Also referred as Multidelay Filter (MDF)
%
% [w,y] = cDFTBpNLMS( u,d,m,b,B,N )
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
% p. 597
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
uB = zeros(2*B,1);
uBF = zeros(2*B,ceil(N/B));
p = zeros(2*B,1);
lkc = zeros(2*B,ceil(N/B));
y = zeros(1,uLen);
fn = 1/sqrt(2*B);
% Filtering ---------------------------------------------------------------
for n=1:ceil(uLen/B)
    uB = [ [zeros(1,n*B-uLen),u(min(n*B,uLen):-1:(n-1)*B+1)].';uB(1:B) ];
    dB = [zeros(1,n*B-uLen),d(min(n*B,uLen):-1:(n-1)*B+1)].';
    uBF = [ fft(uB),uBF(:,1:length(uBF(1,:))-1) ];
    p = b*p + (1-b)*abs(uBF(:,1)).^2;
    yB = ifft( sum(uBF.*lkc,2) )*2*B;
    y(min(n*B,uLen):-1:(n-1)*B+1) = yB(max(n*B-uLen+1,1):B);
    ek = fft([dB-yB(1:B);zeros(B,1)]);
    for i=1:2*B
        lkc(i,:) = lkc(i,:) + (m*fn)*(ek(i)/p(i))*conj(uBF(i,:));
    end;
    g = fft(lkc);
    % Put a constraint
    lkc = ifft( [g(1:B,:);zeros(B,ceil(N/B))] );
end;
% Coefficients transformation ---------------------------------------------
g = fft(lkc);
g = g(1:B,:);
for i=1:ceil(N/B)
    w((i-1)*B+1:i*B) = g(:,i); 
end;
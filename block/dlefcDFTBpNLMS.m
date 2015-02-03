function [w,y] = dlefcDFTBpNLMS( u,d,m,b,B,R,N )
% dlefcDFTBpNLMS - Delayless constrained DFT Block Power Normalized LMS
%                 with error estimate by efficient convolution in fullband
%
% [w,y] = dlefcDFTBpNLMS( u,d,m,b,B,R,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% b         - smoothing factor;
% B         - block length;
% R         - convolution block length;
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
% p. 628
%

% Ensure row vectors ------------------------------------------------------
s = size(u); if s(1)>s(2), u = u.'; end;
s = size(d); if s(1)>s(2), d = d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
iB = 1;
P = zeros( 2*R,ceil(N/R) );
L = zeros( 2*R,ceil(N/R) );
uR = zeros(2*R,1);
uRF = zeros(2*R,ceil(N/R));
y = zeros( 1,uLen );
eB = zeros(B,1);
uB = zeros(2*B,1);
uBF = zeros(2*B,ceil(N/B));
p = zeros(2*B,1);
lkc = zeros(2*B,ceil(N/B));
gdir = zeros(R,1);
uFb = zeros(1,R);
fn = 1/sqrt(2*B);
% Filtering ---------------------------------------------------------------
for n=1:ceil(uLen/R)
    % Polyphase components
    if imag(iB)
        for i=1:ceil(N/R)
            P(:,i) = [ w( (i-1)*R+1:min(i*R,N) ).',zeros(1,i*R-N),zeros(1,R) ];
        end;
        P(1:R-1,1) = zeros(R-1,1);
        L = ifft(P);
        gdir = [ g(1:R-1).';0 ];
        iB = iB - sqrt(-1);
    end;
    % Direct convolution with B-1 coefficients
    for i=1:R
        uFb = [ (((n-1)*R+i)<=uLen)*u(min((n-1)*R+i,uLen)),uFb(1:R-1) ];    
        y(min((n-1)*R+i,uLen)) = uFb*gdir;
    end;
    uR = [ [zeros(1,n*R-uLen),u(min(n*R,uLen):-1:(n-1)*R+1)].';uR(1:R) ];
    uRF = [ fft(uR),uRF(:,1:length(uRF(1,:))-1) ];
    yR = ifft( sum(uRF.*L,2) )*2*R;
    y(min(n*R,uLen):-1:(n-1)*R+1) = y(min(n*R,uLen):-1:(n-1)*R+1) + ...
    yR(max(n*R-uLen+1,1):R).';
    % Block adaptation error
    if n*R>=iB*B 
        eB(B:-1:1) = [d((iB-1)*B+1:min(iB*B,uLen))-y((iB-1)*B+1: ...
                     min(iB*B,uLen)),zeros(1,iB*B-uLen)];
        uB = [ [zeros(1,iB*B-uLen),u(min(iB*B,uLen):-1:(iB-1)*B+1)].';uB(1:B) ];
        uBF = [ fn*fft(uB),uBF(:,1:length(uBF(1,:))-1) ];
        ek = fn*fft( [eB;zeros(B,1)] );
        p = b*p + (1-b)*abs(uBF(:,1)).^2;
        for i=1:2*B
            lkc(i,:) = lkc(i,:) + m*(ek(i)/p(i))*conj(uBF(i,:));
        end;
        % Subband/fullband mapping
        g = fn*fft(lkc);
        g = g(1:B,:);
        for i=1:ceil(N/B)
            w((i-1)*B+1:i*B) = g(:,i); 
        end;
        iB=iB+1+sqrt(-1); 
    end;
end;
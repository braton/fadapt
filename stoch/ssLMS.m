function [w,y] = ssLMS( u,d,m,N )
% ssLMS - Sign-sign Least Mean Square
%
% [w,y] = ssLMS( u,d,m,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - step size;
% N         - number of taps.
%
% w         - Adapted (or not) taps;
% y         - output signal.
%
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 07-Mar-2006$
% $Revision: 02-Nov-2006$
%
% Reference:
% A.H. SAYED, "Fundamentals of Adaptive Filtering", John Wiley & Sons 2003
% p. 253
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end;
s=size(d); if s(1)>s(2), d=d.'; end;
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    ui = [ u(i),ui(1:N-1) ];
    y(i) = ui*w;
    w = w + m*(sign(real(ui'))+sqrt(-1)*sign(imag(ui')))* ...
    (  sign(real(d(i)-y(i)))+sqrt(-1)*sign(imag(d(i)-y(i)))  );
end;
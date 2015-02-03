function [w,y] = SWRLS( u,d,e,L,N )
% SWRLS - Sliding Window Recursive Least Squares
%         ( RLS with finite memory )
%
% [w,y] = SWRLS( u,d,e,L,N )
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
% p. 750
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
ui = zeros(1,N);
uiL = zeros(1,N);
Pu = (1/e)*eye(N);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    uiL = [ (i-L>0)*u(max(i-L,1)),uiL(1:N-1) ];
    ui = [ u(i),ui(1:N-1) ];
    % Downdate
    eL = (i-L>0)*d(max(i-L,1)) - uiL*w;
    gammad = (1-uiL*Pu*uiL')^(-1);
    gd = Pu*uiL'*gammad;
    w = w - gd*eL;
    Pd = Pu + (gd*gd')/gammad;
    % Update
    gammau = (1+ui*Pd*ui')^(-1);
    gu = Pd*ui'*gammau;
    y(i) = ui*w;
    w = w + gu*(d(i)-ui*w);
    Pu = Pd - (gu*gu')/gammau;
end
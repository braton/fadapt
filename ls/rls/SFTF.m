function [w,y] = SFTF( u,d,a,e,N )
% SFTF - Stabilized Fast Transversal Filters
%
% [w,y] = SFTF( u,d,a,e,N )
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
% p. 842
%

% Ensure row vectors ------------------------------------------------------
s=size(u); if s(1)>s(2), u=u.'; end
s=size(d); if s(1)>s(2), d=d.'; end
% Initialization ----------------------------------------------------------
w=zeros(N,1);
uLen = length(u);
y = zeros(1,uLen);
uN = zeros(1,N);
wfpos = zeros(N,1);
wb = zeros(N,1);
gamma = 1;
rfpos = e*a^-2;
cNe = zeros(N+1,1);
rb = e*a^(-N-2);
wfpri = zeros(N,1);
% Filtering ---------------------------------------------------------------
for i=1:uLen
    fpri = u(i)-uN*wfpri;
    fpos = gamma*fpri;
    rfpri = a*rfpos+fpri'*fpos;
    gamma = gamma*a*rfpos/rfpri;
    wfpri = wfpos+fpos*cNe(1:N);
    cNe = [0;cNe(1:N)]+(fpri'/(a*rfpos))*[1;-wfpos];
    bpri = a*rb*cNe(N+1)';
    gamma = gamma/(1-bpri*gamma*cNe(N+1));
    cNe = cNe-cNe(N+1)*[-wb;1];
    bpos = gamma*bpri;
    rb = a*rb+bpri'*bpos;
    wb = wb+bpos*cNe(1:N);
    uN = [u(i),uN(1:N-1)];
    y(i) = uN*w;
    epos = gamma*(d(i)-y(i));
    w = w+epos*cNe(1:N);
    wfpos = wfpri;
    rfpos = rfpri;
end
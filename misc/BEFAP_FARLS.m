function [y] = BEFAP_FARLS( u,d,s,e,p,L,N )
% BEFAP_FARLS - Block Exact Fast Affine Projection Algorithm
%               with Fast Array Recursive Least Squares prediction
%
% [y] = BEFAP_FARLS( u,d,m,e,p,L,N )
%
% u         - Input signal;
% d         - desired signal;
% m         - stepsize;
% e         - regularization factor;
% p         - projection order;
% L         - block length;
% N         - filter order.
%
% y         - Output signal.
%
% Bartosz Zator
% August 2006
%

function [c,s] = CircGivens(a,b)
if b==0, c=1;s=0; elseif a==0, c=0;s=1; else
    v = b/a; c=1/sqrt(1+abs(v)^2); s=v*c;
end
end

function [ch,sh] = HypGivens(a,b)
if b==0, ch=1;sh=0; else
    if abs(b)<abs(a), v = b/a; else v = (a/b)'; end
    ch=1/sqrt(1-abs(v)^2); sh=v*ch;
end
end

% -------------------------------------------------------------------------
% Initialization ----------------------------------------------------------
% -------------------------------------------------------------------------
% FARLS section -----------------------------------------------------------
uLen=length(u);y=zeros(1,uLen);w=zeros(N,1);
wp=zeros(p-1,1);
uip=zeros(1,p-1);
ulp=zeros(1,p-1);
gammaups=1; ggammaups=zeros(p-1,1);
gammadps=1; ggammadps=zeros(p-1,1);
M=zeros(p,2);
M(1,1)=sqrt(1/e);M(p,2)=sqrt(1/e);
% APA section -------------------------------------------------------------
ruu=zeros(L+p-2,1);ur=zeros(1,N+L+p-2);ul=zeros(1,N+L-1);
E=zeros(p+L-1,1);ei=zeros(p,1);gu=zeros(p-1,1);
% -------------------------------------------------------------------------
% Filtering ---------------------------------------------------------------
% -------------------------------------------------------------------------
for i=1:uLen
    ul=[u(i),ul(1:N+L-2)];
    if ~mod(i,L)
        f2=fft([w.',zeros(1,L)]);f1=fft([0,fliplr(ul)]);
        Y=real(ifft(f1.*f2));Y=flipud(Y(N+1:N+L)');
        for r=i-L+1:i
            % SWFARLS prediction ------------------------------------------
            % Update
            EM = [ [gammaups;0;ggammaups],[[u(r),uip]*M;M] ];
            [cp,sp] = CircGivens(EM(1,1),EM(1,2));
            EM = EM*[cp,sp,0;conj(sp),-cp,0;0,0,1];
            [chp,shp] = HypGivens(EM(1,1),EM(1,3));
            EM = EM*[chp,0,-shp;0,1,0;-conj(shp),0,chp];
            gammaups = EM(1,1); ggammaups = EM(2:p,1);
            M = EM(2:p+1,2:3);
            uip = [ u(r),uip(1:p-2) ];
            ypu = uip*wp;
            wp = wp + (ggammaups/gammaups)*(d(r)-ypu);
            % Downdate
            EM = [ [gammadps;0;ggammadps],[[ (r-N>0)*u(max(r-N,1)),ulp ]*M;-M] ];
            [chp,shp] = HypGivens(EM(1,1),EM(1,2));
            EM = EM*[chp,shp,0;-conj(shp),-chp,0;0,0,1];
            [cp,sp] = CircGivens(EM(1,1),EM(1,3));
            EM = EM*[cp,0,-sp;0,1,0;conj(sp),0,cp];
            gammadps = EM(1,1); ggammadps = EM(2:p,1);
            M = -EM(2:p+1,2:3);
            ulp=[ (r-N>0)*u(max(r-N,1)),ulp(1:p-2) ];
            wp = wp - (ggammadps/gammadps)*((r-N>0)*d(max(r-N,1))-ulp*wp);
            % Affine Projection -------------------------------------------
            ruu=ruu+u(r)*ur(1:L+p-2)'-ur(N)*ur(N+1:N+L+p-2)';
            ur=[u(r),ur(1:N+L+p-3)];
            k=(mod(r,L)+(mod(r,L)==0)*L);
            auxerr=d(r)-Y(L-k+1);
            err=auxerr-E(1:k+p-2)'*ruu(1:k+p-2);
            y(r)=d(r)-err;
            ei=[d(r)-y(r);(1-s)*ei(1:p-1)];
            Bi=M(:,2)*M(:,2)'*ei;
            Fi=M(:,1)*M(:,1)'*ei;
            gi=[0;(1-s)*gu]+Fi;
            gu=gi(1:p-1)-Bi(1:p-1);
            E=[0;E(1:L+p-2)]+[s*gi;zeros(L-1,1)];
        end%for r=i-L+1:i
        f1=fft([0,fliplr(ur(p:p+N+L-2))]);
        f2=fft([E(p:p+L-1).',zeros(1,N)]);
        wf=real(ifft(f1.*f2));
        w=w+flipud(wf(L+1:L+N)');
    end%if ~mod(i,L)
end%for i=1:uLen
end
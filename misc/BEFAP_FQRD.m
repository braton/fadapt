function [y] = BEFAP_FQRD( u,d,s,e,p,L,N )
% BEFAP_FQRD -  Block Exact Fast Affine Projection Algorithm
%               with Fast Array Recursive Least Squares Lattice prediction
%
% [y] = BEFAP_FQRD( u,d,m,e,p,L,N )
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

% -------------------------------------------------------------------------
% Initialization ----------------------------------------------------------
% -------------------------------------------------------------------------
uLen=length(u);y=zeros(1,uLen);w=zeros(N,1);
% FRLS section ------------------------------------------------------------
wfLd=zeros(p,1);wfLd(1)=sqrt(1/e);wbLd=zeros(p,1);wbLd(p)=sqrt(1/e);
Pfd=(sqrt(e^(-1))*ones(1,p));Pbdp=(sqrt(e^(-1))*ones(1,p));
qbd=(zeros(1,p));qfd=qbd;
bposdp=(zeros(1,p+1));bposup=bposdp;
gammaup=(ones(1,p+1));gammadp=gammaup;
Ggain=(zeros(p,1));Ggaini=Ggain;
% -------------------------------------------------------------------------
Pfu=(zeros(1,p));Pbup=Pfu;
qbu=(zeros(1,p));qfu=qbu;
bposd=(zeros(1,p+1));bposu=bposd;fposd=bposd;fposu=bposd;
gammad=(zeros(1,p+1));gammau=gammad;
cNf=rand(1,p);sNf=cNf;chNf=cNf;shNf=cNf;
% APA section -------------------------------------------------------------
ruu=zeros(L+p-2,1);ur=zeros(1,N+L+p-2);ul=zeros(1,N+L-1);
E=zeros(p+L-1,1);ei=zeros(p,1);gu=zeros(p-1,1);
% -------------------------------------------------------------------------
% Filtering ---------------------------------------------------------------
% -------------------------------------------------------------------------
for i=1:uLen
    %if ~mod(i,1000), disp(i); end;
    ul=[u(i),ul(1:N+L-2)];
    if ~mod(i,L)
        f2=fft([w.',zeros(1,L)]);f1=fft([0,fliplr(ul)]);
        Y=real(ifft(f1.*f2));Y=flipud(Y(N+1:N+L)');
        for r=i-L+1:i
            % Forward and Backward Prediction -----------------------------
            gammad(1)=(1);gammau(1)=(1);
            bposd(1)=(r-N+1>0)*u(max(r-N+1,1));bposu(1)=u(r);
            fposd(1)=(r-N+1>0)*u(max(r-N+1,1));fposu(1)=u(r);
            %disp([(r-N+1>0)*u(max(r-N+1,1)),u(r)]);pause;
            for m=1:p
                %-------------------------------------------%
                % Backward Update
                Pbup(m)=Pbdp(m)/sqrt(1+Pbdp(m)^2*bposup(m)^2);
                cNb=Pbup(m)/Pbdp(m);sNb=bposup(m)*Pbup(m);
                qfu(m)=qfd(m)*cNb+fposu(m)*sNb;
                fposu(m+1)=fposu(m)*cNb-qfd(m)*sNb;
                % Backward Downdate
                Pbdp(m)=Pbup(m)/sqrt(1-Pbup(m)^2*bposdp(m)^2);
                chNb=Pbdp(m)/Pbup(m);shNb=bposdp(m)*Pbdp(m);
                qfd(m)=qfu(m)*chNb-fposd(m)*shNb;
                fposd(m+1)=fposd(m)*chNb-qfu(m)*shNb;
                %-------------------------------------------%
                % Forward Update
                Pfu(m)=Pfd(m)/sqrt(1+Pfd(m)^2*fposu(m)^2);
                cNf(m)=Pfu(m)/Pfd(m);sNf(m)=fposu(m)*Pfu(m);
                qbu(m)=qbd(m)*cNf(m)+bposup(m)*sNf(m);
                bposu(m+1)=bposup(m)*cNf(m)-qbd(m)*sNf(m);
                gammau(m+1)=gammaup(m)*cNf(m);
                % Forward Downdate
                %Pfdp(m)=1/Pfd(m);
                Pfd(m)=Pfu(m)/sqrt(1-Pfu(m)^2*fposd(m)^2);
                chNf(m)=Pfd(m)/Pfu(m);shNf(m)=fposd(m)*Pfd(m);
                qbd(m)=qbu(m)*chNf(m)-bposdp(m)*shNf(m);
                bposd(m+1)=bposdp(m)*chNf(m)-qbu(m)*shNf(m);
                gammad(m+1)=gammadp(m)*chNf(m);
                %-------------------------------------------%                
            end%m=1:p
            % Update coefficients -----------------------------------------
            Pbusc=sqrt(1/Pbdp(p)^2+abs(bposu(p))^2);
            Pbdsc=sqrt(Pbusc^2-abs(bposd(p))^2);
            wfLu=wfLd*cNf(p)-[0;Ggain(1:p-1)]*sNf(p);
            Ggain=[0;Ggain(1:p-1)]*cNf(p)+wfLd*sNf(p);
            wbLu=wbLd*(Pbusc*Pbdp(p))-Ggain*(bposu(p)*Pbdp(p));
            Ggain=Ggain*(Pbusc*Pbdp(p))-wbLd*(bposu(p)*Pbdp(p));
            wfLd=[0;Ggaini(1:p-1)]*shNf(p)+wfLu*chNf(p);
            Ggaini=[0;Ggaini(1:p-1)]*chNf(p)+wfLu*shNf(p);
            wbLd=wbLu*(Pbdsc/Pbusc)+Ggaini*(bposd(p)/Pbusc);
            Ggaini=Ggaini*(Pbdsc/Pbusc)-wbLu*(bposd(p)/Pbusc);
            % Delay management --------------------------------------------
            bposdp=bposd;bposup=bposu;gammaup=gammau;gammadp=gammad;
            % Affine Projection -------------------------------------------
            ruu=ruu+u(r)*ur(1:L+p-2)'-ur(N)*ur(N+1:N+L+p-2)';
            ur=[u(r),ur(1:N+L+p-3)];
            k=(mod(r,L)+(mod(r,L)==0)*L);
            auxerr=d(r)-Y(L-k+1);
            err=auxerr-E(1:k+p-2)'*ruu(1:k+p-2);
            y(r)=d(r)-err;
            ei=[d(r)-y(r);(1-s)*ei(1:p-1)];
            Bi=wbLu*wbLu'*ei;
            Fi=wfLu*wfLu'*ei;
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
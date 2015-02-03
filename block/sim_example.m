% Block algorithms simulation example
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 03-Nov-2006$
%--------------------------------------------------------------------------
% Number of filter taps
f_len=64;
% Input signals length
s_len=2^16;
% Number of independent averaging
avr_len=5;
% AWGN Variance [dB]
Var=-40;
% TRUE if input signals are complex
cpl=0;
% Smoothing window length
smth_len=5;
% Shaping filter
a=0.95;
sNum=sqrt(1-abs(a)^2);
sDen=[ 1,-a ];
%--------------------------------------------------------------------------
% Filter taps initialize (FIR)
fNum=rand(1,f_len)+sqrt(-1)*(cpl~=0)*rand(1,f_len);
% Average output signal
avry=zeros(1,s_len);
% Average relative filter mismatch
av_mis=0;
% Independent averaging...
for k=1:avr_len,
    disp( sprintf('Iteration %.0f \\ %.0f',k,avr_len) );
    u=normrnd(0,1,1,s_len)+sqrt(-1)*(cpl~=0)*normrnd(0,1,1,s_len);
    e=normrnd(0,10^(0.05*Var),1,s_len)+...
        sqrt(-1)*(cpl~=0)*normrnd(0,10^(0.05*Var),1,s_len);
    % Input signal
    u = filter(sNum,sDen,u);
    % Desired signal
    d=filter(fNum,1,u)+e;
    % Filtering...
    % ---------------------------------------------------
    % ######### Choose algorithm unmarking one #########
    %[w,y]=cDFTBpNLMS( u,d,(0.1/64),0.90,32,f_len );
    %[w,y]=fcDFTBpNLMS( u,d,(0.1/64),0.90,32,f_len );
    %[w,y]=efcDFTBpNLMS( u,d,(0.1/64),0.90,32,8,f_len );
    %[w,y]=dlefcDFTBpNLMS( u,d,(0.1/64),0.90,32,8,f_len );
    %[w,y]=uDFTBpNLMS( u,d,(0.1/64),0.90,32,f_len );
    [w,y]=cDCTBpNLMS( u,d,0.1,0.90,32,f_len );
    %[w,y]=fcDCTBpNLMS( u,d,0.1,0.90,32,f_len );
    %[w,y]=uDCTBpNLMS( u,d,0.1,0.90,32,f_len );
    %[w,y]=cDHTBpNLMS( u,d,0.1,0.90,32,f_len );
    %[w,y]=fcDHTBpNLMS( u,d,0.1,0.90,32,f_len );
    %[w,y]=uDHTBpNLMS( u,d,0.1,0.90,32,f_len );
    % ---------------------------------------------------    
    avry=avry+abs(d-y).^2;
    av_mis=av_mis+norm(w-fNum.')^2/norm(fNum)^2;
end;
disp('Done!');
avry=avry/avr_len;
disp('Smoothing...');
avry=conv(avry,ones(1,max(smth_len,1)));
figure;hold on;
plot( 1:s_len-smth_len,10*log10(avry(max(smth_len,1):...
    length(avry)-smth_len)/max(smth_len,1)),'LineWidth',1,...
    'Color',[0.90,0.0,0.0] );
plot( 1:s_len-smth_len,10*log10(10^(0.1*Var)*ones(1,s_len-smth_len)),...
    'k','LineWidth',2 );
set(gca,'Box','on','FontName','Sylfaen','FontSize',12);
legend('Learning curve','Minimum Mean Square Error');
v=axis;axis( [v(1),s_len,v(3),v(4)] );
title(sprintf(strcat('Adaptive algorithm learning curve\n',...
'(%d independent averaging, smoothing window - %d samples)\n',...
'Average relative filter mismatch = %.2f [dB]'),...
avr_len,max(smth_len,1),10*log10(av_mis/avr_len)));
xlabel('n [Ts]');
ylabel('E|e(n)|^2 [dB]');
grid;
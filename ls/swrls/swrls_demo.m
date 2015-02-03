% Sliding Window recursive least squares algorithms demo
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 04-Nov-2006$
%--------------------------------------------------------------------------
% Number of filter taps
f_len=16;
% Input signals length
s_len=2^12;
% Number of independent averaging
avr_len=10;
% AWGN Variance [dB]
Var=-40;
% TRUE if input signals are complex
cpl=0;
% Smoothing window length
smth_len=4;
% Shaping filter
a=0.95;
sNum=sqrt(1-abs(a)^2);
sDen=[ 1,-a ];
%--------------------------------------------------------------------------
% Filter taps initialize (FIR)
fNum=rand(1,f_len)+sqrt(-1)*(cpl~=0)*rand(1,f_len);
% Average output signal
avry=zeros(9,s_len);
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
    disp('SWRLS');[w,y]=SWRLS( u,d,1e-3,128,f_len );
    avry(1,:)=avry(1,:)+abs(d-y).^2;
    disp('SWIQRRLS');[w,y]=SWIQRRLS( u,d,1e-3,128,f_len );
    avry(2,:)=avry(2,:)+abs(d-y).^2;
    disp('SWQRRLS');[w,y]=SWQRRLS( u,d,1e-3,128,f_len );
    avry(3,:)=avry(3,:)+abs(d-y).^2;
    disp('SWFARLS');[w,y]=SWFARLS( u,d,1e-3,128,f_len );
    avry(4,:)=avry(4,:)+abs(d-y).^2;
    disp('SWRLSLpos');[y]=SWRLSLpos( u,d,1e-3,128,f_len );
    avry(5,:)=avry(5,:)+abs(d-y).^2;
    disp('SWRLSLposf');[y]=SWRLSLposf( u,d,1e-3,128,f_len );
    avry(6,:)=avry(6,:)+abs(d-y).^2;
    disp('SWRLSLpri');[y]=SWRLSLpri( u,d,1e-3,128,f_len );
    avry(7,:)=avry(7,:)+abs(d-y).^2;
    disp('SWRLSLprif');[y]=SWRLSLprif( u,d,1e-3,128,f_len );
    avry(8,:)=avry(8,:)+abs(d-y).^2;
    disp('SWFARLSL');[y]=SWFARLSL( u,d,1e-3,128,f_len );
    avry(9,:)=avry(9,:)+abs(d-y).^2;
    % ---------------------------------------------------    
end;
disp('Done!');
ALG={'SWRLS';'SWIQRRLS';'SWQRRLS';'SWFARLS';'SWRLSLpos';'SWRLSLposf';...
    'SWRLSLpri';'SWRLSLprif';'SWFARLSL'};
avry=avry/avr_len;
disp('Smoothing...');
figure;
for i=1:9,
    subplot(3,3,i);
    tmp=conv(avry(i,:),ones(1,max(smth_len,1)));
    plot( 1:s_len-2*smth_len+1,10*log10(tmp(max(smth_len,1):length(avry)-...
        smth_len)/max(smth_len,1)),'LineWidth',1,...
    'Color',[0.90,0.0,0.0] );
    v=axis;axis( [v(1),s_len,-45,0] );
    title(ALG{i});
    grid;
end;
% Recursive least squares algorithms demo
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 04-Nov-2006$
%--------------------------------------------------------------------------
% Number of filter taps
f_len=16;
% Input signals length
s_len=2^10;
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
avry=zeros(14,s_len);
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
    disp('RLS');[w,y]=RLS( u,d,0.99,1e-3,f_len );
    avry(1,:)=avry(1,:)+abs(d-y).^2;
    disp('QRRLS');[w,y]=QRRLS( u,d,0.99,1e-3,f_len );
    avry(2,:)=avry(2,:)+abs(d-y).^2;
    disp('IQRRLS');[w,y]=IQRRLS( u,d,0.99,1e-3,f_len );
    avry(3,:)=avry(3,:)+abs(d-y).^2;
    disp('FARLS');[w,y]=FARLS( u,d,0.999,1e-3,f_len );
    avry(4,:)=avry(4,:)+abs(d-y).^2;
    disp('FTF');[w,y]=FTF( u,d,0.999,1e-3,f_len );
    avry(5,:)=avry(5,:)+abs(d-y).^2;
    disp('SFTF');[w,y]=SFTF( u,d,0.999,1e-3,f_len );
    avry(6,:)=avry(6,:)+abs(d-y).^2;
    disp('FAEST');[w,y]=FAEST( u,d,0.999,1e-3,f_len );
    avry(7,:)=avry(7,:)+abs(d-y).^2;
    disp('FKF');[w,y]=FKF( u,d,0.999,1e-3,f_len );
    avry(8,:)=avry(8,:)+abs(d-y).^2;
    disp('RLSLpos');[y]=RLSLpos( u,d,0.999,1e-3,f_len );
    avry(9,:)=avry(9,:)+abs(d-y).^2;
    disp('RLSLposf');[y]=RLSLposf( u,d,0.999,1e-3,f_len );
    avry(10,:)=avry(10,:)+abs(d-y).^2;
    disp('RLSLpri');[y]=RLSLpri( u,d,0.999,1e-3,f_len );
    avry(11,:)=avry(11,:)+abs(d-y).^2;
    disp('RLSLprif');[y]=RLSLprif( u,d,0.999,1e-3,f_len );
    avry(12,:)=avry(12,:)+abs(d-y).^2;
    disp('NRLSL');[y]=NRLSL( u,d,0.999,1e-3,f_len );
    avry(13,:)=avry(13,:)+abs(d-y).^2;
    disp('FARLSL');[y]=FARLSL( u,d,0.999,1e-3,f_len );
    avry(14,:)=avry(14,:)+abs(d-y).^2;
    % ---------------------------------------------------    
end;
disp('Done!');
ALG={'RLS';'QRRLS';'IQRRLS';'FARLS';'FTF';'SFTF';'FAEST';'FKF';...
    'RLSLpos';'RLSLposf';'RLSLpri';'RLSLprif';'NRLSL';'FARLSL'};
avry=avry/avr_len;
disp('Smoothing...');
figure;
for i=1:14,
    subplot(4,4,i);
    tmp=conv(avry(i,:),ones(1,max(smth_len,1)));
    plot( 1:s_len-2*smth_len+1,10*log10(tmp(max(smth_len,1):length(avry)-...
        smth_len)/max(smth_len,1)),'LineWidth',1,...
    'Color',[0.90,0.0,0.0] );
    v=axis;axis( [v(1),s_len,-45,0] );
    title(ALG{i});
    grid;
end;
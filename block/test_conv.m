% Block convolution algorithms example
% (C) Bartosz Zator (braton@gmail.com)
% $Date: 03-Nov-2006$
% -------------------------------------------------------------------------
% Input signal length
s_len = 4096;
% Long filter length (FIR)
g_len = 500;
% TRUE if input signal is complex
cpl=0;
% Plotting boundaries
l_bound=1000;
r_bound=2600;
% -------------------------------------------------------------------------
s=rand(1,s_len)+sqrt(-1)*(cpl~=0)*rand(1,s_len);
g=rand(1,g_len)+sqrt(-1)*(cpl~=0)*rand(1,g_len);
% Filtering...
disp( 'Matlab CONV routine:' );
tic; r=conv( s,g ); toc;
% ---------------------------------------------------
% ######### Choose algorithm unmarking one #########
%v=zeros(1,s_len+g_len-1);
disp( 'Block convolution via DFT:' );
tic; v=Bconv( g,s,97 ); toc;
%disp( 'Delayless block convolution via DFT:' );
%tic; v=dlBconv( g,s,97 ); toc;
%disp( 'Block convolution via DCT:' );
%tic; v=BDCTconv( g,s,97 ); toc;
%disp( 'Block convolution via DHT:' );
%tic; v=BDHTconv( g,s,97 ); toc;
% ---------------------------------------------------
if cpl, r=abs(r); v=abs(v); end;
figure;hold on;
plot( l_bound:r_bound,real(r(l_bound:r_bound)),'b','LineWidth',1 );
plot( l_bound:r_bound,real(v(l_bound:r_bound)),'g','LineWidth',1 );
set(gca,'Box','on','FontName','Sylfaen','FontSize',12);
legend('Matlab CONV routine','Block convolution algorithm');
title('Block convolution algorithms');
xlabel('n [Ts]');
ylabel('Convolved signals');
grid;
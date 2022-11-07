%Verification of Sampling Theorem
%% AUTHOR : ABHILASH SINGH 
clc                      %to clear the comand promt
clear
amplitude=input('Enter the Msg amplitude:\n'); % Amplitude of the signal
fm=input('Enter the Msg frequency(fm):\n');        % Frequency of the  msg signal
t=0:0.01:1;                    %range of the time axis with interval .00002 
xa1=amplitude*sin(2*pi*fm*t);     %our signal
subplot(2,2,1);                  %creating four plot,we are working in the 1st plot
plot(t,xa1);                     % ploting our signal using plot funtion
title('Continuous sinusoidal signal');%giving title to the 1st plot
xlabel('t');                     %labeling the independent axis as 'time' of 1st plot 
ylabel('x(t)');                  %labeling the dependent axis as 'x(t)' of 1st plot
fs1=input('Enter the sampling frequency greater than 2xfm :\n');
n=0:1/fs1:1;                     %sampling time(i.e sampling frequency of 1000Hz) which is much much greater than nyquist rate
xa2=amplitude*sin(2*pi*n*fm);            %Sampling our signal
subplot(2,2,2);                  %working in the  2nd plot out of the four
stem(n,xa2);                     %ploting the signal in discrete mode 
title('Above niquist rate');     %giving title of the  2nd plot
xlabel('n');                     %labeling the independent axis as 'n' of 2nd plot
ylabel('x(n)');                  %labeling the dependent axis as 'x[n]' of 2nd plot
fs2=input('enter the sampling frequency less than  2xfm :\n');
n1=0:1/fs2:1;
xa3=amplitude*sin(2*pi*n1*fm);
subplot(2,2,3);                  %working in the 3rd plot
stem(n1,xa3);                     %ploting the signal in discrete mode
title('      Under sampled case');%giving title to the 3rd plot
xlabel('n');                     %labeling the independent axis as 'n' of 3rd plot 
ylabel('x(n)');                  %labeling the dependent axis as 'x[n]' of 3rd plot
%ts2=1/50;                        %sampling time(i.e sampling frequency of 50Hz) which is equal to the  nyquist rate
fs3=input('Enter the sampling frequency equal 2xfm :\n');
n2=0:1/fs3:1
clc;
xa4=amplitude*sin(2*pi*n2*fm);         %sampling our signal
subplot(2,2,4);                  %working in the 4th plot
stem(n2,xa4);                     %ploting the signal in discrete mode
title('At Niquist rate');        %giving title to the 4th plot
xlabel('n');                     %labeling the independent axis as 'n' of 4th plot
ylabel('x(n)');                  %labeling the dependent axis as 'x[n]' of 4th plot
%**************************************************************************
%**************************************************************************A
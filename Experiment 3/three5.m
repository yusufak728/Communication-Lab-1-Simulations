fs = 1000; 
fc = 200;  
t = (0:1/fs:0.2)';
x = 5*sin(2*pi*30*t);
fDev = 5;
y = fmmod(x,fc,fs,fDev);
plot(t,x,'black',t,y,'b')
xlabel('Time (s)');
ylabel('Amplitude');
legend('Original Signal','Frequency Modulated Signal');
subplot(2,1,1);
plot(t,x)
title('Original Signal');
xlabel('time(t)');
ylabel('Amplitude');


subplot(2,1,2); 
plot(t,y)
title('Frequency Modulated Signal');
xlabel('time(t)');
ylabel('Amplitude');
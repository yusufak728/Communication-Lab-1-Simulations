% A simple sampling and reconstruction model for students
% beginners of Digital Signal Processing
% by Mukhtar Hussain (Email: mukhtarhussain@ciitlahore.edu.pk)

% f - The frequency of analog sinosoid signal
% F - Sampling Rate
% qbits - Number of Quantizations bits
% A - Amplitude of sinusoid signal
% L - Number of quantization levels based on qbits
% I - Quantization Interval
% sim_time - Simultaion Time
% span - x-axis range of frequency plot 1 & 3 (spectrum scope 1 & 3)
% span1 - x-axis range of frequency plot 2 (spectrum scope 2)
% NFFT - Number of FFT points

clc;
clear;
close all;

f = input('Enter the frequency of signal = ');
F = input('Enter the sampling frequency = ');
A = input('Enter max amplitude of signal = ');
qbits = input('Enter the number of quantization bits = ');
fc = input('Enter the lowpass filter cutoff frequency = ');

L = 2^qbits;
I = 2*A/(L-1); 

% Settings for Spectrum Scope
span = 8*F;
span1 = 8*F;
NFFT = 256;

% To run simulink model
t = 1/f;
sim_time = 10*t; 
sim('sampling.slx');
%sampling_demo   Interactively demonstrates Nyquist's Sampling Theorem.
%--------------------------------------------------------------------------
%   sampling_demo demonstrates Nyquist's Sampling Theorem, by sampling a
%   continuous-time sinusoidal signal of a frequency f = 50 Hz to 3 kHz,
%   with a fixed sampling frequency fs = 2 kHz. Therefore, only signals
%   with frequencies f <= fs/2 = 1 kHz can be faithfully reconstructed by
%   their samples, while those with frequencies f > 1 kHz will exhibit
%   aliasing effects (i.e., reconstructing the original signal will
%   generate new components that were not part of the original signal).
%
%   The generated graph shows both the original signal (in red) and the
%   reconstructed signal (in blue). The reconstructed signal is linearly
%   interpolated and shown together with the samples (circle markers). The
%   user can interactively change the original signal's frequency and
%   phase, and observe how this affects the reconstructed signal, as the
%   frequency approaches and crosses the fs/2 boundary.
%
%   In order to have a better appreciation of the aliasing effects on the
%   signal being sampled, the user can listen to both the original and the
%   reconstructed signals and compare the corresponding sounds. This is
%   only possible on sound-capable computers.
%
%   The simulation can be interactively controlled in real-time using the
%   following keys:
%
%   [Q]/[A] increases/decreases the original signal's frequency (coarse).
%   [W]/[S] increases/decreases the original signal's frequency (fine).
%   [E]/[D] increases/decreases the original signal's phase.
%   [I]     plays a one-second sound of the original signal.
%   [O]     plays a one-second sound of the reconstructed signal.
%   [R]     resets all parameters to their initial values.
%   [ESC]   exits the script.

%   Version:    1.50
%   Programmer: Costas Vlachos
%   Date:       07-Dec-2016
%   Updated:    14-Jan-2017

% Initialize frequency parameters
fs_in = 50000;
fs_out = 2000;
f_min = 50;
f_max = 3000;
f_factor_coarse = 0.002;
f_factor_fine = 0.00005;
f = f_min;

% Initialize phase parameters
phi_min = 0;
phi_max = 360;
phi_step = 1;
phi = phi_min;

% Initialize time parameters
t_sound = 1;
t_width = 1/f_min;
t_long_in = linspace(0,t_sound,fs_in*t_sound+1);
t_in = t_long_in(1:fs_in*t_width+1);
t_long_out = linspace(0,t_sound,fs_out*t_sound+1);
t_out = t_long_out(1:fs_out*t_width+1);

% Create figure & event handlers
figure('Name','Interactive Demo of Nyquist''s Sampling Theorem', ...
    'KeyPressFcn',@(~,event) setappdata(0,'c',event.Key), ...
    'KeyReleaseFcn',@(~,~) setappdata(0,'c',0), ...
    'ResizeFcn',@(~,~) setappdata(0,'s',get(gcf,'Position')));

% Infinite loop
while 1

    % Compute input/output signals
    s_long_in = sin(2*pi*f*t_long_in+phi*pi/180);
    s_in = s_long_in(1:fs_in*t_width+1);
    s_long_out = sin(2*pi*f*t_long_out+phi*pi/180);
    s_out = s_long_out(1:fs_out*t_width+1);

    % Plot input/output signals & sample markers
    plot(1000*t_in,s_in,'r','LineWidth',3);
    hold on
    plot([0 1000*t_width],[0 0],'k');
    plot(1000*t_out,s_out,'b','LineWidth',3);
    plot(1000*t_out,s_out,'ko','LineWidth',3,'MarkerSize',10);
    hold off

    % Adjust figure aesthetics
    s = getappdata(0,'s');
    axis([0 1000*t_width -1.1 1.1])
    title(['{\itf} = ' num2str(f,'%.1f') ' Hz     \phi = ' ...
        num2str(phi) '^\circ     {\itf_s} = ' num2str(fs_out) ...
        ' Hz     {\itf_s} / 2 = ' num2str(fs_out/2) ' Hz'])
    xlabel('Time (ms)')
    ylabel('Signal value')
    set(gca,'FontSize',s(3)/50)
    drawnow

    % Control simulation parameters on keyboard input
    c = getappdata(0,'c');
    if c == 'q', f = min(f*(1+f_factor_coarse),f_max); end
    if c == 'a', f = max(f*(1-f_factor_coarse),f_min); end
    if c == 'w', f = min(f*(1+f_factor_fine),f_max); end
    if c == 's', f = max(f*(1-f_factor_fine),f_min); end
    if c == 'e', phi = min(phi+phi_step,phi_max); end
    if c == 'd', phi = max(phi-phi_step,phi_min); end
    if c == 'r', f = f_min; phi = phi_min; end
    if c == 'i', try sound(s_long_in,fs_in), catch, end, end
    if c == 'o', try sound(s_long_out,fs_out), catch, end, end
    if strcmp(c,'escape'), break; end

end

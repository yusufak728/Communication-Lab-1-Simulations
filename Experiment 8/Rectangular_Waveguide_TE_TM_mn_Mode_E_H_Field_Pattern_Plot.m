% E-H Field Pattern plot for Rectangular waveguide for TEmn and TMmn mode
% SEP 20, 2018
% Author: AJEET KUMAR

clc;
close all;

% Waveguide dimensions
a = 2.286;  % Length in cm in x-direction
b = a/2;    % Length in cm in y-direction
f = 45*10.^9;   % Frequency of operation 45GHz
c = 3*10.^8;    % Velocity of light 
% m = 1;    % Mode number in X-Direction
% n = 0;    % Mode number in Y-Direction
choice = input('Enter choice: 1 for TE and 2 for TM:    ');
if choice == 1
     m = input('Enter mode value m:');
     n = input('Enter mode value n:');
    
elseif choice == 2
     m = input('Enter mode value m:');
     n = input('Enter mode value n:');    
else
    sprintf('Alert!!! Wrong choice!!!')
    
end




Amn = 1;    % Particular mode Constant
% A10 = 1;    % for example

% Wave propagation in Z-Direction
%********************************%
fc = c*100/2*sqrt((m/a).^2+(n/b).^2);    % Cutoff frequency calculation in GHz
% lambda = 2*a;                 %for TE10 mode
lambda = c*100/fc;              % Wavelength in cm
epsilon = 8.8540e-12;           % Permittivity constant
epsilon_r = 1;                  % Relative Permittivity constant

mu1 = 4*pi*10e-7;               % Permeability constant
mu1_r = 1;                      % Relative Permeability constant
omega = 2*pi*f;                 % Frequency of operation in rad/s
M = 40;                         % Number of points to be poltted

beta = omega*(sqrt(mu1*epsilon));  %Propagation constant
Bx = m*pi/a;    %Beta(x)
By = n*pi/b;    %Beta(y)
Bc = sqrt(Bx.^2+By.^2); %Beta(c), cutoff wavenumber
Bz = sqrt(beta.^2-Bc.^2);



if choice ==1
    if m == 0 && n == 0
        fprintf(['TE_',num2str(m),num2str(n), ' mode doesnot exist']);
    elseif fc>f
        fprintf(['TE_',num2str(m),num2str(n), ' mode cutoff frequency exceeds frequency of operation; hence mode does not porpagate\n']);
        sprintf('The frequency of operation is up to: %0.5g',f)
        sprintf('The cutoff frequency is: %0.5g',fc)
        
    else
        sprintf('The frequency of operation is up to: %0.5g',f)
        sprintf('The cutoff frequency is: %0.5g',fc)
% Front View
z = 0;
x = linspace(0,a,M);
y = linspace(0,b,M);
[x,y] = meshgrid(x,y);
% z = linspace(0,2*lambda,M);

%Field Expression for TEmn
% Ex = Amn*(By/epsilon)*cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
% Ex = Amn*(By/epsilon)*cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
Ex = cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);       
% Ey = -Amn*(Bx/epsilon)*sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
Ez = 0;

% Hx = Amn*(Bx*Bz/(omega*mu1*epsilon))*sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
% Hy = Amn*(Bx*Bz/(omega*mu1*epsilon))*cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
% Hz = -1i*Amn*(Bc.^2/(omega*mu1*epsilon))*cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
Hz = -cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);


figure();
quiver(x,y,real(Ex),real(Ey));
title(['Plot of front view for TE_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
figure();
quiver(x,y,real(Hx),real(Hy));
title(['Plot of front view for TE_',num2str(m),'_',num2str(n),' H-Field']);
legend('H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
figure();
quiver(x,y,real(Ex),real(Ey));
hold on
quiver(x,y,real(Hx),real(Hy));
grid on
title(['Plot of front view for TE_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');

% Top View for TEmn
y = b;      % Position of x-z plane
x = linspace(0,a,M);
% y = linspace(0,b,M);
z = linspace(0,lambda,M);
[x,z] = meshgrid(x,z);      % Create Mesh grid in x-z

% Field Expression for TEmn
% Ex = Amn*(By/epsilon)*cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
Ex = cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);    
Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
% Ez = 0;
Ez = zeros(size(real(Ey)));

Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1j*Bz*z);
% Hx = A10*(Bz/(omega*mu1*epsilon))*pi/a.*sin(pi.*x./a).*exp(-j*Bz*z);
Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-1j*Bz*z);
Hz = -cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-1j*Bz*z);


figure();
quiver(z,x,real(Ez),real(Ex));
title(['Plot of Top view for TE_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');
figure();
quiver(z,x,real(Hz),real(Hx));
title(['Plot of Top view for TE_',num2str(m),'_',num2str(n),' H-Field']);
legend('H-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');
figure();
quiver(z,x,real(Ez),real(Ex));
hold on
quiver(z,x,real(Hz),real(Hx));
grid on
title(['Plot of TOP view of E-H for TE_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');

% Side View for TEmn
x = a/2;
% x = linspace(0,a,M);
y = linspace(0,b,M);
z = linspace(0,2*lambda,M);
[y,z] = meshgrid(y,z);

% Field Expressions for TEmn
Ex = cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
Ez = 0;
Ez = zeros(size(real(Ey)));

Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
Hz = -cos(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);


figure();
quiver(z,y,real(Ez),real(Ey));
title(['Plot of Side view for TE_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
figure();
quiver(z,y,real(Hz),real(Hy));
title(['Plot of Side view for TE_',num2str(m),'_',num2str(n),' H-Field']);
legend('E-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
figure();
quiver(z,y,real(Ez),real(Ey));
hold on
quiver(z,y,real(Hz),real(Hy));
grid on
title(['Plot of Side view of E-H for TE_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
    end

elseif choice == 2
    
     if m == 0 || n == 0
        fprintf(['TM_',num2str(m),num2str(n), ' mode doesnot exist']);
    elseif fc>f
        fprintf(['TM_',num2str(m),num2str(n), ' mode cutoff frequency exceeds frequency of operation; hence mode does not porpagate\n']);
        sprintf('The frequency of operation is up to: %0.5g',f)
        sprintf('The cutoff frequency is: %0.5g',fc)
     else
        sprintf('The frequency of operation is up to: %0.5g',f)
        sprintf('The cutoff frequency is: %0.5g',fc)
% Field Pattern plot for Rectangular wave guide for TMmn mode
%TM_mn mode field expressions

% Front View
x = linspace(0,a,M);
y = linspace(0,b,M);
% z = linspace(0,2*lambda,M);
z = 0;
[x,y] = meshgrid(x,y);

% % % Field Expressions for TMmn
% % Ex = -cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
% % Ez = -sin(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % 
% % Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
% % Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
% % % Hz = 0;
% % Hz = zeros(size(real(Hy)));
tmequation();
%Plot of TMmn E-Field view
figure();
quiver(x,y,real(Ex),real(Ey));
title(['Plot of front view for TM_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
%Plot of TMmn H-Field view
figure();
quiver(x,y,real(Hx),real(Hy));
title(['Plot of front view for TM_',num2str(m),'_',num2str(n),' H-Field']);
legend('H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');
%Plot of TMmn E-Field and H-Field view
figure();
quiver(x,y,real(Ex),real(Ey));
hold on
quiver(x,y,real(Hx),real(Hy));
grid on
title(['Plot of front view for TM_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
xlabel('x-dimension 0 to a');
ylabel('y-dimension 0 to b=a/2');

% Top View
y = b;      %Position of view
x = linspace(0,a,M);
% y = linspace(0,b,M);
z = linspace(0,lambda,M);
[x,z] = meshgrid(x,z);

% % %Field expression for TMmn
% % Ex = -cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
% % Ez = -sin(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % 
% % Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
% % Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
% % % Hz = 0;
% % Hz = zeros(size(real(Hy)));

tmequation();

figure();
quiver(z,x,real(Ez),real(Ex));
title(['Plot of Top view for TM_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');
figure();
quiver(z,x,real(Hz),real(Hx));
title(['Plot of Top view for TM_',num2str(m),'_',num2str(n),' H-Field']);
legend('H-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');
figure();
quiver(z,x,real(Ez),real(Ex));
hold on
quiver(z,x,real(Hz),real(Hx));
grid on
title(['Plot of TOP view of E-H for TM_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
ylabel('x-dimension 0 to a');
xlabel('z-direction');



% Side View
x = a/2;
% x = linspace(0,a,M);
y = linspace(0,b,M);
z = linspace(0,2*lambda,M);
[y,z] = meshgrid(y,z);

% % %Field Expression for TMmn
% % Ex = -cos(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % Ey = -sin(Bx.*x).*cos(By.*y).*exp(-1i*Bz*z);
% % Ez = -sin(Bx.*x).*sin(By.*y).*exp(-1i*Bz*z);
% % 
% % Hx = sin(m*pi.*x./a).*cos(n*pi.*y./b).*exp(-j*Bz*z);
% % Hy = cos(m*pi.*x./a).*sin(n*pi.*y./b).*exp(-j*Bz*z);
% % % Hz = 0;
% % Hz = zeros(size(real(Hy)));
tmequation();

figure();
quiver(z,y,real(Ez),real(Ey));
title(['Plot of Side view for TM_',num2str(m),'_',num2str(n),' E-Field']);
legend('E-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
figure();
% quiver(y,z,real(Hy),real(Hz));
quiver(z,y,real(Hz),real(Hy));
title(['Plot of Side view for TM_',num2str(m),'_',num2str(n),' H-Field']);
legend('E-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
figure();
quiver(z,y,real(Ez),real(Ey));
hold on
quiver(z,y,real(Hz),real(Hy));
grid on
title(['Plot of Side view of E-H for TM_',num2str(m),'_',num2str(n)]);
legend('E-Field','H-Field');
ylabel('y-dimension 0 to b');
xlabel('z-direction');
     end
else
        
        sprintf('Alert!!! Something went wrong, try again!!!');
end

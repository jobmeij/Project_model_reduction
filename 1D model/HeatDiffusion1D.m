%% 1D heat diffusion
clear all; close all; clc;

% Model parameters
L = 0.5;        % Length beam [m]
Tamb = 309;     % Ambient temperature [K]
u_x = 0.4;      % Heat input [m]
rho = 1100;
c = 3890;
kappa = 0.31;
resolution = 100;   % number of points 
x = linspace(0,L,resolution+1)';   % length vector [m]
Nk = 20;
r = 10;

% Simulation settings
tstep = 0.1;    % sec
tend = 10;      % sec
t = linspace(0,tend,(tend/tstep)+1);

% Compute Fourier basis: phi
phi = zeros(Nk,length(x));
phi2 = zeros(Nk,length(x));
for k = 1:Nk
    for i = 1:length(x) % used i to index x
        % First
        if (k == 1)
            phi(k,i) = 1/L;
        else
            phi(k,i) = sqrt(2/L)*cos((pi*k*x(i))/L);
        end
    end
    
    % Computing phi2(k,x)
    phi2(k,:) = ((-pi^2*k^2)/(L^2))*phi(k,:);
end

if (false)
    % Plotting phi
    figure(1)
    plot(x,phi)
    grid on
    title('Phi over distance x')
    xlabel('x')
    ylabel('phi')
    
    figure(2)
    plot(x,phi2)
    grid on
    title('Phi2 over distance x')
    xlabel('x')
    ylabel('phi2')
end

% phi_inner = phi2'*phi;

% Computing A matrix
A1 = zeros(r);
A2 = zeros(r);
for i = 1:r 
    % Calculate A based on slideset 8 sheet 17
    A1(i,i) = ((-i^2*pi^2)/L^2);

    % Calculate A based on slideset 8 sheet 20
    for j = 1:r
       A2(i,j) = phi(i,:)*phi2(j,:)';
    end
end
% Multiply matrix with constants
A1 = (kappa/(rho*c))*A1;
A2 = (kappa/(rho*c))*A2*tstep;

% Combute B matrix (s'*phi(r,:)) where u(x,t)=s(x)v(t)
B = zeros(r,1);
s = x;
for i = 1:r
    B(r,1) = s'*phi(r,:)';
end
% TBD: Multiply B with dx?? slideset 8 sheet 17: no, sheet 20: yes
B = B*(1/(rho*c));
% B = B*(1/(rho*c))*tstep;

% Compute output



% Plot output




%% 5LMA0 - Project heat diffusion

%% Simulation
clear all; close all; clc;

% parameter:                unit:
Lx = 0.2;                   % m
Ly = 0.15;                  % m
rho = [1100, 1000];         % kg/m3
c = [3890,3350];            % J/(kg K)
kappa = [0.31, 0.31];       % W/(m K)
X1 = Lx/4;                  % m
Y1 = Ly/2;                  % m
X2 = 3*Lx/4;                % m
Y2 = Ly/2;                  % m
W = 0.05;                   % m
Tamb = 309;                 % K
Tend = 1;                   % End time of simulation

rx = 5;
ry = 5;

Nx = 50;
Ny = 50;
M = 100;

dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
dt = Tend/(M-1);

x_vect = 0:dx:Lx;
y_vect = 0:dy:Ly;
t_vect = 0:dt:Tend;

% Model is assumed to be isotropic, so:
K = kappa*eye(2);

% Computing phi x
phix = zeros(Nx,rx);
phix(:,1) = 1/sqrt(Lx);
for k = 2:rx
    phix(:,k) = sqrt(2/Lx)*cos(((k-1)*pi*x_vect)/(Lx));
end

% Compute phi y
phiy = zeros(Ny,ry);
phiy(:,1) = 1/sqrt(Ly);
for k = 2:ry
    phiy(:,k) = sqrt(2/Ly)*cos(((k-1)*pi*y_vect)/(Ly));
end

% Compute phi k,l
for x = 1:length(phix)
    for y = 1:length(phiy)
        for k = 1:rx
            for l = 1:ry
                phi_kl(x,y,k,l) = phix(x,k)*phiy(y,l);
            end
        end
    end    
end

% Determine a_k,l


% Exercise 1
% the model is called homogeneous if rho, c and K do not depend on (x,y)
% the model is called isotropic if K = kappa*I with I = eye(2). In that
% case kappa is called the thermal conductivity and has [W/m K] as unit.


%% Plotting
% Create grid
n = 100;    % amount of points
x = 0:Lx/n:Lx;
y = 0:Ly/n:Ly;
[X,Y] = meshgrid(x,y);

figure(1)
Temperature = sin(X).^2+cos(Y).^2;  % random, replace with result
Temperature(1:101,1:101) = Tamb;
surf(X,Y,Temperature)
colorbar
title('Temperature of plate [K]')
view(0, 90)
xlabel('x axis')
ylabel('y axis')


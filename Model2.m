clear all;
close all hidden;
clc;

Lx = 0.2;
Ly = 0.15;
rho = 1100;
c = 3890;
k = 0.31;
p1 = [(1*Lx)/4,Ly/2];
p2 = [(3*Lx)/4,Ly/2];
W = 0.05;
Tamb = 309;

Tend = 1;

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

phix = zeros(Nx,rx);
phix(:,1) = 1/sqrt(Lx);
for k = 2:rx
    phix(:,k) = sqrt(2/Lx)*cos(((k-1)*pi*x_vect)/(Lx));
end

phiy = zeros(Ny,ry);
phiy(:,1) = 1/sqrt(Ly);
for k = 2:ry
    phiy(:,k) = sqrt(2/Ly)*cos(((k-1)*pi*y_vect)/(Ly));
end

phi = phix*phiy;



%% Plotting
figure(1)

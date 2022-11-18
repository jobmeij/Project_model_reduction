clear all;
close all hidden;
clc;

% MODEL PARAMETERS
disp("Set model parameters");
model.Lx = 0.2;
model.rho = 1100;
model.c = 3890;
model.k = 0.31;
model.p1 = (1*model.Lx)/4;
model.p2 = (3*model.Lx)/4;
model.W = 0.05;
model.Tamb = 309;

% SIMULATION PARAMETERS
disp("Set simulation parameters");
par.Tend = 600;
par.rx = 200;
par.Nx = 1000;
par.M = 1800;
par.dx = model.Lx/(par.Nx-1);
par.dt = par.Tend/(par.M-1);

vect.x = (0:par.dx:model.Lx)';
vect.t = (0:par.dt:par.Tend)';

% INPUT PARAMETERS
disp("Set input parameters");
u1.Starttime = 60;
u1.Duration = 300;
u1.Amplitude = 50000;
u1.InitTemp = 0;
u2.Starttime = 0;
u2.Duration = 0;
u2.Amplitude = 0;
u2.InitTemp = 4;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,vect,u1,u2);

% BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate basis functions and gradients");
phi = CalculateBasis(model,par,vect);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,U,10,1);

% POD
T = T - model.Tamb;
%%
[POD.U,POD.S,POD.Y] = svd(T*T');
% [POD.U,POD.S,POD.Y] = svd(T);
figure(2);
plot(log10(diag(POD.S)),'.r');
ylabel("Log10(Hankel singular values)");
par.rx = input('Enter trunctation value: ');

phi.x = POD.U(:,1:par.rx)/sqrt(par.dx);

phi.xdd = zeros(par.Nx,par.rx); % Clear phi.xdd
for i = 1:par.rx;
    phi.xdd(:,i) = gradient(gradient(phi.x(:,i),par.dx),par.dx);
end

sys = SolveEquation(model,par,init,U,phi);

T_POD = CalculateOutput(model,phi,sys);

PlotTemperature(T_POD,model,par,vect,init,U,10,3);

%% EVD


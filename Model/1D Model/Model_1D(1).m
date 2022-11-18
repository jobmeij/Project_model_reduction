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
par.rx = 40;
par.Nx = 200;
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
input = InputFunction(model,vect,u1,u2);

% BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate basis functions and gradients");
phi = CalculateBasis(model,par,vect);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,input,phi);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,input,10);

%% POD

par.r = 25;
[POD.U,POD.S,POD.Y] = svd(T'*T);

POD.Ur = POD.U(:,1:par.r);

POD.Vr = POD.Ur;
POD.A = POD.Vr'*sys.A*POD.Ur;
% Sys.A = 40x40, POD.Ur=Pod.Vr=1800x1
POD.B = POD.Vr'*sys.B;

POD.a = zeros(size(POD.A,1),par.M);
POD.a_dot = zeros(size(POD.A,1),par.M);
for t = 1:par.M-1
    POD.a_dot(:,t) = POD.A*POD.a(:,t) + POD.B*input.g(t);
    POD.a(:,t+1) = POD.a(:,t+1) + POD.a_dot(:,t)*par.dt; 
end

% OUTPUT TEMPERATURE PROFILE
T = CalculateOutput(model,par,phi,POD);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,input,1);


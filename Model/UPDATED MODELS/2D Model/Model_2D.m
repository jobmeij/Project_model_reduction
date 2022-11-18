clear all;
close all hidden;
clc;

% MODEL PARAMETERS
disp("Set model parameters");
model.Lx = 0.2;
model.Ly = 0.15;
model.rho = [1100 1000];
model.c = [3890 3350];
model.k = 0.31;
model.p1 = [(1*model.Lx)/4,model.Ly/2];
model.p2 = [(3*model.Lx)/4,model.Ly/2];
model.W = 0.05;
model.p3 = [model.Lx/2,model.Ly/2];
model.lx = 0.03;
model.ly = 0.08;
model.Tamb = 309;
model.Homogeneous = true;
model.Analytic = true;
model.ShowProfile = false;

% SIMULATION PARAMETERS
disp("Set simulation parameters");
par.Tend = 1800;
par.rx = 20;
par.ry = 20;
par.r = par.rx*par.ry;
par.rPOD = 5;
par.Nx = 100;
par.Ny = 100;
par.M = 3600;
par.dx = model.Lx/(par.Nx-1);
par.dy = model.Ly/(par.Ny-1);
par.dt = par.Tend/(par.M-1);

vect.x = (0:par.dx:model.Lx)';
vect.y = (0:par.dy:model.Ly)';
vect.t = (0:par.dt:par.Tend)';

% INPUT PARAMETERS
disp("Set input parameters");
u1.Starttime = 60;
u1.Duration = 1500;
u1.Amplitude = 0;
u1.InitTemp = 5;
u1.InitWidth = 0.02;
u2.Starttime = 60;
u2.Duration = 300;
u2.Amplitude = 0;
u2.InitTemp = 0;
u2.InitWidth = 0.02;
p3.InitTemp = 0;
p3.InitWidth = 0.02;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,par,vect,u1,u2,p3);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,par,vect,u1,u2);

% BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate basis functions and gradients");
phi = CalculateBasis(model,par,vect);

% CALCULATE MATERIAL PROPERTIES
disp("Calculate material properties");
e = CalculateMaterialProperties(model,par,vect,phi);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi,e);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,U,1000,1);

% BASIS FUNCTION'S PLOT
PlotHarmonics(phi,par,vect,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% COMPUTE POD BASIS
disp("Calculate POD basis functions and gradients");
Tsnap = T - model.Tamb;
[phi, par] = CalculatePODbasis(par,Tsnap);

% INPUT PARAMETERS
disp("Set input parameters");
% u1.Starttime = 0;
% u1.Duration = 0;
% u1.Amplitude = 0;
% u1.InitTemp = 0;
% u1.InitWidth = 0.01;
% u2.Starttime = 0;
% u2.Duration = 0;
% u2.Amplitude = 0;
% u2.InitTemp = 10;
% u2.InitWidth = 0.01;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,par,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,par,vect,u1,u2);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi,e);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T_POD = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T_POD,model,par,vect,init,U,100,2);
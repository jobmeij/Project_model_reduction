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
model.Tamb = 309;   %~36 degrees Celsius

% SIMULATION PARAMETERS
disp("Set simulation parameters");
par.Tend = 600;
par.rx = 50;        % # harmonics in x direction (K)
par.ry = par.rx;    % # harmonics in y direction (L)
par.r = par.rx*par.ry;
par.rPOD = 3;
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
u1.Starttime = 30;
u1.Duration = 330;      % Heat for 5 minutes
u1.Amplitude = 60000;   % Reaching about 40.5 degC
u1.InitTemp = -5;
u1.InitWidth = 0.03;
u2.Starttime = u1.Starttime;
u2.Duration = u1.Duration;
u2.Amplitude = u1.Amplitude;
u2.InitTemp = 5;
u2.InitWidth = 0.03;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,par,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,par,vect,u1,u2);

% BASIS FUNCTION'S AND ITS GRADIENTS
disp("Calculate basis functions and gradients");
phi = CalculateBasis(model,par,vect);

% CALCULATE MATERIAL PROPERTIES
disp("Calculate material properties");
e = CalculateMaterialProperties(model,par,vect,phi,false);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi,e,false,true);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,U,100,1);

% BASIS FUNCTION'S PLOT
PlotHarmonics(phi,par,vect,false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% COMPUTE POD BASIS
disp("Calculate POD basis functions and gradients");
Tsnap = T - model.Tamb;
[phi, par] = CalculatePODbasis(par,Tsnap);

% changing input parameter
model.W = 0.075

% INPUT PARAMETERS
disp("Set input parameters");
u1.Starttime = 30;
u1.Duration = 330;
% u1.Amplitude = 0;
% u1.InitTemp = 5;
% u1.InitWidth = 0.03;
% u2.Starttime = 0;
% u2.Duration = 0;
% u2.Amplitude = 0;
u2.InitTemp = 5;
u2.InitWidth = 0.03;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,par,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,par,vect,u1,u2);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi,e,true,false);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T_POD = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T_POD,model,par,vect,init,U,100,2);
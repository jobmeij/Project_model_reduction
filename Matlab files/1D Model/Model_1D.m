clear all;
close all hidden;
clc;

% MODEL PARAMETERS
disp("Set model parameters");
model.Lx = 0.2;
model.rho = [1100 1000];
%model.rho = [1100 200];
model.c = [3890 3350];
%model.c = [3890 3000];
model.k = 0.31;
model.p1 = (1*model.Lx)/4;
model.p2 = (3*model.Lx)/4;
model.W = 0.05;
model.p3 = 4*model.Lx/8;
model.lx = 0.03;
model.Tamb = 309;
model.Homogeneous = true;
model.Analytic = true;

% SIMULATION PARAMETERS
disp("Set simulation parameters");
par.Tend = 1800;
par.rx = 50;
par.rPOD = 5;
par.Nx = 1000;
par.M = 7200;
par.dx = model.Lx/(par.Nx-1);
par.dt = par.Tend/(par.M-1);

vect.x = (0:par.dx:model.Lx)';
vect.t = (0:par.dt:par.Tend)';

% INPUT PARAMETERS
disp("Set input parameters");
u1.Starttime = 60;
u1.Duration = 300;
u1.Amplitude = 0;
u1.InitTemp = 0;
u1.InitWidth = 0.01;
u2.Starttime = 0;
u2.Duration = 0;
u2.Amplitude = 0;
u2.InitTemp = 5;
u2.InitWidth = 0.01;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,vect,u1,u2);

% BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate basis functions and gradients");
phi = CalculateBasis(model,par,vect);

X1 = ceil((model.p3-model.lx/2)/par.dx);
X2 = ceil((model.p3+model.lx/2)/par.dx);
F = ones(par.Nx,1)*model.rho(1)*model.c(1);
F(X1:X2) = model.rho(2)*model.c(2);
e = zeros(1,par.rx);
S = zeros(1,par.Nx);
for i = 1:par.rx
    e(i) = phi.x(:,i)'*F*par.dx;
    S = S + phi.x(:,i)'*e(i);
end

figure(2);
plot(vect.x, S, 'b');
hold on;
plot(vect.x, F, 'r');
grid minor;

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sysP = SolveEquation(model,par,init,U,phi,e);
model.Homogeneous = false;
sys = SolveEquation(model,par,init,U,phi,e);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T(:,:,1) = CalculateOutput(model,phi,sys);
T(:,:,2) = CalculateOutput(model,phi,sysP);

% TEMPERATURE PLOT
PlotTemperature(T,model,par,vect,init,U,300,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% POD BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate POD basis functions and gradients");
Tsnap = T - model.Tamb;
[phi,par] = CalculatePODbasis(par, Tsnap);

% INPUT PARAMETERS
disp("Set input parameters");
% u1.Starttime = 600;
% u1.Duration = 600;
% u1.Amplitude = 0;
% u1.InitTemp = 0;
% u1.InitWidth = 0.01;
% u2.Starttime = 0;
% u2.Duration = 0;
% u2.Amplitude = 0;
% u2.InitTemp = 4;
% u2.InitWidth = 0.01;

% INITIAL PROFILE
disp("Calculate initial profile");
init = InitialProfile(model,vect,u1,u2);

% INPUT FUNCTION
disp("Calculate input function");
U = InputFunction(model,vect,u1,u2);

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct and solve linear equations");
sys = SolveEquation(model,par,init,U,phi,e);

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T_POD = CalculateOutput(model,phi,sys);

% TEMPERATURE PLOT
PlotTemperature(T_POD,model,par,vect,init,U,1000,2);


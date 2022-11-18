function init = InitialProfile(model,par,vect,u1,u2)
% mu_Tx = model.Lx*0.2;
% sigma_Tx = model.Lx*0.05;
% Tx = (exp((-(vect.x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
% Tx = init.Tx/max(Tx);
% mu_Ty = model.Ly*0.2;
% sigma_Ty = model.Ly*0.1;
% Ty = (exp((-(vect.y-mu_Ty).^2)/(2*sigma_Ty^2))/(sigma_Ty*sqrt(2*pi)));
% Ty = init.Ty/max(Ty);
% init.T0 = 5*Tx*Ty';

u1x = (abs(vect.x-model.p1(1))<=model.W/2);
u1y = (abs(vect.y-model.p1(2))<=model.W/2);
u1xy = u1y*u1x';
u1xy = reshape(u1xy,[par.Nx*par.Ny 1]);
u2x = (abs(vect.x-model.p2(1))<=model.W/2);
u2y = (abs(vect.y-model.p2(2))<=model.W/2);
u2xy = u2y*u2x';
u2xy = reshape(u2xy,[par.Nx*par.Ny 1]);
init.T0 = u1.InitTemp*u1xy+u2.InitTemp*u2xy;
end
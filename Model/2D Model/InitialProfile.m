function init = InitialProfile(model,par,vect,u1,u2)
mu_Tx = model.p1(1);
sigma_Tx = u1.InitWidth;
Tx = (exp((-(vect.x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);
mu_Ty = model.p1(2);
sigma_Ty = u1.InitWidth;
Ty = (exp((-(vect.y-mu_Ty).^2)/(2*sigma_Ty^2))/(sigma_Ty*sqrt(2*pi)));
Ty = Ty/max(Ty);
Txy = reshape(Ty*Tx',[par.Nx*par.Ny 1]);
init.T0 = u1.InitTemp*Txy;

mu_Tx = model.p2(1);
sigma_Tx = u2.InitWidth;
Tx = (exp((-(vect.x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);
mu_Ty = model.p2(2);
sigma_Ty = u2.InitWidth;
Ty = (exp((-(vect.y-mu_Ty).^2)/(2*sigma_Ty^2))/(sigma_Ty*sqrt(2*pi)));
Ty = Ty/max(Ty);
Txy = reshape(Ty*Tx',[par.Nx*par.Ny 1]);

init.T0 = init.T0+u2.InitTemp*Txy;

% u1x = (abs(vect.x-model.p1(1))<=model.W/2);
% u1y = (abs(vect.y-model.p1(2))<=model.W/2);
% u1xy = u1y*u1x';
% u1xy = reshape(u1xy,[par.Nx*par.Ny 1]);
% u2x = (abs(vect.x-model.p2(1))<=model.W/2);
% u2y = (abs(vect.y-model.p2(2))<=model.W/2);
% u2xy = u2y*u2x';
% u2xy = reshape(u2xy,[par.Nx*par.Ny 1]);
% init.T0 = u1.InitTemp*u1xy+u2.InitTemp*u2xy;
end
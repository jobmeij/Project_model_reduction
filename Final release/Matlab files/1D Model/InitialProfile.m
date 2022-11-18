function init = InitialProfile(model,vect,u1,u2)
% u1x = (abs(vect.x-model.p1(1))<=model.W/2);
% u2x = (abs(vect.x-model.p2(1))<=model.W/2);
% init.T0 = u1.InitTemp*u1x+u2.InitTemp*u2x;

mu_Tx = model.p1(1);
sigma_Tx = u1.InitWidth;
Tx = (exp((-(vect.x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);
init.T0 = u1.InitTemp*Tx;

mu_Tx = model.p2(1);
sigma_Tx = u2.InitWidth;
Tx = (exp((-(vect.x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);

init.T0 = init.T0+u2.InitTemp*Tx;
end
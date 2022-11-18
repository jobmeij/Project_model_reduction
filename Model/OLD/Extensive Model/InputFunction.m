function input = InputFunction(model,par,vect,u1,u2)
% mu_fx = model.Lx*0.7;
% sigma_fx = model.Lx*0.05;
% fx = (exp((-(vect.x-mu_fx).^2)/(2*sigma_fx^2))/(sigma_fx*sqrt(2*pi)));
% fx = fx/max(fx);
% mu_fy = model.Ly*0.5;
% sigma_fy = model.Ly*0.05;
% fy = (exp((-(vect.y-mu_fy).^2)/(2*sigma_fy^2))/(sigma_fy*sqrt(2*pi)));
% fy = fy/max(fy);
% input.f = 300*fx*fy';
% 
% mu_g = par.Tend*0.3;
% sigma_g = par.Tend*0.1;
% input.g = (exp((-(vect.t-mu_g).^2)/(2*sigma_g^2))/(sigma_g*sqrt(2*pi)));
% input.g = input.g/max(input.g);
% 
% input.q = zeros(par.Nx,par.Ny,par.M);
% for t = 1:par.M
%     input.q(:,:,t) = input.g(t)*input.f;
% end

u1x = (abs(vect.x-model.p1(1))<=model.W/2);
u1y = (abs(vect.y-model.p1(2))<=model.W/2);
input.u1xy = u1y*u1x';
input.u1xy = reshape(input.u1xy,[par.Nx*par.Ny 1]);
input.u1t = u1.Amplitude*((vect.t>=u1.Starttime) .* (vect.t<=u1.Starttime+u1.Duration));
u1 = input.u1xy*input.u1t';
% for t = 1:par.M
%     input.u1(:,:,t) = input.u1t(t)*input.u1xy;
% end
u2x = (abs(vect.x-model.p2(1))<=model.W/2);
u2y = (abs(vect.y-model.p2(2))<=model.W/2);
input.u2xy = u2y*u2x';
input.u2xy = reshape(input.u2xy,[par.Nx*par.Ny 1]);
input.u2t = u2.Amplitude*((vect.t>=u2.Starttime) .* (vect.t<=u2.Starttime+u2.Duration));
u2 = input.u2xy*input.u2t';
% for t = 1:par.M
%     input.u2(:,:,t) = input.u2t(t)*input.u2xy;
% end
input.u = u1+u2;
end
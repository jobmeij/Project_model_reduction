function input = InputFunction(model,par,vect,u1,u2)
u1x = (abs(vect.x-model.p1(1))<=model.W/2);
u1y = (abs(vect.y-model.p1(2))<=model.W/2);
input.u1xy = u1y*u1x';
input.u1xy = reshape(input.u1xy,[par.Nx*par.Ny 1]);
input.u1t = u1.Amplitude*((vect.t>=u1.Starttime) .* (vect.t<=u1.Starttime+u1.Duration));
u1 = input.u1xy*input.u1t';

u2x = (abs(vect.x-model.p2(1))<=model.W/2);
u2y = (abs(vect.y-model.p2(2))<=model.W/2);
input.u2xy = u2y*u2x';
input.u2xy = reshape(input.u2xy,[par.Nx*par.Ny 1]);
input.u2t = u2.Amplitude*((vect.t>=u2.Starttime) .* (vect.t<=u2.Starttime+u2.Duration));
u2 = input.u2xy*input.u2t';

input.u = u1+u2;
end
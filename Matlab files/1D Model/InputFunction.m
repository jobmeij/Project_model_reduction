function input = InputFunction(model,vect,u1,u2)
input.u1x = (abs(vect.x-model.p1(1))<=model.W/2);
input.u1t = u1.Amplitude*((vect.t>=u1.Starttime) .* (vect.t<=u1.Starttime+u1.Duration));
u1xt = input.u1x*input.u1t';

input.u2x = (abs(vect.x-model.p2(1))<=model.W/2);
input.u2t = u2.Amplitude*((vect.t>=u2.Starttime) .* (vect.t<=u2.Starttime+u2.Duration));
u2xt = input.u2x*input.u2t';

input.u = u1xt+u2xt;
end
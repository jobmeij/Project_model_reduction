function init = InitialProfile(model,vect,u1,u2)
u1x = (abs(vect.x-model.p1(1))<=model.W/2);
u2x = (abs(vect.x-model.p2(1))<=model.W/2);
init.T0 = u1.InitTemp*u1x+u2.InitTemp*u2x;
end
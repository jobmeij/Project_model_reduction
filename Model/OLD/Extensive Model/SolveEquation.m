%% TODO check function

function sys = SolveEquation(model,par,init,input,phi)
sys.A = phi.xy'*phi.x2norm + phi.xy'*phi.y2norm;
%sys.A = phi.xy*phi.x2norm' + phi.xy*phi.y2norm';
sys.A = (model.k/(model.rho*model.c))*sys.A;
%sys.A = (model.k/(model.rho*model.c))*sys.A';

sys.B = zeros(par.r,2);
for i = 1:par.r
    sys.B(i,1) = phi.xy(:,i)'*input.u1xy;
    sys.B(i,2) = phi.xy(:,i)'*input.u2xy;
end
sys.B = (1/(model.rho*model.c))*sys.B;

sys.a = zeros(par.r,par.M);
sys.a_dot = zeros(par.r,par.M);
for i = 1:par.r
    sys.a(i,1) = init.T0'*phi.xynorm(:,i);
end
for t = 1:par.M-1
    sys.a_dot(:,t) = sys.A*sys.a(:,t) + sys.B*[input.u1t(t); input.u2t(t)];
    sys.a(:,t+1) = sys.a(:,t) + sys.a_dot(:,t)*par.dt;
end
end
function sys = SolveEquation(model,par,init,input,phi)
sys.A = phi.x'*phi.xdd*par.dx;
%sys.A = diag(-((pi^2*(1:par.rx).^2)/model.Lx^2));
sys.A = (model.k/(model.rho*model.c))*sys.A;

sys.B = zeros(par.rx,2);
for i = 1:par.rx
    sys.B(i,1) = phi.x(:,i)'*input.u1x;
    sys.B(i,2) = phi.x(:,i)'*input.u2x;
end
sys.B = (1/(model.rho*model.c))*sys.B*par.dx;

sys.a = zeros(par.rx,par.M);
sys.a_dot = zeros(par.rx,par.M);
for i = 1:par.rx
    sys.a(i,1) = init.T0'*phi.x(:,i);
end
sys.a = sys.a*par.dx;

for t = 1:par.M-1
    sys.a_dot(:,t) = sys.A*sys.a(:,t) + sys.B*[input.u1t(t); input.u2t(t)];
    sys.a(:,t+1) = sys.a(:,t) + sys.a_dot(:,t)*par.dt;
end
end
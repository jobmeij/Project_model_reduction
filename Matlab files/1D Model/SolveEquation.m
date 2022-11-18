function sys = SolveEquation(model,par,init,input,phi,e)
if model.Homogeneous
    pc = zeros(par.rx,par.rx);
    for i2 = 1:par.rx
        pc = pc + (e(i2)*((phi.x(:,i2).*phi.x)'*phi.x*par.dx));
    end
    %pc = model.rho(2)*model.c(2);
    invpc = pc^-1;
    kpc = model.k*invpc;
else
    pc = model.rho(1)*model.c(1);
    invpc = pc^-1;
    kpc = model.k*invpc;
end

if model.Analytic
    sys.A = diag([0 -((pi^2*(2:par.rx).^2)/model.Lx^2)]); 
else
    sys.A = phi.x'*phi.xdd*par.dx;
end
sys.A = kpc*sys.A;

sys.B = zeros(par.rx,2);
for i = 1:par.rx
    sys.B(i,1) = phi.x(:,i)'*input.u1x*par.dx;
    sys.B(i,2) = phi.x(:,i)'*input.u2x*par.dx;
end

sys.B(:,1) = (invpc*sys.B(:,1));
sys.B(:,2) = (invpc*sys.B(:,2));

sys.a = zeros(par.rx,par.M);
sys.a_dot = zeros(par.rx,par.M);
for i = 1:par.rx
    sys.a(i,1) = init.T0'*phi.x(:,i)*par.dx;
end
sys.a = sys.a;

for t = 1:par.M-1
    sys.a_dot(:,t) = sys.A*sys.a(:,t) + sys.B*[input.u1t(t); input.u2t(t)];
    sys.a(:,t+1) = sys.a(:,t) + sys.a_dot(:,t)*par.dt;
end
end
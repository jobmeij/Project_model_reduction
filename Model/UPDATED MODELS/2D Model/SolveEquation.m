function sys = SolveEquation(model,par,init,input,phi,e)
if model.Homogeneous
    pc = zeros(par.r,par.r);
    for i2 = 1:par.r
        pc = pc + ((e(i2)*((phi.xy(:,i2).*phi.xy)'*phi.xy*par.dx*par.dy)));
    end
%     pc = model.rho(2)*model.c(2);
    invpc = pc^-1;
    kpc = invpc*model.k;
else
    pc = model.rho(1)*model.c(1);
    invpc = 1/pc;
    kpc = invpc*model.k;
end

if model.Analytic
    % Fill A matrix with pre calculated values
    ki = mod(((2:par.r)-1),par.rx);
    li = floor(((2:par.r)-1)/par.rx);
    sys.A = diag([0, -(((pi*ki)/(model.Lx)).^2 +((pi*li)/(model.Ly)).^2)]);
else
    % Calculate A matrix via inner product with basis function
    sys.A = phi.xy'*phi.x2*par.dx*par.dy + phi.xy'*phi.y2*par.dx*par.dy;
end
sys.A = kpc*sys.A;

sys.B = zeros(par.r,2);
for i = 1:par.r
    sys.B(i,1) = phi.xy(:,i)'*input.u1xy*par.dx*par.dy;
    sys.B(i,2) = phi.xy(:,i)'*input.u2xy*par.dx*par.dy;
end
sys.B = invpc*sys.B;

sys.a = zeros(par.r,par.M);
sys.a_dot = zeros(par.r,par.M);
for i = 1:par.r
    sys.a(i,1) = init.T0'*phi.xy(:,i)*par.dx*par.dy;
end

for t = 1:par.M-1
    sys.a_dot(:,t) = sys.A*sys.a(:,t) + sys.B*[input.u1t(t); input.u2t(t)];
    sys.a(:,t+1) = sys.a(:,t) + sys.a_dot(:,t)*par.dt;
end
end
function sys = SolveEquation(model,par,init,input,phi,e,profile,analytic)
if profile
%     pc = zeros(par.r,1);
%     for n = 1:par.r
%         for k = 1:par.rx
%             for l = 1:par.ry
%                 i = ((k-1)*par.ry)+l;
%                 if k < 2 || l < 2
%                     pc(n) = pc(n) + e(i);    
%                 else
%                     pc(n) = pc(n) + e(i)*(phi.xy(:,n)'*phi.xy(:,i).^2*par.dx*par.dy);
%                 end
%             end
%         end
%     end
    pc = model.rho(2)*model.c(2);
    kpc = model.k./pc;
else
    pc = model.rho(1)*model.c(1);
    kpc = model.k/pc;
end

if ~analytic
    % Calculate A matrix via inner product with basis function 
    sys.A = phi.xy'*phi.x2*par.dx*par.dy + phi.xy'*phi.y2*par.dx*par.dy;
else
    % Fill A matrix with pre calculated values
    ki = mod(((2:par.r)-1),par.rx);
    li = floor(((2:par.r)-1)/par.rx);
    sys.A = diag([0, -(((pi*ki)/(model.Lx)).^2 +((pi*li)/(model.Ly)).^2)]);
end
sys.A = kpc.*sys.A;

sys.B = zeros(par.r,2);
for i = 1:par.r
    sys.B(i,1) = phi.xy(:,i)'*input.u1xy*par.dx*par.dy;
    sys.B(i,2) = phi.xy(:,i)'*input.u2xy*par.dx*par.dy;
end
sys.B = sys.B./pc;

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
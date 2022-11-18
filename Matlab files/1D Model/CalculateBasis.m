function phi = CalculateBasis(model,par,vect)
phi.x = zeros(par.Nx,par.rx);
phi.x(:,1) = 1/sqrt(model.Lx);
phi.xdd = zeros(par.Nx,par.rx);
for k = 2:par.rx
    phi.x(:,k) = sqrt(2/model.Lx)*cos(((k-1)*pi*vect.x)/(model.Lx));
    phi.xdd(:,k) = -((pi*(k-1))/(model.Lx))^2*phi.x(:,k);
end
phi.xdd([1 par.Nx],:) = 0;
end
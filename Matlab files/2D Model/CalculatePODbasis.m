function [phi,par] = CalculatePODbasis(par_, data)
par = par_;
par.r = par.rPOD;
[POD.U, POD.S, POD.Y] = svd(data,0);

% figure(2);
% set(2,'Position',[0 442 650 554]);
% plot(log10(diag(POD.S)),'.r');
% ylabel("Log10 Hankel singular values");
% par.r = input('Enter trunctation value: ');

phi.xy = POD.U(:,1:par.r)/sqrt(par.dx*par.dy);
phi.x2 = zeros(par.Nx,par.Ny,par.r);
phi.y2 = zeros(par.Nx,par.Ny,par.r);

phi.xy = reshape(phi.xy,[par.Nx, par.Ny, par.r]);
for i = 1:par.r
[X, Y] = gradient(phi.xy(:,:,i),par.dx,par.dy);
[phi.x2(:,:,i), ~] = gradient(X,par.dx,par.dy);
[~, phi.y2(:,:,i)] = gradient(Y,par.dx,par.dy);
end
phi.xy = reshape(phi.xy,[par.Nx*par.Ny, par.r]);
phi.x2 = reshape(phi.x2,[par.Nx*par.Ny, par.r]);
phi.y2 = reshape(phi.y2,[par.Nx*par.Ny, par.r]);

end
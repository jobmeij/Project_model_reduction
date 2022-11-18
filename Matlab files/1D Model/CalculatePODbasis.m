function [phi,par] = CalculatePODbasis(par_, data)
par = par_;
par.rx = par.rPOD;
[POD.U, POD.S, POD.Y] = svd(data,0);

% figure(2);
% set(2,'Position',[0 442 650 554]);
% plot(log10(diag(POD.S)),'.r','Markersize',15);
% ylabel("Log10 Hankel singular values");
% par.rx = input('Enter trunctation value: ');
% grid on;

phi.x = POD.U(:,1:par.rx)/sqrt(par.dx);
phi.xdd = zeros(par.Nx, par.rx);
for i = 1:par.rx
phi.xdd(:,i) = gradient(gradient(phi.x(:,i),par.dx),par.dx);
end
end
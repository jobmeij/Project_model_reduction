function e = CalculateMaterialProperties(model,par,vect,phi)
X1 = round((model.p3(1)-model.lx/2)/par.dx);
X2 = round((model.p3(1)+model.lx/2)/par.dx);
Y1 = round((model.p3(2)-model.ly/2)/par.dy);
Y2 = round((model.p3(2)+model.ly/2)/par.dy);
Areax = X1:X2;
Areay = Y1:Y2;
Material = ones(par.Ny, par.Nx)*model.rho(1)*model.c(1);
Material(Areay, Areax) = model.rho(2)*model.c(2);
F = reshape(Material,[par.Nx*par.Ny 1]);

e = zeros(1,par.r);
S = zeros(1,par.Nx*par.Ny);
for i = 1:par.r
    e(i) = phi.xy(:,i)'*F*par.dx*par.dy;
    S = S + phi.xy(:,i)'*e(i);
end

S = reshape(S,[par.Ny par.Nx]);

if model.ShowProfile
    figure(2);
    set(2,'Position',[620 42 650 954]);
    subplot(211);
    surf(vect.x, vect.y, Material,'linestyle','none','FaceColor','interp');
    view(0,90);
    colorbar;
    xlabel("X");
    ylabel("Y");
    grid on;
    subplot(212);
    surf(vect.x, vect.y, S,'linestyle','none','FaceColor','interp');
    caxis([min(F) max(F)]);
    view(0,90);
    colorbar;
    xlabel("X");
    ylabel("Y");
    grid on;
end
end
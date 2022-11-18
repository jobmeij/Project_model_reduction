function PlotTemperature(T,model,par,vect,init,input,res)
T = reshape(T,[par.Nx par.Ny par.M]);
input.u = reshape(input.u,[par.Nx par.Ny par.M]);
init.T0 = reshape(init.T0,[par.Nx par.Ny])+model.Tamb;

Tmin = min(T,[],'all');
Tmax = max([max(T,[],'all') Tmin+1]);
Inputmin = min(input.u,[],'all');
Inputmax = max([max(input.u,[],'all'), Inputmin +1]);
figure(1);
set(1,'Position',[960 42 960 954]);
dt = par.Tend/par.M;
for tl = 1:ceil(res/dt):par.M+ceil(res/dt)
    t = min(tl,par.M);
    figure(1);
    subplot(211);
    surf(vect.x, vect.y, T(:,:,t),'linestyle','none','FaceColor','interp');
    xlim([0 model.Lx]);
    ylim([0 model.Ly]);
    zlim([Tmin Tmax]);
    caxis([Tmin Tmax]);
%     view(0,90);
    xlabel("X");
    ylabel("Y");
    zlabel("Temperature");
    title("Temperature profile @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    grid on;
    colorbar;
    colormap(jet);
    subplot(212);
    surf(vect.x, vect.y, input.u(:,:,t),'linestyle','none');
    xlim([0 model.Lx]);
    ylim([0 model.Ly]);
    zlim([Inputmin Inputmax]);
    caxis([Inputmin Inputmax]);
%     view(0,90);
    xlabel("X");
    ylabel("Y");
    zlabel("Input");
    title("Input power @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    grid on;
    colorbar;
    colormap(jet);
    drawnow;
end

subplot(212)
surf(vect.x, vect.y, T(:,:,1),'linestyle','none','FaceColor','interp');
caxis([Tmin Tmax]);
view(0,90);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
title("Initial temperature profile");
grid on;
colorbar;
colormap(jet);

end
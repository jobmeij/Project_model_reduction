function PlotTemperature(T,model,par,vect,init,input,speed,fig)
T = reshape(T,[par.Ny par.Nx par.M]);
input.u = reshape(input.u,[par.Ny par.Nx par.M]);
init.T0 = reshape(init.T0,[par.Ny par.Nx])+model.Tamb;
Tmin = min(T,[],'all');
Tmax = max([max(T,[],'all') Tmin+1]);
Inputmin = min(input.u,[],'all');
Inputmax = max([max(input.u,[],'all'), Inputmin +1]);

figure(fig);
set(fig,'Position',[1270 42 650 954]);

T_axis = subplot(211);
T_surf = surf('linestyle','none','FaceColor','interp');
xlim([0 model.Lx]);
ylim([0 model.Ly]);
zlim([Tmin Tmax]);
caxis([Tmin Tmax]);
view(0,90);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
grid on;
colorbar;
colormap(jet);

Tinput_axis = subplot(212);
Tinput_surf = surf('linestyle','none','FaceColor','interp');
xlim([0 model.Lx]);
ylim([0 model.Ly]);
zlim([Inputmin Inputmax]);
caxis([Inputmin Inputmax]);
view(0,90);
xlabel("X");
ylabel("Y");
zlabel("Input");
grid on;
colorbar;
colormap(jet);

fps = 30;
dt = par.Tend/par.M;
for tl = 1:ceil((speed/fps)/dt):par.M+ceil((speed/fps)/dt)
    t = min(tl,par.M);
    set(T_surf,'XData', vect.x, 'YData', vect.y, 'ZData', T(:,:,t));
    title(T_axis,"Temperature profile @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    set(Tinput_surf,'XData', vect.x, 'YData', vect.y, 'ZData', input.u(:,:,t));
    title(Tinput_axis,"Input power @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    drawnow();
    pause(1/fps);
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
drawnow();
end
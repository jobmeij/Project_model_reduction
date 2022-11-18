function PlotTemperature(T,model,par,vect,init,input,speed,fig)
init.T0 = init.T0+model.Tamb;
Tmin = min(T,[],'all');
Tmax = max([max(T,[],'all') Tmin+1]);
Inputmin = min(input.u,[],'all');
Inputmax = max([max(input.u,[],'all'), Inputmin +1]);

figure(fig);
set(fig,'Position',[1320 300 600 700]);
T_axis = subplot(211);
T_plot = plot(vect.x,squeeze(T(:,1,:)));
set(T_plot,{'color'},{'r';'b'});
xlim([0 model.Lx]);
ylim([Tmin Tmax]);
xlabel("X");
ylabel("Temperature");
grid on;

Tinput_axis = subplot(212);
Tinput_plot = plot(vect.x,input.u(:,1),'-b');
xlim([0 model.Lx]);
ylim([Inputmin Inputmax]);
xlabel("X");
ylabel("Input");
grid on;

fps = 30;
dt = par.Tend/par.M;
for tl = 1:ceil((speed/fps)/dt):par.M+ceil((speed/fps)/dt)
    t = min(tl,par.M);
    set(T_plot, {'YData'}, num2cell(squeeze(T(:,t,:)),1)');
    title(T_axis,"Temperature profile @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    set(Tinput_plot, 'YData', input.u(:,t));
    title(Tinput_axis,"Input power @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    drawnow;
    pause(1/fps);
end

figure(fig);
subplot(212)
plot(vect.x, squeeze(T(:,1,:)),'r');
xlim([0 model.Lx]);
ylim([Tmin Tmax]);
xlabel("X");
ylabel("Temperature");
title("Initial temperature profile");
grid on;
drawnow();
end
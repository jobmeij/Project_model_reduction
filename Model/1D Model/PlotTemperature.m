function PlotTemperature(T,model,par,vect,init,input,res,fig)
init.T0 = init.T0+model.Tamb;

Tmin = min(T,[],'all');
Tmax = max([max(T,[],'all') Tmin+1]);
Inputmin = min(input.u,[],'all');
Inputmax = max([max(input.u,[],'all'), Inputmin +1]);
figure(fig);
set(fig,'Position',[1320 300 600 700]);
dt = par.Tend/par.M;
for tl = 1:ceil(res/dt):par.M+ceil(res/dt)
    t = min(tl,par.M);
    figure(fig);
    subplot(211);
    plot(vect.x, T(:,t),'r');
    xlim([0 model.Lx]);
    ylim([Tmin Tmax]);
    xlabel("X");
    ylabel("Temperature");
    title("Temperature profile @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    grid on;
    subplot(212);
    plot(vect.x, input.u(:,t),'b');
    xlim([0 model.Lx]);
    ylim([Inputmin Inputmax]);
    xlabel("X");
    ylabel("Input");
    title("Input power @ time: "+num2str(round(vect.t(t),0))+" [sec]");
    grid on;
    drawnow;
end
figure(fig);
subplot(212)
plot(vect.x, T(:,1),'r');
xlim([0 model.Lx]);
ylim([Tmin Tmax]);
xlabel("X");
ylabel("Temperature");
title("Initial temperature profile");
grid on;
end
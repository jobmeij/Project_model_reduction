clear all;
close all hidden;
clc;

k = 0.31;
c = 3890;
p = 1100;

Lx = -10;
Ly = -10;
La = -0.1;

Ax = [0 1;Lx 0];
Ay = [0 1;Ly 0];
Aa = La;

X = 100;
dx = 0.01;

Y = 100;
dy = 0.01;

T = 100;
dt = 0.1;

Fx = zeros(2,X);
Fy = zeros(2,Y);
a = zeros(1,T);

Fx(:,1) = [1;50];
Fy(:,1) = [1;50];
a(1) = 1;

for x = 1:X-1
    Fx(:,x+1) = Fx(:,x) + Ax*Fx(:,x)*dx;
end
for y = 1:Y-1
    Fy(:,y+1) = Fy(:,y) + Ay*Fy(:,y)*dy;
end
for t = 1:T-1
    a(t+1) = a(t) + Aa*a(t)*dt;
end

s = Fx(1,:)'*Fy(1,:);

s = reshape(s,[1,X*Y]);

Temp = a'*s;

Temp = reshape(Temp,[T X Y]);
Temp = permute(Temp,[2 3 1]);
Tmax = max(Temp,[],'all');
Tmin = min(Temp,[],'all');
figure(1);
for t=1:T
    surf(0:dx:X*dx-dx, 0:dy:Y*dy-dy, Temp(:,:,t),'linestyle','none','FaceColor','interp');
    grid on;
    colorbar;
    colormap(jet);
    xlim([0 X*dx-dx]);
    ylim([0 Y*dy-dy]);
    zlim([Tmin Tmax]);
    caxis([Tmin Tmax]);
    %view(0,90);
    xlabel("X");
    ylabel("Y");
    zlabel("Temperature");
    title("Temperature @ "+t*dt+"[sec]");
    drawnow();
end
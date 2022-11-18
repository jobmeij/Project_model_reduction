clear all;
close all hidden;
clc;

% MODEL PARAMETERS
disp("Set model parameters");
Lx_par = 0.2;
Ly_par = 0.15;
rho_par = 1100;
c_par = 3890;
k_par = 0.31;
p1_par = [(1*Lx_par)/4,Ly_par/2];
p2_par = [(3*Lx_par)/4,Ly_par/2];
W_par = 0.05;
Tamb_par = 309;

% REDUCTION PARAMETERS
disp("Set reduction parameters");
Tend = 2000 %600;

rx = 27 %11;
ry = 27 %11;
r = rx*ry;

Nx = 80 %80;
Ny = 80 %80;
M = 600;

dx = Lx_par/(Nx-1);
dy = Ly_par/(Ny-1);
dt = Tend/(M-1);

x_vect = (0:dx:Lx_par)';
y_vect = (0:dy:Ly_par)';
t_vect = (0:dt:Tend)';

% INITIAL PROFILE
disp("Calculate initial profile");
mu_Tx = Lx_par*0.2;
sigma_Tx = Lx_par*0.05;
Tx = (exp((-(x_vect-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);
mu_Ty = Ly_par*0.2;
sigma_Ty = Ly_par*0.1;
Ty = (exp((-(y_vect-mu_Ty).^2)/(2*sigma_Ty^2))/(sigma_Ty*sqrt(2*pi)));
Ty = Ty/max(Ty);
T0 = 3*Tx*Ty';

% INPUT FUNCTION
disp("Calculate initial input");
mu_fx = Lx_par*0.7;
sigma_fx = Lx_par*0.05;
fx = (exp((-(x_vect-mu_fx).^2)/(2*sigma_fx^2))/(sigma_fx*sqrt(2*pi)));
fx = fx/max(fx);
mu_fy = Ly_par*0.5;
sigma_fy = Ly_par*0.05;
fy = (exp((-(y_vect-mu_fy).^2)/(2*sigma_fy^2))/(sigma_fy*sqrt(2*pi)));
fy = fy/max(fy);
f = 300*fx*fy';

mu_g = Tend*0.3;
sigma_g = Tend*0.1;
g = (exp((-(t_vect-mu_g).^2)/(2*sigma_g^2))/(sigma_g*sqrt(2*pi)));
g = g/max(g);
q = zeros(Nx,Ny,M);
for t = 1:M
    q(:,:,t) = g(t)*f;
end

% BASIS FUNCTION'S AND IT'S GRADIENTS
disp("Calculate basis functions and gradients");
phix = zeros(Nx,rx);
phix2 = zeros(Nx,rx);
phix(:,1) = 1/sqrt(Lx_par);
for k = 2:rx
    phix(:,k) = sqrt(2/Lx_par)*cos(((k-1)*pi*x_vect)/(Lx_par));
    phix2(:,k) = -((pi^2*k^2)/(Lx_par^2))*phix(:,k);
end

phiy = zeros(Ny,ry);
phiy2 = zeros(Ny,ry);
phiy(:,1) = 1/sqrt(Ly_par);
for l = 2:ry
    phiy(:,l) = sqrt(2/Ly_par)*cos(((l-1)*pi*y_vect)/(Ly_par));
    phiy2(:,l) = -((pi^2*l^2)/(Ly_par^2))*phiy(:,l);
end

phi = zeros(Nx,Ny,r);
phinorm = zeros(Nx,Ny,r);
phi2x = zeros(Nx,Ny,r);
phi2y = zeros(Nx,Ny,r);
for k = 1:rx
    for l = 1:ry
        i = ((k-1)*ry)+l;
        phi(:,:,i) = phix(:,k)*phiy(:,l)';
    end
end

for k = 1:r
    [phi1x, phi1y] = gradient(phi(:,:,k));
    [phi2x(:,:,k),~] = gradient(phi1x);
    [~,phi2y(:,:,k)] = gradient(phi1y);
    phinorm(:,:,k) = phi(:,:,k)/norm(phi(:,:,k),2);
end

% LINEAR EQUATIONS CONSTRUCTION AND SOLVER
disp("Construct ans solve linear equations");
A = zeros(r,r);
for i = 1:r
    for j = 1:r
        for x = 1:Nx
            for y = 1:Ny
                A(i,j) = A(i,j) + phi(x,y,i)*phi2x(x,y,j) + phi(x,y,i)*phi2y(x,y,j);
            end
        end
    end
end
A = (k_par/(rho_par*c_par))*A;

B = zeros(r,1);
for i = 1:r
    for x = 1:Nx
        for y = 1:Ny
            B(i) = B(i) + phi(x,y,i)*f(x,y);
        end
    end
end
B = (1/(rho_par*c_par))*B;

a = zeros(r,M);
for k = 1:r
    for x = 1:Nx
        for y = 1:Ny
            a(k,1) = a(k,1) + T0(x,y)*phinorm(x,y,k);
        end
    end
end
for t = 1:M-1
    a_dot = A*a(:,t) + B*g(t);
    a(:,t+1) = a(:,t) + a_dot*dt;
end

% OUTPUT TEMPERATURE PROFILE
disp("Calculate output");
T = zeros(Nx,Ny,M);
for k = 1:r
    for x = 1:Nx
        for y = 1:Ny
            for t = 1:M
                T(x,y,t) = T(x,y,t) + a(k,t)*phinorm(x,y,k);
            end
        end
    end
end
T = T + Tamb_par;

% TEMPERATURE PLOT
disp("Plot output");
Tmin = min(T,[],'all');
Tmax = max([max(T,[],'all') Tmin+1]);
Inputmin = min(q(:,:,:),[],'all');
Inputmax = max([max(q(:,:,:),[],'all'), Inputmin +1]);

figure(1);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0, 0.5, 1]);

for t = 1:5:M
    figure(1);
    subplot(211);
    surf(x_vect, y_vect, T(:,:,t),'linestyle','none');
    xlim([0 Lx_par]);
    ylim([0 Ly_par]);
    zlim([Tmin Tmax]);
    caxis([Tmin Tmax]);
    view(0,90);
    xlabel("X");
    ylabel("Y");
    zlabel("Temperature");
    title("Temperature profile @ time: "+num2str(round(t_vect(t),0))+" [sec]");
    grid on;
    colorbar;
    colormap(jet);
    subplot(212);
    surf(x_vect, y_vect, q(:,:,t),'linestyle','none');
    xlim([0 Lx_par]);
    ylim([0 Ly_par]);
    zlim([Inputmin Inputmax]);
    caxis([Inputmin Inputmax]);
    view(0,90);
    xlabel("X");
    ylabel("Y");
    zlabel("Input");
    title("Input power @ time: "+num2str(round(t_vect(t),0))+" [sec]");
    grid on;
    colorbar;
    colormap(jet);
    drawnow;
end

figure(2)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.5, 0.5, 0.5]);
surf(x_vect, y_vect, T0+Tamb_par,'linestyle','none');
caxis([Tmin Tmax]);
view(0,90);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
title("Initial temperature profile");
grid on;
colorbar;
colormap(jet);


% BASIS FUNCTION'S PLOT
if r <= 25 && false
    for i = 1:r
        figure(3);
        subplot(rx,ry,i);
        surf(x_vect, y_vect, phi(:,:,i));
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
        figure(4);
        subplot(rx,ry,i);
        surf(x_vect, y_vect, phi2x(:,:,i));
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
        figure(5);
        subplot(rx,ry,i);
        surf(x_vect, y_vect, phi2y(:,:,i));
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
    end
    
    figure(3);
    sgtitle("Phi");
    figure(4);
    sgtitle("Phi 2X");
    figure(5);
    sgtitle("Phi 2Y");
    
    for i = 1:rx
        figure(6);
        subplot(rx,1,i);
        plot(x_vect,phix(:,i));
        xlabel("X");
        ylabel("Amplitude");
        grid on;
        figure(7);
        subplot(rx,1,i);
        plot(x_vect,phix2(:,i));
        xlabel("X");
        ylabel("Amplitude");
        grid on;
    end
    figure(6);
    sgtitle("Phi X");
    figure(7);
    sgtitle("Phi X2");
    
    for i = 1:ry
        figure(8);
        subplot(ry,1,i);
        plot(y_vect,phiy(:,i));
        xlabel("Y");
        ylabel("Amplitude");
        grid on;
        figure(9);
        subplot(ry,1,i);
        plot(y_vect,phiy2(:,i));
        xlabel("Y");
        ylabel("Amplitude");
        grid on;
    end
    figure(8);
    sgtitle("Phi Y");
    figure(9);
    sgtitle("Phi Y2");
    
end


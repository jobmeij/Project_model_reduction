clear all;
close all hidden;
clc;

% d4w/dx4 + mu d2w/dt2 - q = 0

% Declaring parameters
mu = 210e6;
L = 1;
T = 5;

r = 20;
N = 100;
M = 500;

dx = L/(N-1);
dt = T/(M-1);

% Initial vectors
x_vect = [0:dx:L]';
t_vect = 0:dt:T;

sigma_f = 0.05;
mu_f = 0.2;
f = (exp((-(x_vect-mu_f).^2)/(2*sigma_f^2))/(sigma_f*sqrt(2*pi)));
sigma_g = 1;
mu_g = 0.2;
g = (exp((-(t_vect-mu_g).^2)/(2*sigma_g^2))/(sigma_g*sqrt(2*pi)));
q = f*g;

w = zeros(N,M+1);

% Creating Fourier basis
phi = zeros(N,r);
for k = 1:r
    phi(:,k) = exp(i*(k-1)*pi*x_vect);
end

% Creating phi(4) (4th spacial derivative)
phi4 = zeros(N,r);
for j = 1:r
    phi4(:,j) = gradient(gradient(gradient(gradient(phi(:,j)))));
end

% Creating A(i,j) matrix <phi_i,phi(4)_j>
A = zeros(r,r);
for i = 1:r
    for j = 1:r
        A(i,j) = -(1/mu)*dot(phi(:,i),phi4(:,j));
    end
end
Ae = [zeros(r,r) eye(r,r);A zeros(r,r)];

% Creating B(i) <phi_r,f(x)>*g
B = zeros(r,1);
for i = 1:r
    B(i) = (1/mu)*dot(phi(:,i),f);
end
Be = [zeros(r,1);B];

% Creating ak
e = zeros(r*2,M);
for t = 1:M-1
    e_dot = Ae*e(:,t) + Be*g(t);
    e(:,t+1) = e(:,t) + e_dot*dt;
end
a = e(1:r,:);

% Computing wr - Sum:
for k = 1:r
    % Computing space
    for x = 1:N
        % computing time
        for t = 1:M
            w(x,t) = w(x,t) + a(k,t)*phi(x,k);
        end
    end
end

% Plotting
for t = 1:M
    figure(1);
    subplot(211);
    plot(x_vect,q(:,t));
    grid on;
    xlim([0 L]);
    ylim([-10 10]);
    xlabel("X");
    ylabel("Load");
    title("Input");
    
    subplot(212);
    plot(x_vect,real(w(:,t)));
    grid on;
    xlim([0 L]);
    ylim([-1e-5 1e-5]);
    xlabel("X");
    ylabel("Deflection");
    title("Output");
    drawnow;
end

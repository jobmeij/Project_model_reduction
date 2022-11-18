clear all;
close all hidden;
clc;

% MODEL PARAMETERS
model.Lx = 0.2;
model.Ly = 0.15;

% SIMULATION PARAMETERS
par.rx = 20;
par.ry = 20;
par.r = par.rx*par.ry;
par.Nx = 100;
par.Ny = 100;
par.dx = model.Lx/(par.Nx-1);
par.dy = model.Ly/(par.Ny-1);

vect.x = (0:par.dx:model.Lx)';
vect.y = (0:par.dy:model.Ly)';

%x harmonics
phix = zeros(par.Nx,par.rx);
phix(:,1) = 1/sqrt(model.Lx);
for k = 2:par.rx
    phix(:,k) = sqrt(2/model.Lx)*cos(((k-1)*pi*vect.x)/(model.Lx));
end

%y harmonics
phiy = zeros(par.Ny,par.ry);
phiy(:,1) = 1/sqrt(model.Ly);
for l = 2:par.ry
    phiy(:,l) = sqrt(2/model.Ly)*cos(((l-1)*pi*vect.y)/(model.Ly));
end

%xy harmonics
phi.xy = zeros(par.Ny,par.Nx,par.r);
for l = 1:par.ry
    for k = 1:par.rx
        i = ((l-1)*par.rx)+k;
        phi.xy(:,:,i) = phiy(:,l)*phix(:,k)';
    end
end

phi.xy = reshape(phi.xy,[par.Nx*par.Ny,par.r]);

Orthonormal = phi.xy'*phi.xy*par.dx*par.dy;
Orthonormal(Orthonormal <= 0.1) = 0;
spy(Orthonormal);

phi.xy = reshape(phi.xy,[par.Nx,par.Ny,par.r]);
figure(2);
set(2,'Position',[100 200 800 800]);
subplot(221);
surfploti = surf('linestyle','none','FaceColor','interp');
set(surfploti,'XData', vect.x, 'YData', vect.y ,'ZData', phi.xy(:,:,1));
xlabel('X');
ylabel('Y');
zlabel('PHI i');
title('PHI i');
view(0,90);

subplot(222);
surfplotj = surf('linestyle','none','FaceColor','interp');
set(surfplotj,'XData', vect.x, 'YData', vect.y ,'ZData', phi.xy(:,:,1));
xlabel('X');
ylabel('Y');
zlabel('PHI j');
title('PHI j');
view(0,90);

combinedaxis = subplot(223);
surfplotij = surf('linestyle','none','FaceColor','interp');
set(surfplotij,'XData', vect.x, 'YData', vect.y ,'ZData', phi.xy(:,:,1));
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Combined');
view(0,90);

while 1
    i = input('Basisfunction index i = ');
    j = input('Basisfunction index j = ');
    disp('-----------------------------');
    set(surfploti,'ZData', phi.xy(:,:,i));
    set(surfplotj,'ZData', phi.xy(:,:,j));
    set(surfplotij,'ZData', phi.xy(:,:,i).*phi.xy(:,:,j)*par.dx*par.dy);
    combinedsum = string(sum(phi.xy(:,:,i).*phi.xy(:,:,j)*par.dx*par.dy,'all'));
    title(combinedaxis,'Combined = '+combinedsum);
    drawnow();
end


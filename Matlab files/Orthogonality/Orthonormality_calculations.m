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

%% Determining if solution is orthogonal for given i and j
clear all;
close all hidden;
clc;

% Insert (natural) values for K and L
K = 3;
Lxc = 0.2;
Lyc = 0.15;
% Creating symbols for usage in formulas
syms Lx Ly x y ;

% Insert values for i and j
A = zeros(9,9);
for i = 0:size(A,1)-1
    for j = 0:size(A,2)-1
        ki = mod(i,K);
        li = floor(i/K);
        kj = mod(j,K);
        lj = floor(j/K);
        
        % Distinguish different cases for ki,kj,li and lj TODO check functions!
        if ki == 0 && li == 0
            Fi = 1/sqrt(Lx*Ly);
        elseif ki == 0 && li > 0
            Fi = (1/sqrt(Lx))*sqrt(2/Ly)*cos((li*pi*y)/Ly);
        elseif ki > 0 && li == 0
            Fi = (1/sqrt(Ly))*sqrt(2/Lx)*cos((ki*pi*x)/Lx);
        elseif ki > 0 && li > 0
            Fi = (2/sqrt(Lx*Ly))*cos((ki*pi*x)/Lx)*cos((li*pi*y)/Ly);
        end
        if kj == 0 && lj == 0
            Fj = 1/sqrt(Lx*Ly);
        elseif kj == 0 && lj > 0
            Fj = (1/sqrt(Lx))*sqrt(2/Ly)*cos((lj*pi*y)/Ly);
        elseif kj > 0 && lj == 0
            Fj = (1/sqrt(Ly))*sqrt(2/Lx)*cos((kj*pi*x)/Lx);
        elseif kj > 0 && lj > 0
            Fj = (2/sqrt(Lx*Ly))*cos((kj*pi*x)/Lx)*cos((lj*pi*y)/Ly);
        end
        F = Fi*Fj;
        
        F = subs(F,[Lx Ly],[Lxc Lyc]);
        
        % Integrate function
        %S = int(int(F,y,0,Ly),x,0,Lx);
        S = int(int(F,y,0,Lyc),x,0,Lxc);
        disp('i='+string(i)+' j='+string(j));
        %disp('Answer = '+string(S));
        %S = subs(S,[Lx Ly],[Lxc Lyc]);
        A(i+1,j+1) = double(S);
    end
end
A
%% Determining if solution is orthogonal for given i and j
%%%%%%%%%%% VARIATION EXERCISE 10%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all hidden;
clc;

% Insert (natural) values for K and L
K = 3;
Lxc = 0.2;
Lyc = 0.15;
% Creating symbols for usage in formulas
syms Lx Ly x y ;

% Insert values for i and j
A = zeros(9,9,9);
for i = 0:size(A,1)-1
    for j = 0:size(A,2)-1
        for z = 0:size(A,3)-1
            ki = mod(i,K);
            li = floor(i/K);
            kj = mod(j,K);
            lj = floor(j/K);
            kz = mod(z,K);
            lz = floor(z/K);
            
            % Distinguish different cases for ki,kj,li and lj TODO check functions!
            if ki == 0 && li == 0
                Fi = 1/sqrt(Lx*Ly);
            elseif ki == 0 && li > 0
                Fi = (1/sqrt(Lx))*sqrt(2/Ly)*cos((li*pi*y)/Ly);
            elseif ki > 0 && li == 0
                Fi = (1/sqrt(Ly))*sqrt(2/Lx)*cos((ki*pi*x)/Lx);
            elseif ki > 0 && li > 0
                Fi = (2/sqrt(Lx*Ly))*cos((ki*pi*x)/Lx)*cos((li*pi*y)/Ly);
            end
            if kj == 0 && lj == 0
                Fj = 1/sqrt(Lx*Ly);
            elseif kj == 0 && lj > 0
                Fj = (1/sqrt(Lx))*sqrt(2/Ly)*cos((lj*pi*y)/Ly);
            elseif kj > 0 && lj == 0
                Fj = (1/sqrt(Ly))*sqrt(2/Lx)*cos((kj*pi*x)/Lx);
            elseif kj > 0 && lj > 0
                Fj = (2/sqrt(Lx*Ly))*cos((kj*pi*x)/Lx)*cos((lj*pi*y)/Ly);
            end
            if kz == 0 && lz == 0
                Fz = 1/sqrt(Lx*Ly);
            elseif kz == 0 && lz > 0
                Fz = (1/sqrt(Lx))*sqrt(2/Ly)*cos((lz*pi*y)/Ly);
            elseif kz > 0 && lz == 0
                Fz = (1/sqrt(Ly))*sqrt(2/Lx)*cos((kz*pi*x)/Lx);
            elseif kz > 0 && lz > 0
                Fz = (2/sqrt(Lx*Ly))*cos((kz*pi*x)/Lx)*cos((lz*pi*y)/Ly);
            end
            F = Fi*Fz*Fj;
            %F = Fi*Fj;
            
            F = subs(F,[Lx Ly],[Lxc Lyc]);
            
            % Integrate function
            %S = int(int(F,y,0,Ly),x,0,Lx);
            S = int(int(F,y,0,Lyc),x,0,Lxc);
            
            disp('i='+string(i)+' j='+string(j)+' z='+string(z));
            %disp('Answer = '+string(S));
            
            %S = subs(S,[Lx Ly],[Lxc Lyc]);
            
            A(i+1,j+1,z+1) = double(S);
        end
    end
end
A
%%
clear all;
close all hidden;
clc;

Lx = 0.2;
Ly = 0.15;

Nx = 100;
Ny = 100;
dx = Lx/(Nx-1);
dy = Ly/(Ny-1);
x = (0:dx:Lx)';
y = (0:dy:Ly)';

%%%%%%%%%%%%%%%%INITIAL PROFILE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mu_Tx = 0.1;
sigma_Tx = 0.03;
Tx = (exp((-(x-mu_Tx).^2)/(2*sigma_Tx^2))/(sigma_Tx*sqrt(2*pi)));
Tx = Tx/max(Tx);
mu_Ty = 0.075;
sigma_Ty = 0.03;
Ty = (exp((-(y-mu_Ty).^2)/(2*sigma_Ty^2))/(sigma_Ty*sqrt(2*pi)));
Ty = Ty/max(Ty);
Txy = reshape(Ty*Tx',[Nx*Ny 1]);
T0 = 5*Txy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phi0xy = ones(Nx,Ny)*1/sqrt(Lx*Ly);
phi0xy = reshape(phi0xy,[Nx*Ny,1]);

a0 = T0'*phi0xy*dx*dy;

T0 = reshape(T0,[Nx,Ny]);

Tend = ones(Nx,Ny)*a0*(1/sqrt(Lx*Ly));

figure(1);
surf(x,y,T0,'linestyle','none','FaceColor','interp');
hold on;
surf(x,y,Tend,'linestyle','none','FaceColor','interp');
xlim([0 Lx]);
ylim([0 Ly]);
xlabel("X");
ylabel("Y");
zlabel("Temperature");
grid on;
colorbar;
colormap(jet);


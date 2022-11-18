function phi = CalculateBasis(model,par,vect)

%x harmonics
phix = zeros(par.Nx,par.rx);
phix(:,1) = 1/sqrt(model.Lx);
phixdd = zeros(par.Nx,par.rx);
for k = 2:par.rx
    phix(:,k) = sqrt(2/model.Lx)*cos(((k-1)*pi*vect.x)/(model.Lx));
    phixdd(:,k) = -((pi*(k-1))/(model.Lx))^2*phix(:,k);
end

%y harmonics
phiy = zeros(par.Ny,par.ry);
phiy(:,1) = 1/sqrt(model.Ly);
phiydd = zeros(par.Ny,par.ry);
for l = 2:par.ry
    phiy(:,l) = sqrt(2/model.Ly)*cos(((l-1)*pi*vect.y)/(model.Ly));
    phiydd(:,l) = -((pi*(l-1))/(model.Ly))^2*phiy(:,l);
end

%xy harmonics
phi.xy = zeros(par.Ny,par.Nx,par.r);
phi.x2 = zeros(par.Ny,par.Nx,par.r);
phi.x2 = zeros(par.Ny,par.Nx,par.r);

for l = 1:par.ry
    for k = 1:par.rx
        i = ((l-1)*par.rx)+k;
        phi.xy(:,:,i) = phiy(:,l)*phix(:,k)';
        [X, Y] = gradient(phi.xy(:,:,i),par.dx,par.dy);
        [phi.x2(:,:,i), ~] = gradient(X,par.dx,par.dy);
        [~, phi.y2(:,:,i)] = gradient(Y,par.dx,par.dy);
%         if l == 1
%             phi.x2(:,:,i) = ones(par.Ny,1)*phixdd(:,k)';
%         else
%             phi.x2(:,:,i) = -phiydd(:,l)*phixdd(:,k)';
%         end
%         if k == 1
%             phi.y2(:,:,i) = phiydd(:,l)*ones(1,par.Nx);
%         else
%             phi.y2(:,:,i) = -phiydd(:,l)*phixdd(:,k)';
%         end
    end
end

% phi.x2([1 par.Nx],:,:) = 0;
% phi.y2(:,[1 par.Ny],:) = 0;

phi.xy = reshape(phi.xy,[par.Nx*par.Ny,par.r]);
phi.x2 = reshape(phi.x2,[par.Nx*par.Ny,par.r]);
phi.y2 = reshape(phi.y2,[par.Nx*par.Ny,par.r]);

% phi.xynorm = zeros(par.Nx*par.Ny,par.r);
% phi.x2norm = zeros(par.Nx*par.Ny,par.r);
% phi.y2norm = zeros(par.Nx*par.Ny,par.r);
% 
% phi.xynorm = normalize(phi.xy,'norm');
% phi.xynorm(isnan(phi.xynorm)) = 0;
% phi.x2norm = normalize(phi.x2,'norm');
% phi.x2norm(isnan(phi.x2norm)) = 0;
% phi.y2norm = normalize(phi.y2,'norm');
% phi.y2norm(isnan(phi.y2norm)) = 0;

end
function phi = CalculateBasis(model,par,vect)

%x harmonics
phix = zeros(par.Nx,par.rx);
phix(:,1) = 1/sqrt(model.Lx);
phixdd = zeros(par.Nx,par.rx);
for k = 2:par.rx
    phix(:,k) = sqrt(2/model.Lx)*cos(((k-1)*pi*vect.x)/(model.Lx));
    phixdd(:,k) = -((pi*k)/(model.Lx))^2*phix(:,k);
end

%y harmonics
phiy = zeros(par.Ny,par.ry);
phiy(:,1) = 1/sqrt(model.Ly);
phiydd = zeros(par.Ny,par.ry);
for l = 2:par.ry
    phiy(:,l) = sqrt(2/model.Ly)*cos(((l-1)*pi*vect.y)/(model.Ly));
    phiydd(:,l) = -((pi*l)/(model.Ly))^2*phiy(:,l);
end

%xy harmonics
phi.xy = zeros(par.Nx,par.Ny,par.r);
phi.xynorm = zeros(par.Nx,par.Ny,par.r);
phi.x2norm = zeros(par.Nx,par.Ny,par.r);
phi.y2norm = zeros(par.Nx,par.Ny,par.r);

for k = 1:par.rx
    for l = 1:par.ry
        i = ((k-1)*par.ry)+l;
        j = ((l-1)*par.rx)+k;
        phi.xy(:,:,i) = phix(:,k)*phiy(:,l)';
        normxy = norm(phi.xy(:,:,i));
        if normxy ~= 0
            phi.xynorm(:,:,i) = phi.xy(:,:,i)/normxy;
        end
        if l == 1
            phix2 = ones(par.Ny,1)*phixdd(:,k)';
        else
            phix2 = -phiydd(:,l)*phixdd(:,k)';
        end
        normx2 = norm(phix2);
        if normx2 ~= 0
            phi.x2norm(:,:,j) = phix2/normx2;
        end
        if k == 1
            phiy2 = phiydd(:,l)*ones(1,par.Nx);
        else
            phiy2 = -phiydd(:,l)*phixdd(:,k)';
        end
        normy2 = norm(phiy2);
        if normy2 ~= 0
            phi.y2norm(:,:,j) = phiy2/normy2;
        end
    end
end

phi.x2norm([1 par.Nx],:,:) = 0;
phi.y2norm(:,[1 par.Ny],:) = 0;

phi.xy = reshape(phi.xy,[par.Nx*par.Ny,par.r]);
phi.xynorm = reshape(phi.xynorm,[par.Nx*par.Ny,par.r]);
phi.x2norm = reshape(phi.x2norm,[par.Nx*par.Ny,par.r]);
phi.y2norm = reshape(phi.y2norm,[par.Nx*par.Ny,par.r]);
end


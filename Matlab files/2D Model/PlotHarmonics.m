function PlotHarmonics(phi,par,vect,bool)
if par.r <= 25 && bool

    phi.xy = reshape(phi.xy,[par.Nx par.Ny par.r]);
    phi.x2norm = reshape(phi.x2norm,[par.Nx par.Ny par.r]);
    phi.y2norm = reshape(phi.y2norm,[par.Nx par.Ny par.r]);

    w = 600;
    h = 600;
    b = 42;
    sw = 1920;
    sh = 1080;

    figure(3);
    sgtitle("Phi");
    set(3,'Position',[sw-3*w sh-(h+b)*1 w h-b]);
    figure(4);
    sgtitle("Phi X2");
    set(4,'Position',[sw-2*w sh-(h+b)*1 w h-b]);
    figure(5);
    sgtitle("Phi Y2");
    set(5,'Position',[sw-1*w sh-(h+b)*1 w h-b]);

    for i = 1:par.r
        figure(3);
        subplot(par.rx,par.ry,i);
        surf(vect.x, vect.y, phi.xy(:,:,i),'linestyle','none');
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
        figure(4);
        subplot(par.rx,par.ry,i);
        surf(vect.x, vect.y, phi.x2norm(:,:,i),'linestyle','none');
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
        figure(5);
        subplot(par.rx,par.ry,i);
        surf(vect.x, vect.y, phi.y2norm(:,:,i),'linestyle','none');
        view(0,90);
        xlabel("X");
        ylabel("Y");
        zlabel("Amplitude");
        grid on;
    end
end
end
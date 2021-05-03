function [Z,Phi,fres,Zmax,QZ] = Impedance2D(a,b,c,d,ff)

    % Impedance computation for a 2D linear system (analytical)
    %   resonant frequency:     fres
    %   maximal impedance:      zmax
    %   QZ-value:               QZ 
    
    
    deltaf = ff(2)-ff(1);
    omega = 2*pi*ff/1000;
    Z = sqrt((d^2+omega.^2)./((a*d-b*c-omega.^2).^2+(a+d)^2*omega.^2));
    A2 = (omega.*(a*d-b*c-omega.^2)-(a+d)*omega*d);
    A1 = (-d*(a*d-b*c-omega.^2)-(a+d)*omega.^2);
    phiaux = A2./A1;
    Phi = -atan(phiaux);
    for j=1:length(ff)
        if A1(j) < 0 && A2(j) > 0 
             Phi(j) = Phi(j) - pi;
        elseif A1(j) < 0 && A2(j) < 0 
             Phi(j) = Phi(j)+pi;
        end
    end
    
    [Zmax,jfres] = max(Z);
    fres = (jfres-1)*deltaf;
    QZ = Zmax-Z(1);
    
    
    GRPH = 0;
    if GRPH == 1
        figure
        hold on
        plot(ff,Z,'-b','linewidth',2);   
        plot([0 max(ff)],[0 0],':');
        axis([0 ff(end) 0 max(Z)*1.2])    
        set(gca,'fontsize',24);
        xlabel('f  [Hz]');
        ylabel('Z');
        title('Impedance profile')

        figure
        hold on
        plot(ff,Phi,'-b','linewidth',2);    
        plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
        plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
        plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
        plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
        axis([0 ff(end) -pi/2-0.2 pi/2+0.2])
        set(gca,'fontsize',24);
        xlabel('f  [Hz]');
        ylabel('\Phi');
        title('Phase profile');
    end
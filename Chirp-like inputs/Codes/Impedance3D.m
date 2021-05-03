function [Z,Phi,fres,fares,fphas,faphas,Zmax,QZ,Qo] = Impedance3D(a,b,c,d,alp,bet,gma,ff)

    % Impedance computation for a 3D linear system (analytical)
    %   resonant frequency:     fres
    %   maximal impedance:      zmax
    %   QZ-value:               QZ 
    
  
     
    omega = 2*pi*ff/1000;
    
    Preal = bet*d-omega.^2;
    Pimag = -omega*(bet+d);
    Qreal = omega.^2*(a+d+bet)-bet*(a*d-b*c)+d*alp*gma;
    Qimag = omega.*(a*d-b*c-alp*gma+a*bet+d*bet-omega.^2);
    Z = sqrt((Preal.^2+Pimag.^2)./(Qreal.^2+Qimag.^2));
    A1 = (Preal.*Qreal+Pimag.*Qimag);
    A2 = (-Preal.*Qimag+Pimag.*Qreal);
    A1 = A1./(Qreal.^2+Qimag.^2);
    A2 = A2./(Qreal.^2+Qimag.^2);
    Phiaux = A2./A1;
    Phi = -atan(Phiaux);
    for j=1:length(ff)
        if A1(j) < 0 && A2(j) >= 0 
             Phi(j) = Phi(j) + pi;
        elseif A1(j) < 0 && A2(j) < 0 
             Phi(j) = Phi(j)-pi;
        end
    end
    
    fres=0;
    fares=0;
    Zmax=Z(1);
    Zmin=Z(1);
    if Z(2)>Z(1)
        for j=2:length(ff)
            if Z(j)<=Z(j-1)
                jres=j;
                fres=ff(jres);
                Zmax=Z(jres);
                break;
            end
        end
        Zmin=Z(1);
    elseif Z(2)<Z(1)
        for j=2:length(ff)
            if Z(j)>Z(j-1)
                jares=j;
                fares=ff(jares);
                Zmin=Z(jares);
                break;
            end
        end
        if j<length(ff)
            for j=jares+1:length(ff)
                if Z(j)<Z(j-1)
                    jres=j;
                    fres=ff(jres);
                    Zmax=Z(jres);
                    break;
                end
            end
        end
    end

    QZ=Zmax-Zmin;
    Qo=Zmax-Z(1);

    fphas=0;
    faphas=0;
    if Phi(2)<0 
        for j=2:length(ff)
            if Phi(j)>=0
                jphas=j;
                fphas=ff(j);
                break;
            end        
        end
    elseif Phi(2)>0
        for j=2:length(ff)
            jcount=j;
            if Phi(j)<=0
                japhas=j;
                faphas=ff(japhas);
                break;
            end
        end
        if j<ff(end)
            for j=japhas+1:length(ff)
                if Phi(j)>=0
                    jphas=j;
                    fphas=ff(jphas);
                    break;
                end
            end
        end
    end
    
    fgraphmax = 100;
    
    figure
    plot(ff,Z,'-b','linewidth',2);
    hold on
    plot([0 max(ff)],[0 0],':');
    plot([0 0],[-10 100],'--b','linewidth',1);
    axis([-10 fgraphmax 0 max(Z)*1.2])
    set(gca,'fontsize',20);
    xlabel('f [Hz]');
    ylabel('Z');
    title('Impedance profile')

    figure
    plot(ff,Phi,'-b','linewidth',2);
    hold on
    plot([0 max(ff)],[0 0],'--','Color',[.7 .7 .7]);
    plot([0 max(ff)],[pi/2 pi/2],'--','Color',[.7 .7 .7]);
    plot([0 max(ff)],[-pi/2 -pi/2],'--','Color',[.7 .7 .7]);
    plot([0 max(ff)],[2*pi 2*pi],'--','Color',[.7 .7 .7]);
    plot([0 0],[-10 10],'--b','linewidth',1);
    axis([-10 fgraphmax -pi pi])
    set(gca,'fontsize',20);
    %axis([0 1000 -pi/2 pi/2])
    xlabel('f [Hz]');
    ylabel('\Phi');
    title('Phase profile')
    
    
    
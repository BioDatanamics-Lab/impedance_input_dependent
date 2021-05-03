function [r,mu,fnat] = Eigenvalues3D(a,b,c,d,alp,bet,gma)

    % Eigenvalues computation for a 3D linear system
    
    polynch1 = 1;
    polynch2 = -(a+d+bet);
    polynch3 = (a*d+a*bet+d*bet-b*c-alp*gma);
    polynch4 = (b*c-a*d)*bet+gma*alp*d;
    r = roots([polynch1 polynch2 polynch3 polynch4]);

    mu = abs(imag(r))*1000/(2*pi);
    for j=1:3
        if mu(j) > 0
            fnat = mu(j);
        end
    end
    if mu(1)==0 && mu(2) == 0
        fnat = 0;
    end

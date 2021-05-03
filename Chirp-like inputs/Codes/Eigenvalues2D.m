function [r,mu,fnat] = Eigenvalues2D(a,b,c,d)

    % Eigenvalues computation for a 2D linear system
    
    r(1) = (a+d + sqrt((a-d)^2+4*b*c))/2;
    r(2) = (a+d - sqrt((a-d)^2+4*b*c))/2;
    if (a-d)^2+4*b*c<0
        mu = sqrt(-(4*b*c+(a-d)^2))/2;
    else
        mu = 0;
    end
    fnat = mu*1000/(2*pi);
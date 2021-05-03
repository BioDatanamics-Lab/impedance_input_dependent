function [peakt,Vmax,cnt] = PeaksOsc(V,t,tmin,tmax,dt)

jmin = floor(tmin/dt);
jmax = floor(tmax/dt);
peakt = zeros(1,1);
Vmax = zeros(1,1);
cnt = 0;
for j=jmin+2:jmax-1          
    if V(j)>V(j-1) && V(j)>=V(j+1)
        cnt=cnt+1;
        peakt(cnt)=t(j); 
        Vmax(cnt) = V(j);
    end
end
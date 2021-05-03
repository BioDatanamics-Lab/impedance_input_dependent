function [trought,Vmin,cnt] = TroughsOsc(V,t,tmin,tmax,dt)

jmin = floor(tmin/dt);
jmax = floor(tmax/dt);
trought = zeros(1,1);
Vmin = zeros(1,1);
cnt = 0;
for j=jmin+2:jmax-1          
    if V(j)<V(j-1) && V(j)<=V(j+1)
        cnt=cnt+1;
        trought(cnt)=t(j); 
        Vmin(cnt) = V(j);
    end
end
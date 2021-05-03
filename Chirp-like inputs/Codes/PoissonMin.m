function [tspk,ISI,cnt] = PoissonMin(Freq,Tmax,dt,ISImin)

cnt = 0;
t = 0:dt:Tmax;
r = Freq/(length(t)-1)*Tmax/1000;       % [rate] = Spk/bin
spk = double(rand(1,length(t)+1)<r);    % Generated spike trains
spk = spk(1:length(t));
tspkbase = find(spk'>0)*dt;                 % Spike times
ISIbase = diff(tspkbase);
ISI = zeros(1);
i=0;
for j=1:length(ISIbase)
    if ISIbase(j) > ISImin
        i=i+1;
        ISI(i) = ISIbase(j);
    else
        cnt = cnt+ISIbase(j);
    end
end    
tspk = zeros(1);
tspk(1) = tspkbase(1);
for i=2:length(ISI)
    tspk(i) = tspk(i-1)+ISI(i);
end
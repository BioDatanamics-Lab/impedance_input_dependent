function [tspk,ISI] = Poisson(Freq,Tmax,dt)

t = 0:dt:Tmax;
r = Freq/(length(t)-1)*Tmax/1000;          % [rate] = Spk/bin
spk = double(rand(1,length(t))<r);    % Generated spike trains
tspk = find(spk'>0)*dt;                 % Spike times
ISI = diff(tspk);   
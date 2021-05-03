function [tspk,ISI] = Uniform(Freq,Tmax,dt)

t = 0:dt:Tmax;

ISI = 1000/Freq;
tspk = ISI:ISI:Tmax;
spk = zeros(1,length(t));
for j=1:length(tspk)
    spk(floor(tspk/dt))=1;
end
    

function [CycleFreqOrd,vmaxord] = OrderingCycleFreqPeak(CycleFreq,vmax)

CycleFreqOrd = CycleFreq;    
vmaxord = vmax;

L = 0;
for j=1:length(CycleFreqOrd)
    L = L+1;
    for k=L:length(CycleFreqOrd)
        if CycleFreqOrd(k)<CycleFreqOrd(j)          
            aux1 = CycleFreqOrd(j);
            CycleFreqOrd(j) = CycleFreqOrd(k);
            CycleFreqOrd(k) = aux1;
            aux2 = vmaxord(j);
            vmaxord(j) = vmaxord(k);
            vmaxord(k) = aux2;   
        end
    end
end
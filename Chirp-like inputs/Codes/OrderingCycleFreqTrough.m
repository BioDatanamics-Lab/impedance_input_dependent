function [CycleFreqOrd,vminord] = OrderingCycleFreqTrough(CycleFreq,vmin)

CycleFreqOrd = CycleFreq;    
vminord = vmin;

L = 0;
for j=1:length(CycleFreqOrd);
    L = L+1;
    for k=L:length(CycleFreqOrd)
        if CycleFreqOrd(k)<CycleFreqOrd(j)
            aux1 = CycleFreqOrd(j);
            CycleFreqOrd(j) = CycleFreqOrd(k);
            CycleFreqOrd(k) = aux1;  
            aux4 = vminord(j);
            vminord(j) = vminord(k);
            vminord(k) = aux4;
        end
    end
end
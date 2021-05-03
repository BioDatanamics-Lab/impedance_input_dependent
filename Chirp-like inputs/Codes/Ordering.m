A=[3 1 4 2]

L=0;
for j=1:length(A)
    L=L+1;
    for k=L:length(A)
        if A(k)>A(j)
            aux=A(j);
            A(j)=A(k);
            A(k)=aux;
        end
    end
end

A

L=0;
for j=1:length(A)
    L=L+1;
    for k=L:length(A)
        if A(k)<A(j)
            aux=A(j);
            A(j)=A(k);
            A(k)=aux;
        end
    end
end

A
          
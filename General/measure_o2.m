function o = measure_o2(mpsA,o1,o2,N)
% Measure a two site operator across of the system.

o=zeros([1,N-1]);
rho=cell([1,N-1]);
rho{N}=ones([1,1]);



for i=N:-1:3
    rho{i-1}=ncon({mpsA{i},conj(mpsA{i}),rho{i}},{[-1 1 3],[-2 2 3],[1 2]});
end

for i=1:N-1
        o(i)=ncon({mpsA{i},mpsA{i+1},o1,o2,conj(mpsA{i}),conj(mpsA{i+1}),rho{i+1}},{[1 4 2],[4 8 6],[2 3],[6 7],[1 5 3], [5 9 7],[8 9]});
end
%o=o(1);
end


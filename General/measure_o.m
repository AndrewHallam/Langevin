function o = measure_o(mpsA,O,N)
% Measure a single site operator at each site of the system.

o=zeros([1,N]);
rho=cell([1,N-1]);
rho{N}=ones([1,1]);



for i=N:-1:2
    rho{i-1}=ncon({mpsA{i},conj(mpsA{i}),rho{i}},{[-1 1 3],[-2 2 3],[1 2]});
end

for i=1:N
        o(i)=ncon({mpsA{i},O,conj(mpsA{i}),rho{i}},{[1 4 2],[2 3],[1 5 3],[4 5]});
    end
%o=o(1);
end


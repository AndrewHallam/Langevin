function eta = measure_eta(mpsA,mpsB,N)
%
eta=1;
for i=1:N
    eta=ncon({mpsA{i},conj(mpsB{i}),eta},{[1 -1 3],[2 -2 3],[1,2]});
end
eta=trace(eta);
end


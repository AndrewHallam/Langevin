function SE = measure_entanglement(mpsA,n,N)
% Measure bipartitite entanglement entropy,SE between sites 1:n and n+1:N of
% the chain.



rhoR=1;
for i=N:-1:n+1
    rhoR=ncon({mpsA{i},conj(mpsA{i}),rhoR},{[-1 1 3],[-2 2 3],[1 2]});
end

lam=eig(rhoR);
SE=sum(arrayfun(@(x) -x*log(x),lam));

end


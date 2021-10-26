function SE = measure_entanglement_1site(mpsA,n,N)
% Measure bipartitite entanglement entropy,SE between sites 1:n and n+1:N of
% the chain.



rhoR=1;
for i=N:-1:n+1
    rhoR=ncon({mpsA{i},conj(mpsA{i}),rhoR},{[-1 1 3],[-2 2 3],[1 2]});
end
rhoL=1;
for i=1:n-1
     rhoL=ncon({mpsA{i},conj(mpsA{i}),rhoL},{[1 -1 3],[2 -2 3],[1 2]});
end
rho=ncon({rhoL,mpsA{n},conj(mpsA{n}),rhoR},{[1 2],[1 3 -1],[2 4 -2],[3 4]});

lam=eig(rho);
SE=sum(arrayfun(@(x) -x*log(x),lam))/log(2);

end


function mps = finite_init(chi,d,N)

mps=cell([1,N]);

chi2=min([chi,d]);
A1=rand([1,chi2,d])+1i*rand([1,chi2,d]);

M=ncon({A1,conj(A1)},{[2 -1 1],[2 -2 1]});
[u1,d1,~]=svd(M);

A1=ncon({A1, conj(u1)*pinv(sqrtm(d1))},{[-1 1 -3],[1 -2]});

mps{1}=A1;



for q=2:N-1
    chi1=chi2;
    chi2=min([d^q,d^(N-q),chi]);
    Atemp=rand([chi1,chi2,d])+1i*rand([chi1,chi2,d]);
    M=ncon({Atemp,conj(Atemp)},{[1 -1 2],[1 -2 2]});
    [u1,d1,~]=svd(M);
    A2=ncon({Atemp, conj(u1)*pinv(sqrtm(d1))},{[-1 1 -3],[1 -2]});
    mps{q}=A2;
end


AN=rand([chi2,1,d])+1i*rand([chi2,1,d]);
M=ncon({AN,conj(AN)},{[1 2 3],[1 2 3]});
AN=AN/sqrt(M);
mps{N}=AN;


end




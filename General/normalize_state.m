function mps = normalize_state(mps,N)
% Normalize the incoming mps.

C=1;
for q=1:N
    A=mps{q};
    [Dl,Dr,d] = size(A);
    A=ncon({C,A},{[-1 1],[1 -2 -3]});
    A=permute(A,[1 3 2]);
    [Q,C] = qr(reshape(A,d*Dl,Dr),0);
    AL = reshape(Q,Dl,d,Dr);
    AL=permute(AL,[1 3 2]);
    mps{q}=AL;
end

end


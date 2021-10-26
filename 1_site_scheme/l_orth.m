function [C,AL] = l_orth(A)

A=permute(A,[1,3,2]);
[Dl,d,Dr] = size(A);
[Q,C] = qr(reshape(A,d*Dl,Dr),0);
AL = reshape(Q,Dl,d,Dr);
 AL=permute(AL,[1,3,2]);
end


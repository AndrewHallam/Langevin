function [C,AR] = r_orth(A)
A=permute(A,[1,3,2]);
[Dl,d,Dr]=size(A);
[Q,R]=qr(reshape(A,Dl,d*Dr)',0);
C=R';
AR=reshape(Q',Dl,d,Dr);
AR=permute(AR,[1,3,2]);

end


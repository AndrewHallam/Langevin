function [AL,AC] = l_orth_2site(AA,Dmax)

% Take in two sites of an MPS and left orthogonalize the 1st MPS. 
trunctol=eps;

AA=permute(AA,[1 3 4 2]);

[chi1,d1,d2,chi2]=size(AA);
D = min([Dmax, chi1*d1, d2*chi2]);

[U,S,V]=svd(reshape(AA,chi1*d1,d2*chi2),'econ');
S = diag(S);
normS = norm(S);

while norm(S(1:D-1))/normS > 1-trunctol
    D = D - 1;
end
truncerr = 1-norm(S(1:D))/normS;

AL = reshape(U(:,1:D),chi1,d1,D);
AC = reshape(diag(S(1:D))*V(:,1:D)',D,d2,chi2);

AL=permute(AL,[1 3 2]);
AC=permute(AC,[1 3 2]);

end


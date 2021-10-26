function ACout = dmrgupdate_1site(AC,FL,FR,H)
% Find the lowest eigenvalue of Hn in order to update the mps with a better
% estimation of the ground state.


% First reshape AC into a vector and Hn into a matrix.
[chi1,chi2,d]=size(AC);

AC2=reshape(AC,[chi1*chi2*d,1]);

% The initial guess for the eigenvalue problem is AC.
eigopts.v0=AC2;
%eigopts.p=min(chi1*chi2*d,20);

% Solve the eigenvalue problem.
[ACout,e]=eigs(@(x)(apply_H_to_A(x,FL,FR,H)),chi1*chi2*d,  1, 'SA', eigopts);

ACout=reshape(ACout,[chi1, chi2, d]);


end


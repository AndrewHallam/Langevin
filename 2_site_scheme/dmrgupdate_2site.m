function AACout = dmrgupdate_2site(AAC,FL,FR,H1,H2)
% Find the lowest eigenvalue of the 2-site Hn in order to update the mps with a better
% estimation of the ground state.


% First reshape AAC into a vector and Hn into a matrix.

[chi1,chi2,d1,d2]=size(AAC);
AAC2=reshape(AAC,[chi1*chi2*d1*d2,1]);


% The initial guess for the eigenvalue problem is AAC.

eigopts.v0=AAC2;
%eigopts.p=min(chi1*chi2*d^2,20);

% Solve the eigenvalue problem.

[AACout,e]=eigs(@(x)(apply_HH_to_AA(x,FL,FR,H1,H2)), chi1*chi2*d1*d2, 1, 'SA', eigopts);

AACout=reshape(AACout,[chi1, chi2 d1,d2]);

end


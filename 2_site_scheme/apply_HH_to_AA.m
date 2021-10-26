function bb = apply_HH_to_AA(aa,FL,FR,H1,H2)
% Apply the 2site effective Hamiltonian to the two center sites AA. AA is passed
% in as a vector, aa so it can be used for functions such as eigs and lanczos_exp. 
% This function is used so we never explicitly have to store the effective
% Hamiltonian in memory.

chi1 = size(FL, 1);
d1 = size(H1, 1);
d2 = size(H2, 1);

chi2 = size(FR, 1);

AA = reshape(aa, [chi1, chi2, d1,d2]);

% Contract tensor network for computing Heff * AA
HA=ncon({FL,AA,H1,H2,FR},{[1 2 -1],[1 6 3 5],[3 -3 2 4],[5 -4 4 7],[6 7 -2]});




bb = reshape(HA, size(aa));



end


function d = apply_K_to_C(c,FL,FR)
% Apply the 1site effective Hamiltonian to the center site A. A is passed
% in as a vector, a so it can be used for functions such as eigs. 
% This function is used so we never explicitly have to store the effective
% Hamiltonian in memory.

chi1 = size(FL, 1);
chi2 = size(FR, 1);

C = reshape(c, [chi1, chi2]);

% Contract tensor network for computing Heff * A
HA=ncon({FL,C,FR},{[1 3 -1],[1 2],[2 3 -2]});




d = reshape(HA, size(c));



end


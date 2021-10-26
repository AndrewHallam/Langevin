function b = apply_H_to_A(a,FL,FR,H)
% Apply the 1site effective Hamiltonian to the center site A. A is passed
% in as a vector, a so it can be used for functions such as eigs. 
% This function is used so we never explicitly have to store the effective
% Hamiltonian in memory.

chi1 = size(FL, 1);
d = size(H, 1);
chi2 = size(FR, 1);

A = reshape(a, [chi1, chi2, d]);

% Contract tensor network for computing Heff * A
HA=ncon({FL,A,H,FR},{[1 2 -1],[1 4 3],[3 -3 2 5],[4 5 -2]});




b = reshape(HA, size(a));



end


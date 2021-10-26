function FR = update_right_env(A,H,FR)
% The right environment is the MPS A contracted with the Hamiltonian H to
% the right of the site we care about.
       
[c1,~,~]=size(A);
[~,~,c2,~]=size(H);

FR=ncon({A,H,conj(A),FR},{[-1 1 2],[2  5 -2 3],[-3 4 5],[1 3 4]});
FR=reshape(FR,[c1,c2,c1]);
end


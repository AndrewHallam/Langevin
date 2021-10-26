function mpsA = dmrg_finite(mpsA,H,N)

%% DMRG algorithm for N sites, using the single site algorithm. 
% This algorithm takes in an MPS and outputs and MPS of the same bond dimension
% bond dimension. Dmax is the maximum possible bond dimension. 

% mpsA is the mps being evolved. 
% H is the Hamiltonian in MPO form. 

% N is the size of the system. 



% It begins by defining the environment (the MPS and Hamiltonian contracted
% together) from the left all the way to the first site. 
F=cell([1,N+2]);
F{1}=ones([1,1,1]);

for q=1:N
    F{q+1}=update_left_env(mpsA{q},H{q},F{q});
end
F{N+2}=ones([1,1,1]);


%% SWEEP RIGHT
% Start from site N and sweep across.
CA=1;

for n=N:-1:2
    % construct center site by contracting A{n} with C
    AC = ncon({mpsA{n},CA},{[-1,1,-3],[1,-2]});

    
    % find the lowest energy eigenvalue of Hn
    AC=dmrgupdate_1site(AC,F{n},F{n+2},H{n});
    
    % Find the right-canonical MPS.
    
    [CA,AR]=r_orth(AC);
    
    mpsA{n}=AR;
    
    %Update the environment.
    F{n+1}=update_right_env(AR,H{n},F{n+2});
    
    

    
end

%% SITE NUMBER 1:

n = 1;
    AC = ncon({mpsA{n},CA},{[-1,1,-3],[1,-2]});

    

    AC=dmrgupdate_1site(AC,F{n},F{n+2},H{n});

    
    % We now sweep from left to right so we left orthogonalize.
    [CA,AL]=l_orth(AC);
    
    mpsA{n}=AL;

    
    F{n+1}=update_left_env(AL,H{n},F{n});




%% SWEEP LEFT
% We now sweep from the left to the right.
for n=2:N

    

    AC = ncon({CA,mpsA{n}},{[-1,1],[1,-2,-3]});
    % Find the 1-site effective Hamiltonian.

    AC=dmrgupdate_1site(AC,F{n},F{n+2},H{n});

    % We now sweep from left to right so we left orthogonalize.
    [CA,AL]=l_orth(AC);

    
    mpsA{n}=AL;


    F{n+1}=update_left_env(AL,H{n},F{n});


end


end


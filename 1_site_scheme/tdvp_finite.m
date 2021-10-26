function mpsA = tdvp_finite(mpsA,H,dt,N,evolution)

%% TDVP algorithm for N sites, using the single site algorithm. 
% This algorithm takes in an MPS and outputs and MPS of the same bond dimension
%bond dimension. Dmax is the maximum possible bond dimension. 

% mpsA is the mps being evolved. 
% H is the Hamiltonian in MPO form. 
% dt is the time step. 
% N is the size of the system. 
% evolution can be 'real' or 'imag'. 'imag' is used for finding ground
% states. 


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

%     % Find the effective Hamiltonian.
%     Hn=find_Hn(H{n},F{n},F{n+2});
%     % forward evolution of center site: half time step
% 
%     [AC]=apply_Hn(AC,Hn,dt/2,evolution);
    

    AC=update_A(AC,F{n},F{n+2},H{n},dt/2,evolution);
        
    % Find the right-canonical MPS.
    
    [CA,AR]=r_orth(AC);
    
    mpsA{n}=AR;
    
    %Update the environment.
    F{n+1}=update_right_env(AR,H{n},F{n+2});
    
    
    CA=update_C(CA,F{n},F{n+1},-dt/2,evolution);

    
end

%% SITE NUMBER 1:

% We treat the first site in a special way, evolving it a full time step. 
n = 1;
    AC = ncon({mpsA{n},CA},{[-1,1,-3],[1,-2]});

    
    
    % forward evolution of center site: full time step

    AC=update_A(AC,F{n},F{n+2},H{n},dt,evolution);

    
    
    [CA,AL]=l_orth(AC);
    
    mpsA{n}=AL;

    
    F{n+1}=update_left_env(AL,H{n},F{n});




%% SWEEP LEFT
% We now sweep from the left to the right.
for n=2:N
    
    CA=update_C(CA,F{n},F{n+1},-dt/2,evolution);
    
    

    AC = ncon({CA,mpsA{n}},{[-1,1],[1,-2,-3]});
    
    AC=update_A(AC,F{n},F{n+2},H{n},dt/2,evolution);

    % Split the two sites into a left orthogonal AL and AC.
    [CA,AL]=l_orth(AC);

    
    mpsA{n}=AL;


    F{n+1}=update_left_env(AL,H{n},F{n});


end


end


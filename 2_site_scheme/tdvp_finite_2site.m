function [mpsA,Dmaxout] = tdvp_finite_2site(mpsA,H,dt,N,Dmax,evolution)

%% 2-site TDVP algorithm for N sites.
% This algorithm takes in an MPS and outputs and MPS of (potentially)
% higher bond dimension. Dmax is the maximum possible bond dimension. Dmaxout 
% is the maximum value of the bond dimension of the outgoing mps.
% 
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

Dmaxout=1;
%% SWEEP RIGHT
% Start from site N and sweep across.

AC=mpsA{N};
for n=N-1:-1:2
        % construct 2-site center by contracting A{n} with AC

    AAC=ncon({mpsA{n},AC},{[-1 1 -3],[1 -2 -4]});
        
    
    AAC=update_AA(AAC,F{n},F{n+3},H{n},H{n+1},dt/2,evolution);

    
    % Split the two sites into an AC and a right orthogonal AR.
    [AC,AR]=r_orth_2site(AAC,Dmax);
    
    mpsA{n+1}=AR;
    
    % Update the environments. 
    F{n+2} = update_right_env(AR,H{n+1},F{n+3});
    
    % We then backwards evolve the single site AC. This stops us doing
    % twice as much time evolution to some sites by evolving two at a time.
    
    AC=update_A(AC,F{n},F{n+2},H{n},-dt/2,evolution);
    
    Dmaxout=max(Dmaxout,size(AR,1));
    
    
end

%% SITE NUMBER 1 and 2:

% We treat the first and second sites in a special way, evolving them a full time step. 

n = 1;

AAC = ncon({mpsA{n},AC},{[-1,1,-3],[1,-2,-4]});

AAC=update_AA(AAC,F{n},F{n+3},H{n},H{n+1},dt,evolution);




[AL,AC]=l_orth_2site(AAC,Dmax);

mpsA{n}=AL;

F{n+1}=update_left_env(AL,H{n},F{n});
    Dmaxout=max(Dmaxout,size(AL,2));




%% SWEEP LEFT
% We now sweep from the left to the right.

for n=2:N-1
    
    
    AC=update_A(AC,F{n},F{n+2},H{n},-dt/2,evolution);
    

    AAC=ncon({AC,mpsA{n+1}},{[-1 1 -3],[1 -2 -4]});

    AAC=update_AA(AAC,F{n},F{n+3},H{n},H{n+1},dt/2,evolution);
    
    
    [AL,AC]=l_orth_2site(AAC,Dmax);
    
    mpsA{n}=AL;
    
    
    F{n+1}=update_left_env(AL,H{n},F{n});
        Dmaxout=max(Dmaxout,size(AL,2));

    
end

[CA, AL] = l_orth(AC);
mpsA{N} = AL;





end


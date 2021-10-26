function [mpsA,Dmaxout] = dmrg_finite_2site(mpsA,H,N,Dmax)

%% 2-site DMRG algorithm for N sites.
% This algorithm takes in an MPS and outputs and MPS of (potentially)
% higher bond dimension. Dmax is the maximum possible bond dimension. Dmaxout 
% is the maximum value of the bond dimension of the outgoing mps.
% 
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

Dmaxout=1;
%% SWEEP RIGHT
% Start from site N and sweep across.

AC=mpsA{N};
for n=N-1:-1:2
        % construct 2-site center by contracting A{n} with AC

    AAC=ncon({mpsA{n},AC},{[-1 1 -3],[1 -2 -4]});
        
    AAC=dmrgupdate_2site(AAC,F{n},F{n+3},H{n},H{n+1});
    
    % Split the two sites into an AC and a right orthogonal AR.
    [AC,AR]=r_orth_2site(AAC,Dmax);
    
    mpsA{n+1}=AR;
    
    % Update the environments. 
    F{n+2} = update_right_env(AR,H{n+1},F{n+3});
    
    
    Dmaxout=max(Dmaxout,size(AR,1));
    
    
end

%% SITE NUMBER 1 and 2:


n = 1;

AAC = ncon({mpsA{n},AC},{[-1,1,-3],[1,-2,-4]});


      AAC=dmrgupdate_2site(AAC,F{n},F{n+3},H{n},H{n+1});



[AL,AC]=l_orth_2site(AAC,Dmax);

mpsA{n}=AL;

F{n+1}=update_left_env(AL,H{n},F{n});
    Dmaxout=max(Dmaxout,size(AL,2));




%% SWEEP LEFT
% We now sweep from the left to the right.

for n=2:N-1
    
    
    
     % construct 2-site center by contracting AC with A{n+1}

    AAC=ncon({AC,mpsA{n+1}},{[-1 1 -3],[1 -2 -4]});

        % Find the 2-site effective Hamiltonian.

     AAC=dmrgupdate_2site(AAC,F{n},F{n+3},H{n},H{n+1});

    
    
    [AL,AC]=l_orth_2site(AAC,Dmax);
    
    mpsA{n}=AL;
    
    
    F{n+1}=update_left_env(AL,H{n},F{n});
        Dmaxout=max(Dmaxout,size(AL,2));

    
end

[CA, AL] = l_orth(AC);
mpsA{N} = AL;





end


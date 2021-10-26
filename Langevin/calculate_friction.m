function friction = calculate_friction(mpsA,H,Bk,gamma)
% This calculates the friction terms in the Langevin equation.
% The output is a size K (K=3*N for isotropic noise) vector which looks
% like friction = gamma * inv(I+alpha*M)*p
% Where gamma is the coupling of the system to the bath.
% I is the identity (K by K). 
% M is the poisson bracket of the different noise fields M={f_i,f_j}
% p is the poisson bracket of the noise fields with H (H is the standard
% Hamiltonian WITH the noise terms added), p = {f_i,H+\sum_n f_n}.



N=size(mpsA,2);
K=size(Bk,1);


% We will run through the systema find the contributions to p and M on each
% site. 
p=zeros([K,1]);
M=zeros([K,K]);

FH=cell([1,N+2]);
FH{1}=ones([1,1,1]);






    % Run through the system from right to left. 

    % Calculate the environment terms for H (FH) and for the noise fields (FB).

for q=1:N

    
% We start by calculating the AC contributions to the poisson bracket.

FH{q+1}=update_left_env(mpsA{q},H{q},FH{q});
end
FH{N+2}=ones([1,1,1]);

FB=cell([K,N+2]);
for k=1:K
FB{k,1}=ones([1,1,1]);

for q=1:N
FB{k,q+1}=update_left_env(mpsA{q},Bk{k,q},FB{k,q});
end
FB{k,N+2}=ones([1,1,1]);

end


CA=1;

for n=N:-1:2
p_temp_AC=zeros([K,1]);
M_temp_AC=zeros([K,K]);
p_temp_C=zeros([K,1]);
M_temp_C=zeros([K,K]);

% construct center site by contracting A{n} with C
AC = ncon({mpsA{n},CA},{[-1,1,-3],[1,-2]});

% Find the effective Hamiltonian applied to A.
Hn=apply_H_to_A(AC,FH{n},FH{n+2},H{n});

% Find the effective Hamiltonian applied to A for the noise fields. 
%I can probably restrict this to only calculate some of these but I
%haven't done this.

for k=1:K
Bn{k}=apply_H_to_A(AC,FB{k,n},FB{k,n+2},Bk{k,n});
end


% Only calculate the poisson bracket for terms within n+/-lengthscale. the
% function looks strange because we assume 3 noise fields per site.
for k=1:K
p_temp_AC(k)=poisson(Hn,Bn{k});


% Same for calculating M.

% M is necessarily antisymmetric so only calculate the upper triangular
% part.

% Since M is made up of local noise terms {f_i,f_j} = 0 if the fields
% aren't on the same sites. Therefore we only need to calculate a handful
% of terms along the diagonal.
for l=k+1:K
if  abs(floor((k-1)/3))==abs(floor((l-1)/3)) 
M_temp_AC(k,l)=poisson(Bn{k},Bn{l});
else
 M_temp_AC(k,l)=0;
end
end
end
M_temp_AC=M_temp_AC-M_temp_AC.';

[CA,AR]=r_orth(AC);

mpsA{n}=AR;

%Update the environment.
FH{n+1}=update_right_env(AR,H{n},FH{n+2});
for k=1:K
FB{k,n+1}=update_right_env(AR,Bk{k,n},FB{k,n+2});
end


% There are also contributions to the Poisson bracket due to the C parts.
% These are very similar to the AC part.

Kn=apply_K_to_C(CA,FH{n},FH{n+1});


for k=1:K
    Cn{k}=apply_K_to_C(CA,FB{k,n},FB{k,n+1});

end
for k=1:K
p_temp_C(k)=poisson(Kn,Cn{k});

%abs(floor((k-1)/3))==abs(floor((l-1)/3)) &&
for l=k+1:K
if  abs(floor((k-1)/3))==abs(floor((l-1)/3)) 
M_temp_C(k,l)=poisson(Cn{k},Cn{l});
else
    M_temp_C(k,l)=0;
end
end
end

M_temp_C=M_temp_C-M_temp_C.';

p=p-p_temp_C+p_temp_AC;
M=M-M_temp_C+M_temp_AC;


end

%% SITE NUMBER 1:

% We then do the 1st site. There are only AC parts here since we don't
% change the norm of the state with just H + noise.
n = 1;
p_temp_AC=zeros([k,1]);
M_temp_AC=zeros([k,k]);

AC = ncon({mpsA{n},CA},{[-1,1,-3],[1,-2]});



Hn=apply_H_to_A(AC,FH{n},FH{n+2},H{n});


for k=1:K
Bn{k}=apply_H_to_A(AC,FB{k,n},FB{k,n+2},Bk{k,n});

end

for k=1:K
p_temp_AC(k)=poisson(Hn,Bn{k});

for l=k+1:K
if abs(floor((k-1)/3))==abs(floor((l-1)/3))
M_temp_AC(k,l)=poisson(Bn{k},Bn{l});
else
 M_temp_AC(k,l)=0;
end
end
end

M_temp_AC=M_temp_AC-M_temp_AC.';

p=p+p_temp_AC;
M=M+M_temp_AC;

% Finally we calculate the friction term. Using A\b to do matrix inversion
% might not be optimal? 

friction=-1i*gamma*((eye(K)-1i*gamma*M)\p);


end


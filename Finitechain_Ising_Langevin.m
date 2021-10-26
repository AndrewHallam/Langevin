%clear;
close all;


% This code will calculate the time evolution of the transverse and longitudinal field 
% Ising model, quenching from some values for the transverse and longitudinal fields to
% some others. The evolution will be done using Langevin TDVP with noise
% and dissipation. 

N=30; % System size
D=250; % D is the maximum possible bond dimension of our MPS.

Tmax=10  ; % Maximum time.
Num_timestep=500; % Number of timesteps.

Num_dmrg=20; % Number of steps to take in dmrg algorithm to find initial state.


% These parameters determine the initial state we quench from.
hx_init=-1;
J_init=0;
hy_init=0;
hz_init=0;  

% These parameters determine the Hamiltonian we time evolve with.
hx=-1.05;
J=-1.00;
hz=0.5;
hy=0.0;


% These parameters are used for the Langevin dynamics. Num_iter determines
% how many noise iterations to consider. gamma determines the strength of
% the friction and T is the strength of the bath.
Num_iter=1;
gamma=0.0;
T=0.00;

% If you want to work in the dephasing limit set dephasing to be true.
% Otherwise set it false. 
% In the dephasing limit dF/dt is not calculated so the code will be fast.
dephasing=true;
gammaT=0.0;



%The Ising hamiltonian is defined using the following MPO construction.
H=ising_mpo(N,J_init,hx_init,hy_init,hz_init);

X = [0 1;1 0]; Y = [0 -1i;1i 0]; Z = [1 0;0 -1]; I = eye(2); P = (I-Z)/2;


% Define all the noise terms.

bath=cell([1,3]);
bath{1}=X;
bath{2}=Y;
bath{3}=Z;


K=3*N;
Bk=cell([3*N,N]);
H_b=cell([1,3]);
for q=1:3
H_temp=zeros([2,2,1,1]);
H_temp(:,:,1,1)=bath{q};
H_b{q}=H_temp;
end
H_I=zeros([2,2,1,1]);
H_I(:,:,1,1)=eye(2);

for k=1:N
for n=1:N
if n==k
    for k2=1:3
        
Bk{k2+3*(k-1),n}=H_b{k2};

    end
else
    for k2=1:3
Bk{k2+3*(k-1),n}=H_I;


    end

end

end
end



%The MPS is also defined in terms of a cell, one MPS for each site.
% I choose MPS to be defined so the 1st dimension of the array is the left
% virtual leg, the 2nd dimesion is the right virtual leg and the final
% dimension is the physical bond dimension. 

% We initialize the MPS randomly.
mpsA_init=finite_init(1,2,N);

disp('Finding ground state of initial Hamiltonian')
% We then use DMRG to find the ground state of our initial Hamiltonian.
Dout=D;
for n=1:Num_dmrg
% If the bond dimension isn't high enough let it increase with the
% 2-site scheme. Otherwise use the 1 site scheme. 
if Dout<D
[mpsA_init,Dout]=dmrg_finite_2site(mpsA_init,H,N,D);
Dout;
else
[mpsA_init]=dmrg_finite(mpsA_init,H,N);
end

mpsA_init=normalize_state(mpsA_init,N);
end

disp('Finding the time evolution')

% We now do the real time evolution.

% Define the Hamiltonian we will time evolve with.
H=ising_mpo(N,J,hx,hy,hz);


measureZ=zeros([Num_iter,Num_timestep+1,N]);
measureZZ=zeros([Num_iter,Num_timestep+1,N-1]);
measureX=zeros([Num_iter,Num_timestep+1,N]);
measureXX=zeros([Num_iter,Num_timestep+1,N-1]);
measureY=zeros([Num_iter,Num_timestep+1,N]);
measureSE=zeros([Num_iter,Num_timestep+1]);
measureE=zeros([Num_iter,Num_timestep+1]);

tic
for q=1:Num_iter
mpsA=mpsA_init;



% All this is for doing plots. 
sub=subplot(2,2,1);
fplot=animatedline('Marker','.','Color','r');
xlim([0,Tmax]);
title('sigmaZ - sigmaZ correlations')
sub=subplot(2,2,2);
xlim([0,Tmax]);
title('sigmaX - sigmaX correlations')

gplot=animatedline('Marker','.','Color','k');
sub=subplot(2,2,3);
xlim([0,Tmax]);
title('Energy')

hplot=animatedline('Marker','.','Color','g');
sub=subplot(2,2,4);
xlim([0,Tmax]);
title('Entanglement')

jplot=animatedline('Marker','.','Color','c');

xlim([0,Tmax]);
dt=Tmax/Num_timestep;
mpsA_init;


Zo=measure_o(mpsA,Z,N);
Yo=measure_o(mpsA,Y,N);

Xo=measure_o(mpsA,X,N);

ZZo=measure_o2(mpsA,Z,Z,N);
XXo=measure_o2(mpsA,X,X,N);
SE=measure_entanglement(mpsA,round(N/2),N);
E=measure_mpo(mpsA,H,N);




addpoints(fplot,0*dt,real(mean(ZZo)));
addpoints(gplot,0*dt,real(mean(XXo)));
addpoints(hplot,0*dt,real(E));
addpoints(jplot,0*dt,abs(SE));

measureZ(q,1,:)=real(Zo);
measureZZ(q,1,:)=real(ZZo);
measureX(q,1,:)=real(Xo);
measureXX(q,1,:)=real(XXo);
measureY(q,1,:)=real(Yo);
measureSE(q,1)=abs(SE);
measureE(q,1)=real(E);

etak=zeros([K,1]);

Dout=1;
for n=1:Num_timestep
% If the bond dimension isn't high enough let it increase with the
% 2-site scheme. Otherwise use the 1 site scheme. 
for k=1:K
if dephasing
    etak(k)=normrnd(0,sqrt(dt*2*gammaT));
    % in dephasing limit gamma is zero.
    gamma=0;
else
etak(k)=normrnd(0,sqrt(dt*2*gamma*T));
end
end

if Dout<D
[mpsA,Dout]=tdvp_langevin_finite_2site(mpsA,H,dt,N,D,'real',Bk,etak,gamma);
Dout;
else
[mpsA]=tdvp_langevin_finite(mpsA,H,dt,N,'real',Bk,etak,gamma);
end

%Measure observables and plot them.
mpsA=normalize_state(mpsA,N);
Zo=measure_o(mpsA,Z,N);
Xo=measure_o(mpsA,X,N);
Yo=measure_o(mpsA,Y,N);

ZZo=measure_o2(mpsA,Z,Z,N);
XXo=measure_o2(mpsA,X,X,N);
SE=measure_entanglement(mpsA,round(N/2),N);
E=measure_mpo(mpsA,H,N);




measureZ(q,n+1,:)=real(Zo);
measureY(q,n+1,:)=real(Yo);
measureZZ(q,n+1,:)=real(ZZo);
measureX(q,n+1,:)=real(Xo);
measureXX(q,n+1,:)=real(XXo);
measureSE(q,n+1)=abs(SE);
measureE(q,n+1)=real(E);

addpoints(fplot,n*dt,real(mean(ZZo)));

addpoints(gplot,n*dt,real(mean(XXo)));
addpoints(hplot,n*dt,real(E));
addpoints(jplot,n*dt,abs(SE));

drawnow limitrate

end
end






toc
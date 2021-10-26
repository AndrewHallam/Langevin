function [AAdt] = update_AA(AA,FL,FR,H1,H2,dt,time)
% Apply the 2site effective Hamiltonian to the two center sites AA.
% The time step is different for real and imaginary evolutions.


if strcmpi(time,'real')
    DTA=-1i*dt;
    
elseif strcmpi(time,'imag')
    DTA=-1*dt;
    
elseif strcmpi(time,'gradient')
    DTA=-1*dt;
    
end

%Reshape AA to a vector.
[chi1,chi2,d1,d2]=size(AA);
AA2=reshape(AA,[chi1*chi2*d1*d2,1]);


%Using lanczos to calculate the matrix exponential effectively.
AAdt=lanczos_exp(@(x)(apply_HH_to_AA(x,FL,FR,H1,H2)),AA2,DTA);

AAdt=reshape(AAdt,[chi1 chi2 d1 d2]);

end


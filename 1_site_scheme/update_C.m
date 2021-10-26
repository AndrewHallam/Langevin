function [CAdt] = update_C(CA,FL,FR,dt,time)
% Apply the o-site effective Hamiltonian to C.
% The time step is different for real and imaginary evolutions.
%


if strcmpi(time,'real')
    DTA=-1i*dt;
    
elseif strcmpi(time,'imag')
    DTA=-1*dt;
    
elseif strcmpi(time,'gradient')
    DTA=-1*dt;
    
end
%Reshape C to a vector and Kn to a matrix. 
[chi1,chi2]=size(CA);

CA=reshape(CA,[chi1*chi2,1]);

% Use lanczos to calculate the matrix exponential.

CAdt=lanczos_exp(@(x)(apply_K_to_C(x,FL,FR)),CA,DTA);

% Reshape back.
CAdt=reshape(CAdt,[chi1 chi2]);


end


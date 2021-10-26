function [ACdt] = update_A(AC,FL,FR,H,dt,time)
% Apply the 1site effective Hamiltonian to the center site AC.
% The time step is different for real and imaginary evolutions.
if strcmpi(time,'real')
    DTA=-1i*dt;
    
elseif strcmpi(time,'imag')
    DTA=-1*dt;
    
elseif strcmpi(time,'gradient')
    DTA=-1*dt;
    
end

% We reshape AC to a vector.
[chi1,chi2,d]=size(AC);
AC2=reshape(AC,[chi1*chi2*d,1]);

% We use a lanczos method of calculating the matrix exponential because it
% more efficient.

ACdt=lanczos_exp(@(x)(apply_H_to_A(x,FL,FR,H)),AC2,DTA);

% Make AC an MPS again.
ACdt=reshape(ACdt,[chi1 chi2 d]);



end


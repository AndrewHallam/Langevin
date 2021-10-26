function p = poisson(f1,f2)
% This calculates the poisson bracket between two operators on a single
% site of the system. f1 is the effective hamiltonian for one operator F1
% applied to AC. f2 is the effective hamiltonian for the other operator F2
% applied to AC.
f1=reshape(f1,[numel(f1),1]);     
f2=reshape(f2,[numel(f1),1]);     

p=(f2'*f1-f1'*f2);
end


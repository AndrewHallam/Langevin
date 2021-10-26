function FL = update_left_env(A,H,FL)
       
% The left environment is the MPS A contracted with the Hamilotnian H to
% the left of the site we care about.

FL=ncon({FL,A,H,conj(A)},{[1 3 4],[1 -1 2],[2 5 3 -2],[4 -3 5]});

end


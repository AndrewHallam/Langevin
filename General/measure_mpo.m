function o = measure_mpo(mpsA,O,N)
%
o=1;
    for i=1:N
        o=ncon({o,mpsA{i},O{i},conj(mpsA{i})},{[1 3 5],[1 -1 2],[2 4 3 -2],[5 -3 4]});
    end
o=real(o);
end


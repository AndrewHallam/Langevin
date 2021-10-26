function H = ising_mpo(N,J,hx,hy,hz)
X = [0 1;1 0]; Y = [0 -1i;1i 0]; Z = [1 0;0 -1]; 

% H right and H left define the beginning and end of the MPO. 

H_right=zeros([2,2,3,1]);
H_right(:,:,1,1)=eye(2);
H_right(:,:,2,1)=Z;
H_right(:,:,3,1)=hx*X+hy*Y+hz*Z;
H_left=zeros([2,2,1,3]);
H_left(:,:,1,1)=hx*X+hy*Y+hz*Z;
H_left(:,:,1,2)=J*Z;
H_left(:,:,1,3)=eye(2);

% This defines the MPO in the bulk.
H_bulk=zeros([2,2,3,3]);
H_bulk(:,:,1,1)=eye(2);
H_bulk(:,:,2,1)=Z;
H_bulk(:,:,3,1)=hx*X+hy*Y+hz*Z;
H_bulk(:,:,3,2)=J*Z;
H_bulk(:,:,3,3)=eye(2);


% A cell lets me collect all the arrays together in MATLAB. Each element of
% the array corresponds to one site.
H=cell([1,N]);

H{1}=H_left;
H{N}=H_right;
for n=2:N-1
H{n}=H_bulk;
end
end


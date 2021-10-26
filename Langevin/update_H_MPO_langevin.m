function H = update_H_MPO_langevin(H,etak,friction)
K=size(etak,1);
N=size(H,2);
X = [0 1;1 0]; Y = [0 -1i;1i 0]; Z = [1 0;0 -1]; I = eye(2);

if nargin==2
friction=zeros([K,1]);
end

HL=H{1};

HL(:,:,1,1)=HL(:,:,1,1)+(etak(1)+friction(1))*X+(etak(2)+friction(2))*Y+(etak(3)+friction(3))*Z;

H{1}=HL;


HR=H{end};

HR(:,:,end,1)=HR(:,:,end,1)+(etak(end-2)+friction(end-2))*X+(etak(end-1)-1i*friction(end-1))*Y+(etak(end)+friction(end))*Z;

H{end}=HR;

for n=2:N-1
Htemp=H{n};
nk=3*(n-1);
Htemp(:,:,end,1)=Htemp(:,:,end,1)+(etak(nk+1)+friction(nk+1))*X+(etak(nk+2)+friction(nk+2))*Y+(etak(nk+3)+friction(nk+3))*Z;

H{n}=Htemp;
end

end


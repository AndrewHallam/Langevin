function vout = lanczos_exp(H,b,dt)
%
n=size(b,1);
if n>1
m=20;
m=min(m,n);
b=b/norm(b);
v(:,1)=b;
wtemp=H(b);

alpha(1)=wtemp'*v(:,1);

w(:,1)=wtemp-alpha(1)*v(:,1);
for q=2:m
    beta(q-1)=norm(w(:,q-1));
    v(:,q)=w(:,q-1)/beta(q-1);
    wtemp=H(v(:,q));
    alpha(q)=wtemp'*v(:,q);
    w(:,q)=wtemp-alpha(q)*v(:,q)-beta(q-1)*v(:,q-1);
    if abs(beta(q-1))<eps
        m=q-1;
        v=v(:,1:m);
        alpha=alpha(1:m);
        beta=beta(1:m-1);
        break
    end
end
T=diag(beta,1)+diag(beta,1).'+diag(alpha);
basis=zeros([m,1]);
basis(1)=1;
vout=norm(b)*v*expm(dt*T)*basis;


else
    h=H(b)/b;
    vout=exp(h*dt)*b;
end
end


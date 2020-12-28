function [L,U]=getLUinlierscale(X,branch,sigma)
k_max=branch(3);
k_min=branch(1);
b_min=branch(2);
b_max=branch(4);
x=X(1,:);
y=X(2,:);

e_max=k_max*x+b_max-y;
e_min=k_min*x+b_min-y;
e=[abs(e_max);abs(e_min)];
flag=(e_max.*e_min)>0;
r_max=max(e);
r_min=flag.*min(e);%emaxmin Í¬ºÅ

flag1=(sigma<=r_min);rou(flag1==1)=1;
flag2=(sigma>=r_max);rou(flag2==1)=0;
flag3=((sigma<=r_max)+(sigma>=r_min)==2);rou(flag3==1)=0;



L=sum(rou)./size(x,2);
x0=0.5*(branch(1:2)+branch(3:4));
e0=abs(x0(1)*x+x0(2)-y);
rou0((e0<=sigma))=0;
rou0((e0>sigma))=1;
U=sum(rou0)./size(x,2);

end
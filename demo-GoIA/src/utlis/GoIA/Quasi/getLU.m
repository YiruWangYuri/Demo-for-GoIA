function [L,U]=getLU(x,y,z,branch,sigma)
a_min=branch(1);
a_max=branch(4);
b_min=branch(2);
b_max=branch(5);
c_min=branch(3);
c_max=branch(6);

e_max=z-(a_min*x+b_min*y+c_min);
e_min=z-(a_max*x+b_max*y+c_max);

e=[abs(e_max);abs(e_min)];
flag=(e_max.*e_min)>0;
r_max=max(e);
r_min=flag.*min(e);%emaxmin Í¬ºÅ

rou = zeros(1,length(x));
flag1=sigma<=r_min;rou((flag1==1))=1;




L=sum(rou)./size(x,2);

x0=0.5*(branch(1:4)+branch(5:8));
e0=abs(x0(1)*x+x0(2)*y+x0(3)*z+x0(4));
rou0((e0<=sigma))=0;
rou0((e0>sigma))=1;
U=sum(rou0)./size(x,2);

end
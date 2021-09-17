clc
clear all

m=4000/386.22
k=800
W_n=8.7889
kisay=0.4
c=2*m*W_n*kisay
L=3 %in
V=10*1.466
fy=5000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=500;
t_total=5;
t_d=0.2;
delta_t=t_total./n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(1,n);
udot=zeros(1,n);
uddot=zeros(1,n);
t=zeros(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=(6./(delta_t.^2)).*m+(3./delta_t).*c;
a2=(6./delta_t).*m+2.*c;
a3=2.*m+(delta_t).*0.5.*c;
k_hat=k+a1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    j=i+1
    t(1,j)=t(1,i)+delta_t;   
       if t(1,i)<t_d
          p(1,j)=m.*6.*((pi.^2)/(t_d.^2)).*sin((pi./t_d).*t(1,j));
       else
          p(1,j)=0;
       end
       if i==1
          u(1,i)=0;udot(1,i)=0;p(1,i)=0;uddot(1,i)=0;
       end  
      p_hat(1,j)=p(1,j)+a1.*u(1,i)+a2.*udot(1,i)+a3.*uddot(1,i);
      u(1,j)=p_hat(1,j)./k_hat;
      udot(1,j)=(3./delta_t).*(u(1,j)-u(1,i))-2.*udot(1,i)-delta_t.*0.5.*uddot(1,i);
      uddot(1,j)=(6./(delta_t.^2)).*(u(1,j)-u(1,i))-(6./delta_t).*udot(1,i)-2.*uddot(1,i);
end

plot(t,u,'--');grid on; hold on
xlabel('time(sec)','FontSize',12);
ylabel(' max Displacement(in)','FontSize',12);


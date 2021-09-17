clc
clear all

m=4000/386.22
k=800
W_n=sqrt(k/m)
kisay=0.4
c=2*m*W_n*kisay
L=3 %in
V=10*1.466
fy=5000
beta=1/6;
gama=1/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=500;
t_total=5
t_d=L/V
delta_t=t_total./n;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(1,n);
udot=zeros(1,n);
uddot=zeros(1,n);
t=zeros(1,n);
k=zeros(1,n);
fs=zeros(1,n);
 a=m/(beta*delta_t^2)+(c*gama/(beta*delta_t));
a1=((m/(beta*delta_t))+gama*c/beta);
a2=(gama*c*delta_t/(2*beta))+(0.5*m/beta)-c*delta_t;

for i=1:n-1
    j=i+1;
    t(1,j)=t(1,i)+delta_t;
    if abs(fs(1,i))>5000
        k(1,i)=0.0001;
    else
        k(1,i)=800;   
    end
       k_hat=k(1,i)+a1;
       if t(1,i)<t_d
          p(1,j)=m.*6.*((pi.^2)/(t_d.^2)).*sin((pi./t_d).*t(1,j));
       else
          p(1,j)=0;
       end
       k_hat(1,i)=k(1,i)+a;
       delta_U(1,i)=(p(1,j)-p(1,i)+a1*udot(1,i)+a2*uddot(1,i))/k_hat(1,i);
       delta_Udot(1,i)=(1-0.5*gama/beta)*delta_t*uddot(1,i)+(gama/(beta*delta_t))*delta_U(1,i)-gama*udot(1,i)/beta;
       delta_Uddot(1,i)=(1/(beta*(delta_t^2)))*delta_U(1,i)-(1/(beta*delta_t))*udot(1,i)-0.5*uddot(1,i)/beta;
       u(1,j)=u(1,i)+delta_U(1,i);
       udot(1,j)=delta_Udot(1,i)+udot(1,i);
       fs(1,j)=fs(1,i)+k(1,i)*delta_U(1,i);
   if fs(1,j)>fy
       fs(1,j)=fy;
   elseif fs(1,j)<-fy
       fs(1,j)=-fy;
   end
    uddot(1,j)=(p(1,j)-c*udot(1,j)-fs(1,j))/m;
end
plot(t,u);grid on; hold on
xlabel('time(sec)','FontSize',12);
ylabel(' max Displacement(in)','FontSize',12);




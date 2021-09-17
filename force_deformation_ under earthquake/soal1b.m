clc
clear all
close all
m=4000/386.22;
k=0.8;
W_n=sqrt(k*1000/m);
kisay=0.4;
c=2*m*W_n*kisay;
L=3 ;
gama=0.5;
beta=1/6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_t=0.02;
t=0:delta_t:5;
v=1:1:65;
V=1.466*v;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=(m/(beta*(delta_t^2)))+(gama/(beta*delta_t))*c;
a2=(m/(beta*delta_t))+((gama/(beta))-1)*c;
a3=((1/(2*beta))-1)*m+delta_t*((gama/(2*beta))-1)*c;
k_hat=k+a1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(65,5/(delta_t)+1);
udot=zeros(65,5/(delta_t)+1);
uddot=zeros(65,5/(delta_t)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:65
    td=L/(i*1.466);
    q=u(i,1);
    s=0
    for j=1:5/delta_t
          
          if t(j)<=td
           p=m.*0.5*((pi/td)^2).*sin((pi./td).*t(j));
           ug(i,j)=0.5*sin((pi./td).*t(j));
          else
              p=0;
              ug(i,j)=0;
          end
           p_hat=p+a1.*u(i,j)+a2.*udot(i,j)+a3.*uddot(i,j);
           u(i,j+1)=p_hat./k_hat;
           udot(i,j+1)=(gama/(beta*delta_t)).*(u(i,j+1)-u(i,j))+(1-(gama/beta))*udot(i,j)+delta_t.*(1-(gama/(2*beta))).*uddot(i,j);
           uddot(i,j+1)=(1/(beta*(delta_t^2))).*(u(i,j+1)-u(i,j))-(1/(beta*delta_t)).*udot(i,j)-((1/(2*beta))-1).*uddot(i,j);
           utotal(i,j)=ug(i,j)+u(i,j);
          if u(i,j+1)>q
             q=u(i,j+1) 
          end
          if utotal(i,j)>s
             s=utotal(i,j);
          end  
    end
          u_r(i)=q ;
          u_total(i)=s;
end
   
    plot(v,u_r*12,'r',v,u_total*12,'k')
    legend('u_r','u_total');grid on;hold on
    xlabel('velocity(mph)','FontSize',12);
    ylabel(' max Displacement(in)','FontSize',12);
    legend('U_r','U-total')




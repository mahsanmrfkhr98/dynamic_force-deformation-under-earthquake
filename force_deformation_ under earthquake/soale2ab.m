clc
clear all
close all
beta=1/6;
gama=1/2;
delta_t=0.02;
Tn=0.001:0.1:15;
Wn=zeros(1,length(Tn));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kisay=0.02;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elcentro=xlsread('elcentro.xlsx');
t=elcentro(:,1);
p=-elcentro(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(length(Tn),length(t));
udot=zeros(length(Tn),length(t));
uddot=zeros(length(Tn),length(t));
p_hat=zeros(length(Tn),length(t));
uddot_total=zeros(length(Tn),length(t));

for i=1:1:length(Tn)
    q=u(1,1);
    b=udot(1,2);
    a=uddot_total(1,2);
    Wn(i)=(2*pi)/Tn(i);
    c(i)=2*kisay*Wn(i);
    if i==1
        Wn=0
        c=0 
    end
    a1=(1/(beta*(delta_t^2)))+(gama/(beta*delta_t))*c(i);
    a2=(1/(beta*delta_t))+((gama/(beta))-1)*c(i);
    a3=((1/(2*beta))-1)*1+delta_t*((gama/(2*beta))-1)*c(i);
    k_hat(i)=(Wn(i)^2)+a1;
      for j=1:1:(length(t)-1)
           P=p';
           F(i,j+1)=P(1,j+1);
           p_hat(1,j+1)=P(1,j+1)+a1.*u(i,j)+a2.*udot(i,j)+a3.*uddot(i,j);
           u(i,j+1)=p_hat(1,j+1)/k_hat(i);
           udot(i,j+1)=(gama/(beta*delta_t)).*(u(i,j+1)-u(i,j))+(1-(gama/beta)).*udot(i,j)+delta_t.*(1-(gama/(2*beta))).*uddot(i,j);
           uddot(i,j+1)=(1/(beta*(delta_t^2))).*(u(i,j+1)-u(i,j))-(1/(beta*delta_t)).*udot(i,j)-(((1/(2*beta))-1).*uddot(i,j));
           uddot_total(i,j+1)=uddot(i,j+1)-F(i,j+1);
            
          if u(i,j+1)>q
              q=u(i,j+1);
          end
                if udot(i,j+1)>b
                   b=udot(i,j+1);
                end
           if uddot_total(i,j+1)>a
             a=uddot_total(i,j+1);
           end
      end
       umax2darsad(i)=q;
       udotmax2darsad(i)=b;
       uddotmax2darsad(i)=a;
       umax2darsad(1)=0;
       udotmax2darsad(1)=0;
       uddotmax2darsad(1)=0;
   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kisay=0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u=zeros(length(Tn),length(t));
udot=zeros(length(Tn),length(t));
uddot=zeros(length(Tn),length(t));
p_hat=zeros(length(Tn),length(t));

for i=1:1:length(Tn)
    q=u(1,1);
    b=udot(1,2);
    a=uddot_total(1,2);
    Wn(i)=(2*pi)/Tn(i);
    c(i)=2*kisay*Wn(i); 
    if i==1
        Wn=0
        c=0
    end
    a1=(1/(beta*(delta_t^2)))+(gama/(beta*delta_t))*c(i);
    a2=(1/(beta*delta_t))+((gama/(beta))-1)*c(i);
    a3=((1/(2*beta))-1)*1+delta_t*((gama/(2*beta))-1)*c(i);
    k_hat(i)=(Wn(i)^2)+a1;
         for j=1:1:(length(t)-1)
           P=p';
           F(i,j+1)=P(1,j+1);
           p_hat(1,j+1)=P(1,j+1)+a1.*u(i,j)+a2.*udot(i,j)+a3.*uddot(i,j);
           u(i,j+1)=p_hat(1,j+1)/k_hat(i);
           udot(i,j+1)=(gama/(beta*delta_t)).*(u(i,j+1)-u(i,j))+(1-(gama/beta)).*udot(i,j)+delta_t.*(1-(gama/(2*beta))).*uddot(i,j);
           uddot(i,j+1)=(1/(beta*(delta_t^2))).*(u(i,j+1)-u(i,j))-(1/(beta*delta_t)).*udot(i,j)-(((1/(2*beta))-1).*uddot(i,j));
           uddot_total(i,j+1)=uddot(i,j+1)-F(i,j+1);
          if u(i,j+1)>q
              q=u(i,j+1);
          end
                if udot(i,j+1)>b
                   b=udot(i,j+1);
                end
           if uddot_total(i,j+1)>a
             a=uddot_total(i,j+1);
           end
         end
      
       umax10darsad(i)=q;
       udotmax10darsad(i)=b;
       uddotmax10darsad(i)=a;
       umax10darsad(1)=0;
       udotmax10darsad(1)=0;
       uddotmax10darsad(1)=0
   
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seismo=xlsread('damped');
T_n=seismo(1:1500,1);
u2darsad=seismo(1:1500,2);
u10darsad=seismo(1:1500,3);

udot2darsad=seismo(1:1500,7);
udot10darsad=seismo(1:1500,8);

uddot2darsad=seismo(1:1500,13);
uddot10darsad=seismo(1:1500,14);




figure(1)
plot(Tn,umax2darsad*981,'b',Tn,umax10darsad*981,'r');grid on;hold on
plot(T_n,u2darsad,"--",T_n,u10darsad,"--")
legend('umax2darsad','umax10darsad','u2darsadseismo','u10darsadseismo')
xlabel('period(sec)','FontSize',12);
ylabel(' max Displacement','FontSize',12);

figure(2)
plot(Tn,udotmax2darsad*981,'b',Tn,udotmax10darsad*981,'r');grid on;hold on
plot(T_n,udot2darsad,"--",T_n,udot10darsad,"--")
legend('udotmax2darsad','udotmax10darsad','udot2darsadseismo','udot10darsadseismo')
xlabel('period(sec)','FontSize',12);
ylabel(' max velocity','FontSize',12);
figure(3)
plot(Tn,uddotmax2darsad,'b',Tn,uddotmax10darsad,'r');grid on;hold on
plot(T_n,uddot2darsad,"--",T_n,uddot10darsad,"--")
legend('uddotmax2darsad','uddotmax10darsad','uddot2darsadseismo','uddot10darsadseismo')
xlabel('period(sec)','FontSize',12);
ylabel(' max Acceleration','FontSize',12);




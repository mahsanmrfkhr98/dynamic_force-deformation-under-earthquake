clc
clear all
close all
beta=1/6;
gama=1/2;
delta_t=0.02;
t_total=2;
r=t_total/delta_t
v=1:1.466/2:95.3;
n=length(v)
V=v/1.466
td=zeros(n,1);
fy=5000;
m=4000/386.22;
k0=800;
kisay=0.4;
wn=sqrt(k0/m)
c=2*m*wn*kisay;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=zeros(r+1,1);

u=zeros(r+1,n);
udot=zeros(r+1,n);
uddot=zeros(r+1,n);
delta_U=zeros(r,n);
delta_Udot=zeros(r,n);
delta_Uddot=zeros(r,n);
fs=zeros(r+1,n);
k=zeros(r+1,n);
k_hat=zeros(r+1,n);
u_r=zeros(1,n);
U_total=zeros(1,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a1=m/(beta*delta_t^2)+(c*gama/(beta*delta_t));
a2=(m/(beta*delta_t))+gama*c/beta;
a3=(gama*c*delta_t/(2*beta))+(0.5*m/beta)-c*delta_t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:1:n
     td=3/v(j);
    for i=1:1:(t_total/delta_t)+1
        t(i)=(i-1)*delta_t;
         if t(i)<td
          p(i,j)=m*6*((pi/td)^2)*sin((pi/td)*t(i));
          ug(i,j)=6*sin((pi/td)*t(i));
         else 
          p(i,j)=0;
         u(i,j)=0;
        end
    end
end

for j=1:1:n
      uddot(1,j)=(p(1,j)-c*udot(1,j)-k0*u(1,j))/m;
     for i=1:1:r
         if abs(fs(i,j))<fy
             k(i,j)=800;
         else
              k(i,j)=0.0001 ;
         end
         k_hat(i,j)=k(i,j)+a1;
         delta_U(i,j)=(p(i+1,j)-p(i,j)+a2*udot(i,j)+a3*uddot(i,j))/k_hat(i,j);
         delta_Udot(i,j)=(1-0.5*gama/beta)*delta_t*uddot(i,j)+(gama/(beta*delta_t))*delta_U(i,j)-gama*udot(i,j)/beta;
         delta_Uddot(i,j)=(1/(beta*(delta_t^2)))*delta_U(i,j)-(1/(beta*delta_t))*udot(i,j)-0.5*uddot(i,j)/beta;
            u(i+1,j)=u(i,j)+delta_U(i,j);
          udot(i+1,j)=delta_Udot(i,j)+udot(i,j);
             fs(i+1,j)=fs(i)+k(i,j)*delta_U(i,j);
              Ut=u+ug;
             if fs(i+1,j)>fy
                  fs(i+1,j)=fy;
            elseif fs(i+1,j)<-fy
               fs(i+1,j)=-fy;
             end
             uddot(i+1,j)=(p(i+1,j)-c*udot(i+1,j)-fs(i+1,j))/m;
     end
    u_r(j)=max(u(:,j));
    U_total(j)=max(Ut(:,j));
end 
plot(V,u_r,'k',V,U_total,'--');grid on; hold on;
xlabel('velocity(mph)','FontSize',12);
ylabel(' MAX Displacement(in)','FontSize',12);
legend('Ur','Utotal');

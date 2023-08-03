% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-2 Disturbance

%% --------------------------------------------- %%
clear all;
close all;
clc;

%%
%%
format long
z = tf('z');

A = 1-0.59*(z^-1)+0.0564*(z^-2);
Acoef = cell2mat(tfdata(A));
a1 = Acoef(2);
a0 = Acoef(3);

B_prime = 1-0.32*(z^-1);
B_prim_coef = cell2mat(tfdata(B_prime));
b1_prime = B_prim_coef(1);
b0_prime = B_prim_coef(2);

n = 2;
d = 3;
n1 = 1;
N = 500;

%square pulse refrence input
sample_number = zeros(N,1);
for i=1:N
    sample_number(i,1)=i;
end
uc_in=zeros(N,1);
for i=1:N/4
    uc_in(i,1)=1;
end
for i=N/4:N/2
    uc_in(i,1)=-1;
end
for i=N/2:3*N/4
    uc_in(i,1)=1;
end
for i=3*N/4:N
    uc_in(i,1)=-1;
end

v=zeros(N,1); % disturbance
for i=1:N
    if i>=40
       v(i,1)=0.1;
    else
        v(i,1)=0;
    end
end
%  a) zero cancelation
syms f1 f2 g0 g1
eq1=f2*a0+g1==0;
eq2=f1*a0+f2*a1+g0==0;
eq3=a0+f1*a1+f2==0;
eq4=a1+f1==0;

SOL=solve([eq1,eq2,eq3,eq4],[f1,f2,g0,g1]);
f1=double(vpa(SOL.f1));
f2=double(vpa(SOL.f2));
g0=double(vpa(SOL.g0));
g1=double(vpa(SOL.g1));
%
G=g0+g1*(z^-1);
F=1+f1*(z^-1)+f2*(z^-2);
beta=F*B_prime;
beta_coef=cell2mat(tfdata(beta));
beta0=beta_coef(1);
beta1=beta_coef(2);
beta2=beta_coef(3);
beta3=beta_coef(4);
%%
beta_prime=z*(beta-beta0);
beta_prime_coef=cell2mat(tfdata(beta_prime))
beta_prime0=beta_prime_coef(1);
beta_prime1=beta_prime_coef(2);
beta_prime2=beta_prime_coef(3);

%%
y_out=zeros(N,1);
u_control=zeros(N,1);
phi_y=zeros(1,2);
phi_u=zeros(1,2);
RR=zeros(1,3);
SS=zeros(1,2);
landa=0.2;
for k=1:N
    for m=1:2
        if k-m<=0
            y_dummy=0; 
        else 
            y_dummy=-y_out(k-m,1); 
        end
        phi_y(1,m)=y_dummy;
       
    end
    
    for ii=3:4
        if k-ii<=0
            u_dummy=0;
        else 
            u_dummy=u_control(k-ii,1)+v(k-ii,1);
        end
        phi_u(1,ii-2)=u_dummy;
    end
    
    y_out(k,1)=phi_y*[a1;a0]+phi_u*[b1_prime;b0_prime];
    
    
    for ii=0:2
        if k-1-ii<=0
            uu=0;
        else
            uu=-u_control(k-1-ii,1);
        end
        RR(1,ii+1)=uu;
    end
    for ii=0:1
        if k-ii<=0
            yy=0;
        else 
            yy=y_out(k-ii,1);
        end
        SS(1,ii+1)=yy;    
    end
    if k-d<=0
        y_star=0;
    else 
        y_star=uc_in(k-d,1);
    end
        
    
    u_control(k,1)=(beta0*(RR*[beta_prime0;beta_prime1;beta_prime2]+y_star-SS*[g0;g1]))/(beta0^2+landa);
   

end
%%
figure 
plot(sample_number,uc_in, "black" ,'LineWidth',1.5)
hold on 
plot(sample_number,y_out,"green" ,'LineWidth',1.5)
ylim([-2 2])
xlabel('Iteration')
legend('Uc','y')
%%
figure 
plot(sample_number,u_control,"black" ,'LineWidth',1.5)
xlabel('Iteration')
legend('System Control Effort')
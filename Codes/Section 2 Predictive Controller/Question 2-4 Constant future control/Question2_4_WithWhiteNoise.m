% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-3 With white noise

%% --------------------------------------------- %%
clear all;
close all;
clc;

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
n=2;
d0=3;
d=6;
n1=1;

%square pulse refrence input
N=300;
sample_number=zeros(N,1);
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

% Generate white noise zero mean
var_noise=0.01;
white_noise=sqrt(var_noise)*randn(N,1);
white_noise=white_noise- mean(white_noise);
%  a) zero cancelation
syms f1 f2 f3 f4 f5 g0 g1
eq1=f5*a0+g1==0;
eq2=f4*a0+f5*a1+g0==0;
eq3=f5+f4*a1+a0*f3==0;
eq4=f4+a1*f3+a0*f2==0;
eq5=f3+a1*f2+a0*f1==0;
eq6=f2+a1*f1+a0==0;
eq7=f1+a1==0;


SOL=solve([eq1,eq2,eq3,eq4,eq5,eq6,eq7],[f1,f2,f3,f4,f5,g0,g1]);
f1=double(vpa(SOL.f1));
f2=double(vpa(SOL.f2));
f3=double(vpa(SOL.f3));
f4=double(vpa(SOL.f4));
f5=double(vpa(SOL.f5));
g0=double(vpa(SOL.g0));
g1=double(vpa(SOL.g1));
%
G=g0+g1*(z^-1);
F=1+f1*(z^-1)+f2*(z^-2)+f3*(z^-3)+f4*(z^-4)+f5*(z^-5);
left=F*B_prime;
left_coef=cell2mat(tfdata(left));

left0=left_coef(1);
left1=left_coef(2);
left2=left_coef(3);
left3=left_coef(4);
left4=left_coef(5);
left5=left_coef(6);
left6=left_coef(7);
%%
R_star_0=left0;
R_star_1=left1;
R_star_2=left2;
R_star_3=left3;
R_star_4=left4;
R_star_5=left5;
Rbar=left6;
%%
sum_R_star=R_star_0+R_star_1+R_star_2+R_star_3+R_star_4+R_star_5;
y_out=zeros(N,1);
u_control=zeros(N,1);
phi_y=zeros(1,2);
phi_u=zeros(1,2);
RR=zeros(1,1);
SS=zeros(1,2);
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
            u_dummy=u_control(k-ii,1);
        end
        phi_u(1,ii-2)=u_dummy;
    end
    
    y_out(k,1)=phi_y*[a1;a0]+phi_u*[b1_prime;b0_prime]+white_noise(k,1);
    
    
    for ii=1:1
        if k-ii<=0
            uu=0;
        else
            uu=-u_control(k-ii,1);
        end
        RR(1,ii)=uu;
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
        
    
    u_control(k,1)=(Rbar*RR+y_star-SS*[g0;g1])/sum_R_star;
   

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

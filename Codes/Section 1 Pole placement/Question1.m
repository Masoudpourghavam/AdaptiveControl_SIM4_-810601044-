% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 1

%% --------------------------------------------- %%
clear all;
close all;
clc;

%%

format long

z = tf('z');
A = (z^3)*(z-0.12)*(z-0.47);
B = (z-0.32);

Acoef = cell2mat(tfdata(A));
a4 = Acoef(2);
a3 = Acoef(3);

Bcoef = cell2mat(tfdata(B));
b1 = Bcoef(1);
b0 = Bcoef(2);

degA = 5;    degB = 1;

%% desired pole location 
Mp = 0.15;
zeta = ((log(Mp)^2)/(pi^2+log(Mp)^2))^0.5;
ts = 2;
sigma = 4/ts;
wn = sigma/zeta;
s1 = -zeta*wn+i*(wn*(1-zeta^2)^0.5);
s2 = -zeta*wn-i*(wn*(1-zeta^2)^0.5);
s3 = -15*zeta*wn;
s4 = -15*zeta*wn;
s5 = -15*zeta*wn;
Ts = 1;

Pz1 = exp(s1*Ts);
Pz2 = exp(s2*Ts);
Pz3 = exp(s3*Ts);
Pz4 = exp(s4*Ts);
Pz5 = exp(s5*Ts);

Am = (z-Pz1)*(z-Pz2)*(z-Pz3)*(z-Pz4)*(z-Pz5);
Amcoef = cell2mat(tfdata(Am));

am4 = Amcoef(2);
am3 = Amcoef(3);
am2 = Amcoef(4);
am1 = Amcoef(5);
am0 = Amcoef(6);

%% square pulse refrence input
N = 300;
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

%%  a) zero cancelation
degB_plus=degB;
degB_minus=degB-degB_plus;
B_minus=Bcoef(1);
B_plus=B/B_minus;
degBm=degB;
bm0=evalfr(Am,1);
Bm=bm0*z^(degBm);
deg_Bprim_m=degBm-degB_minus;
Bprim_m=(bm0/B_minus)*z^(deg_Bprim_m);
deg_A0=degA-degB_plus-1;
A0=z^(deg_A0);
T=A0*Bprim_m;
Tcoef=cell2mat(tfdata(T));
t3=Tcoef(1);

%%
syms s4 s3 s2 s1 s0 r_prim1 r_prim0 r_prim2
eq0=r_prim2+a4-am4==0;
eq1=r_prim1+a4*r_prim2+a3-am3==0;
eq2=r_prim0+a4*r_prim1+a3*r_prim2-am2==0;
eq3=a4*r_prim0+a3*r_prim1+s4-am1==0;
eq4=a3*r_prim0+s3-am0==0;
eq5=s2==0;
eq6=s1==0;
eq7=s0==0;
SOL=solve([eq0,eq1,eq2,eq3,eq4,eq5,eq6],[s4,s3,s2,s1,s0,r_prim2,r_prim1,r_prim0]);
s4=double(vpa(SOL.s4));
s3=double(vpa(SOL.s3));
s2=double(vpa(SOL.s2));
s1=double(vpa(SOL.s1));
s0=double(vpa(SOL.s0));
r_prim2=double(vpa(SOL.r_prim2));
r_prim1=double(vpa(SOL.r_prim1));
r_prim0=double(vpa(SOL.r_prim0));
%
R_prim=z^3+r_prim2*z^2+r_prim1*z+r_prim0;
R=R_prim*B_plus;
Rcoef=cell2mat(tfdata(R));
r3=Rcoef(2);
r2=Rcoef(3);
r1=Rcoef(4);
r0=Rcoef(5);
S=s4*z^4+s3*(z^3)+s2*(z^2)+s1*z+s0;
Scoef=cell2mat(tfdata(S));
%%
y_out=zeros(N,1);
u_control=zeros(N,1);
phi_y=zeros(1,2);
phi_u=zeros(1,2);
RR=zeros(1,4);
SS=zeros(1,5);
for k=1:N
    for m=1:2
        if k-m<=0
            y_dummy=0; 
        else 
            y_dummy=-y_out(k-m,1); 
        end
        phi_y(1,m)=y_dummy;
    end
    
    for ii=4:5
        if k-ii<=0
            u_dummy=0;
        else 
            u_dummy=u_control(k-ii,1);
        end
        phi_u(1,ii-3)=u_dummy;
    end
    
    y_out(k,1)=phi_y*[a4;a3]+phi_u*[b1;b0];
    
    for ii=1:3
        if k-ii<=0
            uu=0;
        else
            uu=-u_control(k-ii,1);
        end
        RR(1,ii)=uu;
    end
    for ii=0:3
        if k-ii<=0
            yy=0;
        else 
            yy=y_out(k-ii,1);
        end
        SS(1,ii+1)=yy;    
    end
    
    u_control(k,1)=RR*[r3;r2;r1;r0]+t3*uc_in(k,1)-SS*[s4;s3;s2;s1;s0];
end
%%
figure 
plot(sample_number,uc_in, "green",'LineWidth',2)
hold on 
plot(sample_number,1.28*y_out, "black",'LineWidth',1.5)
ylim([-2 2])
xlabel('Iteration')
legend('Uc','y')
%%
figure 
plot(sample_number,1.28*u_control, "black",'LineWidth',1.5)
xlabel('Iter')
legend('System Control Effort')


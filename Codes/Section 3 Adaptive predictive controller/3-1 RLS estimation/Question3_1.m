% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 3-1

%% --------------------------------------------- %%
clear all;
close all;
clc;

%% Parameter extraction
n=6;
a1=0.59;
a2=-0.0564;
b1=0;
b2=0;
b3=1;
b4=-0.32;
theta=[a1;a2;b1;b2;b3;b4];
N = 200;
var_noise=0.056;
white_noise=sqrt(var_noise)*randn(N,1);
white_noise=white_noise- mean(white_noise);

%% Input Specification
%%white noise input
u=sqrt(var_noise)*randn(N,1);

%%pulse input 
% u=zeros(N,1);
% u(1,1)=1;

%%step input
% u=ones(N,1);

%%sinpsoidal input
% u=zeros(N,1);
% omega=0.01;
% for i=1:N
%    u(i,1)=sin(omega*i*Ts);
% end

%%ramp input
% u=zeros(N,1); 
% for i=1:N
%    u(i,1)=i/(N*Ts);
% end
%% calculate output without noise 
y_noiseless=zeros(length(u),1);
y_noiseless(1,1)=0;
y_noiseless(2,1)=[-y_noiseless(1,1) 0 u(1,1) 0 0 0]*theta;
y_noiseless(3,1)=[-y_noiseless(2,1) -y_noiseless(1,1) u(2,1) u(1,1) 0 0]*theta;
y_noiseless(4,1)=[-y_noiseless(3,1) -y_noiseless(2,1) u(3,1) u(2,1) u(1,1) 0]*theta;
for i=5:length(y_noiseless)
    y_noiseless(i,1)=[-y_noiseless(i-1,1) -y_noiseless(i-2,1) u(i-1,1) u(i-2,1) u(i-3,1) u(i-4,1)]*theta;
end
%initial condition
theta_hat=zeros(n,1);
p=1000*eye(n);
phi_t=zeros(1,n);
epsilon=10e-12;

a_hat1=zeros(N,1);
a_hat2=zeros(N,1);
b_hat1=zeros(N,1);
b_hat2=zeros(N,1);
b_hat3=zeros(N,1);
b_hat4=zeros(N,1);

sample_number=zeros(N,1);
for i=1:N
    sample_number(i,1)=i;
end
%%
for k=2:N
    for m=1:4
        if k-m<=0
            u0=0;
        else 
            u0=u(k-m,1);
        end
        phi_t(1,m+2)=u0;
    end
    for f=1:2
        if k-f<=0
            y0=0;
        else 
            y0=-y_noiseless(k-f,1);
        end
        phi_t(1,f)=y0;
    end


    p=p-((p*(phi_t')*phi_t*p)/(1+phi_t*p*(phi_t')));
    gain=p*(phi_t');
    theta_hat=theta_hat+gain*(y_noiseless(k,1)-(phi_t*theta_hat));
    a_hat1(k,1)=theta_hat(1);
    a_hat2(k,1)=theta_hat(2);
    b_hat1(k,1)=theta_hat(3);
    b_hat2(k,1)=theta_hat(4);
    b_hat3(k,1)=theta_hat(5);
    b_hat4(k,1)=theta_hat(6);
    I=(abs(y_noiseless(k,1)-(phi_t*theta_hat)));
end
%%
a1_act=zeros(N,1);
a2_act=zeros(N,1);
b1_act=zeros(N,1);
b2_act=zeros(N,1);
b3_act=zeros(N,1);
b4_act=zeros(N,1);
a1_act(:,1)=a1;
a2_act(:,1)=a2;
b1_act(:,1)=b1;
b2_act(:,1)=b2;
b3_act(:,1)=b3;
b4_act(:,1)=b4;

%% Plot, results
subplot(2,2,1)
plot(sample_number,a1_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,a_hat1,'g' ,'LineWidth',1)
legend('a1-actual','a1-hat')
xlim([0 N])
ylim([-3 3])
xlabel('sample number')

subplot(2,2,2)
plot(sample_number,a2_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,a_hat2,'g' ,'LineWidth',1)
legend('a2-actual','a2-hat')
xlim([0 N])
ylim([-3 3])
xlabel('sample number')


figure
subplot(2,2,1)
plot(sample_number,b1_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,b_hat1,'g' ,'LineWidth',1)
legend('b1-actual','b1-hat')
xlim([0 N])
ylim([-1 1])
xlabel('sample number')
subplot(2,2,2)
plot(sample_number,b2_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,b_hat2,'g' ,'LineWidth',1)
legend('b2-actual','b2-hat')
xlim([0 N])
ylim([-1 1])
xlabel('sample number')
subplot(2,2,3)
plot(sample_number,b3_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,b_hat3,'g' ,'LineWidth',1)
legend('b3-actual','b3-hat')
xlim([0 N])
ylim([-1 1])
xlabel('sample number')
subplot(2,2,4)
plot(sample_number,b4_act,'--k' ,'LineWidth',1.5)
hold on
plot(sample_number,b_hat4,'g' ,'LineWidth',1)
legend('b4-actual','b4-hat')
xlim([0 N])
ylim([-1 1])
xlabel('sample number')
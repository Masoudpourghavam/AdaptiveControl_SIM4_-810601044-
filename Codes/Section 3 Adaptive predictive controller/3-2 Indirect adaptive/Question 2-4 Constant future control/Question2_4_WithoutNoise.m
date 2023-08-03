% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 3-2 Indirect adaptive const. future control without noise

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

%square pulse reference input
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

% RLS estimation parameters
M = 2;            % Number of filter coefficients (order)
lambda = 0.99;    % Forgetting factor
delta = 0.1;      % Small positive constant

% Initialize RLS variables
P = delta^(-1) * eye(M);  % Inverse correlation matrix
theta_hat = zeros(M, 1); % Filter coefficient estimate

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
    
    y_out(k,1)=phi_y*[a1;a0]+phi_u*[b1_prime;b0_prime];
    
    
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
        
    % RLS estimation
    phi = [RR, y_star - SS*theta_hat];
    e = y_out(k) - phi*theta_hat;
    K = P*phi' / (lambda + phi*P*phi');
    theta_hat = theta_hat + K*e;
    P = (P - K*phi*P) / lambda;
    
    % Adaptive control calculation
    g0 = theta_hat(1);
    g1 = theta_hat(2);
    u_control(k) = uc_in(k) - phi*[g0; g1];
end

%%
figure 
plot(sample_number, uc_in, 'black' ,'LineWidth',1.5)
hold on 
plot(sample_number, y_out, 'green' ,'LineWidth',1.5)
ylim([-2 2])
xlabel('Iteration')
legend('Uc', 'y')

%%

figure 
plot(sample_number, u_control, 'black' ,'LineWidth',1.5)
xlabel('Iteration')
legend('System Control Effort')

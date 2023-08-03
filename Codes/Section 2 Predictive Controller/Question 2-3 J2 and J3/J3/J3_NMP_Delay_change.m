% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 2-3 J3 Delay change

%% --------------------------------------------- %%
clear all;
close all;
clc;

%%

%%% One Step Ahead Control, with J3 Cost Functionn, NMP system %%%


%% Initialization 
Ts = 0.035;
tMax = 20;
t = 0:Ts:tMax;

q = tf('q',Ts);
d0 = 4;      % Delay
d = 2;

A = (q-0.12)*(q-0.47);
A = tfdata(A);
A = A{1};


B = (q-0.32);
B = tfdata(B);
B = B{1};

n = numel(A) - 1;
m = numel(B) - 1;
ap=A;
bp=1;
dp=[1 zeros(1,n+d-1)];

%% Controller Parameters
[ G, F ] = Diophantine_solver(A,d,n);

Bprime = B;
alpha = G;

beta = conv(F,Bprime);  %beta(q)

%% bETAPrime-General form of all cases of delay

beta0 = beta(1);
betamin=0.01;       % saturation for beta0
if beta0<betamin
    beta0=betamin;
end
Q = [1 0];     % q
betaPrime = beta(2:end);

% define the Suitable values of lambda for different values of delay-which can be altered when the
% open loop system changes.
if d == 1
    lambda = 3.5;   % Critical value, in order for the closed loop system to be stable
    K = lambda/beta0;
    Ac = plus([Bprime 0],K*A);                % Backward Operator

elseif d == 2
    lambda = -.005;                      
    K = lambda/beta0;
    Ac = plus([Bprime 0],K*A);                
    
elseif d == 3
    lambda = 10;                             
    K = lambda/beta0;
    Ac = plus([Bprime 0],K*A);              
    
else
    lambda = 2.5;                         
    K = lambda/beta0;
    Ac = plus([Bprime 0],K*A);           
    
end
lambda
%% The Roots of the Inverse Model

roots(Ac)
f1 = figure(1);
Ac = Ac(1)*q^2 + Ac(2)*q + Ac(end);
G_z = 1/Ac;
rlocus(G_z);
title(' The closed loop poles locus, for the NMP controller system')
grid on

%% Pulse Reference

Tau = 1;
Tf = tMax;
[u0,~] = gensig('square',Tau,Tf);
Uc  = 2*u0 - 1;

%% Cosine Reference

% Amp = 0.1;
% w = 0.25;
% Uc = Amp*cos(pi*w*t);

%% Ramp Reference

% k = 1.5;  % Slop
% Uc = k*t;

%% White Noise

Var = 0.05;
e = Var*randn(size(t));

%% Disturbance

V = 0.5;                                  % Disurbance Amplitude
ApplyTime = round(length(t)/2);           % Arbitrarly
Disturbance = [zeros(1,ApplyTime) V*ones(1,numel(t) - ApplyTime+1)];

%% Upper and Lower Bounds of the control inputs, assumed to be given.

Lb_u = -.4;
Ub_u = .4;

%% J3 parameters
P = ([1 -.7 -.2]);
R = conv([1 -.8 ],poly(1));
figure
rlocus(tf(1,conv(P,[B,0])+lambda/beta0.*conv(A,R),.1))

u_hat=zeros(size(t));
R_prime = R(2:end);
P_prime = P(2:end);
%% Simulation

u = zeros(size(t));
y = u;
CONSTANT = beta0/(beta0^2 + lambda);        % Constant coefficient of the control signal


for i=d+n:length(t)-d
    
    Predicted_y = Uc(i+d);  % In order for tracking
%     y(i+d) = (-[y(i+d-1) y(i+d-2)]*Ac(2:end)' +...
%         [Uc(i-1) Uc(i-2)]*betaPrime(1:end-1)')/Ac(1); % + e(i)

    y(i)=-A(2:end)*y(i-1:-1:i-n)'+B*u(i-d:-1:i-length(B)+1-d)';        
    if d == 1
          u(i) = CONSTANT*(-u(i-2)*betaPrime + Predicted_y -...
              [y(i) y(i-1)]*alpha'+lambda*P_prime*[u_hat(i-1) u_hat(i-2)]'-lambda*R_prime*[u(i-1) u(i-2)]'); % + Disturbance(i)
          u_hat(i) = -P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]'; 
            % saturation
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
        
    elseif d == 2

        u(i) = CONSTANT*(-[u(i-1) u(i-2)]*betaPrime' + Predicted_y - [y(i) y(i-1)]*alpha'+lambda*P_prime*[u_hat(i-1) u_hat(i-2)]'-lambda*R_prime*[u(i-1) u(i-2)]'); % + Disturbance(i)
        u_hat(i) = -P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]'; 
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
    elseif d == 3
            
        u(i) = CONSTANT*(-[u(i-1) u(i-2) u(i-3)]*betaPrime'...
            + Predicted_y - [y(i) y(i-1)]*alpha'-lambda*(-P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]')); % + Disturbance(i)
        u_hat(i) = -P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]'; 
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
    else
        
        u(i) = CONSTANT*(-[u(i-1) u(i-2) u(i-3) u(i-4)]*betaPrime'...
            + Predicted_y - [y(i) y(i-1)]*alpha'-lambda*(-P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]')); % + Disturbance(i)
        u_hat(i) = -P_prime*[u_hat(i-1) u_hat(i-2)]'+R_prime*[u(i-1) u(i-2)]'; 
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
    end
    
end

%% Plot Results

f2 = figure(2);
subplot(2,1,1)
plot(y,'g','LineWidth',2)
hold on
plot(Uc,'k--','LineWidth',2.15)
xlabel('Sample','Interpreter','Latex','FontWeight','Bold')
title('Output Vs Reference Pulse','Interpreter','Latex','FontWeight','Bold')
grid on
xlim([0 length(t) - 1]);

y(1:d+n) = [];              % Corresponds to the initial conditions
Index = find(y == 0);
y(Index) = Uc(Index);

error = y - Uc(d+n-1:numel(y)-(d+n-1));
MSE = mse(error)

subplot(2,1,2)
plot(u,'k','LineWidth',2)
xlabel('Sample','Interpreter','Latex','FontWeight','Bold')
title('Control Signal','Interpreter','Latex','Fontweight','Bold')
grid on
xlim([0 length(t) - 1]);

movegui(f1,'east')
movegui(f2,'west')

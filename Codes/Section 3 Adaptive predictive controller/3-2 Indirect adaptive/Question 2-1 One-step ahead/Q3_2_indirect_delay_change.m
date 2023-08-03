% Adaptive Control - Simulation 4
% Masoud Pourghavam
% Student Number: 810601044
% Question 3-2 delay change

%% --------------------------------------------- %%
clear all;
close all;
clc;


%% Initialization 

Ts = 0.035;
tMax = 20;
t=0:Ts:tMax;
q = tf('q',Ts);
d0 = 4;      % Delay
d = 4;

A = (q-0.12)*(q-0.47);
A = tfdata(A);
B = (q-0.32);

A = A{1};
B = tfdata(B);
B = B{1};

n = numel(A) - 1;
m = numel(B) - 1;
ap=A;
bp=1;
dp=[1 zeros(1,n+d-1)];

Q = [1 0];     % q
% betaPrime = beta(2:end);

if d == 1
    landa = 3.5;  
%     K = lambda/beta0;
%     Ac = plus([Bprime 0],K*A);                % Backward Operator

elseif d == 2
    landa = 0.8;                      
%     K = lambda/beta0;
%     Ac = plus([Bprime 0],K*A);                
    
elseif d == 3
    landa = 10;                             
%     K = lambda/beta0;
%     Ac = plus([Bprime 0],K*A);              
    
else
    landa = 2.5;                         
%     K = lambda/beta0;
%     Ac = plus([Bprime 0],K*A);           
%     
end

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

% Var = 0.05;
% e = Var*randn(size(t));

%% Disturbance

% V = 0.5;                                  % Disurbance Amplitude
% ApplyTime = round(length(t)/2);           % Arbitrarly
% Disturbance = [zeros(1,ApplyTime) V*ones(1,numel(t) - ApplyTime+1)];

%% Upper and Lower Bounds of the control inputs, assumed to be given.

Lb_u = -.4;
Ub_u = .4;

%% Simulation

u = zeros(size(t));
y = u;
% CONSTANT = beta0/(beta0^2 + lambda);        % Constant coefficient of the control signal

% PA Algorithm parameters
gamma = 0.95;
aalpha = 0.01;
nParams = n+(m+1);                 % Number of Parameters
A_Hat = [1 [-0.1 0]];
B_Hat = [0.5 0.25];
Theta_Hat(:,d+n-1)=[A_Hat(2:end) B_Hat];  

for i=d+n:length(t)-d
    
    %% PA Algorithm
    
    y(i) = -y(i-1:-1:i-n)*A(2:end)' + [u(i-d0) u(i-1-d0)]*B'; % + e(i)
    phi=[-y(i-1:-1:i-n),u(i-d:-1:i-d-n+1)]';
    pa = (gamma*phi)/(aalpha+phi'*phi);
    Theta_Hat(:,i)=Theta_Hat(:,i-1)+pa*(y(i)-phi'*Theta_Hat(:,i-1));
    
    THETAHAT = Theta_Hat(:,i);
    A_Hat = [1 THETAHAT(1:n)'];
    B_Hat = THETAHAT(n+1:end);
    
    %% Controller Parameters Estimation

    [ G, F ] = Diophantine_solver( A , d );
    Bprime = B_Hat';
    alpha = G;
    beta = conv(F,Bprime);
    beta0 = beta(1);
    betamin=0.01;       % saturation for beta0
    if beta0<betamin
        beta0=betamin;
    end
    betaPrime=beta(2:end);
    CONSTANT = 1;        % Unweighted controller

   
    Predicted_y = Uc(i+d);  % In order for tracking
%     y(i+d) = (-[y(i+d-1) y(i+d-2)]*Ac(2:end)' +...
%         [Uc(i-1) Uc(i-2)]*betaPrime(1:end-1)')/Ac(1); % + e(i)

    y(i) = -[y(i-1) y(i-2)]*A(2:end)' + [u(i-d) u(i-1-d)]*Bprime'; % + e(i)
        
    if d == 1
        u(i) = CONSTANT*(-u(i-2)*betaPrime + Predicted_y -...
            [y(i) y(i-1)]*alpha'); % + Disturbance(i)
            % saturation
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
        
    elseif d == 2

        u(1,i) = CONSTANT*(-[u(i-1) u(i-2)]*betaPrime' + Predicted_y - [y(i) y(i-1)]*alpha'); % + Disturbance(i)
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
    elseif d == 3
            
        u(i) = CONSTANT*(-[u(i-1) u(i-2) u(i-3)]*betaPrime'...
            + Predicted_y - [y(i) y(i-1)]*alpha'); % + Disturbance(i)
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end
    else
        
        u(i) = CONSTANT*(-[u(i-1) u(i-2) u(i-3) u(i-4)]*betaPrime'...
            + Predicted_y - [y(i) y(i-1)]*alpha'); % + Disturbance(i)
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
landa

subplot(2,1,2)
plot(u,'k','LineWidth',2)
xlabel('Sample','Interpreter','Latex','FontWeight','Bold')
title('Control Signal','Interpreter','Latex','Fontweight','Bold')
grid on
xlim([0 length(t) - 1]);

% movegui(f1,'east')
% movegui(f2,'west')

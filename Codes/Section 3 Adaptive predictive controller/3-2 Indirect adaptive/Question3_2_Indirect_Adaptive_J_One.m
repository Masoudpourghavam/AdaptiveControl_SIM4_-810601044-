clc;clear;close all;

%{ 
     %%% One Step Ahead, J1, Indirect Adaptive  %%%
%}

%% Initialization 

Ts = 0.025;
tMax = 15;
t = 0:Ts:tMax;

q = tf('q',Ts);
d = 3;   % Delay

% System
A = [1 -0.59 0.0564];
B = [1 -0.32];

n = numel(A) - 1;
m = numel(B) - 1;

%% Pulse Reference

Tau = 1;
Tf = t(end);
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

Var = 0.025;
e = Var*randn(size(t));

%% Disturbance

V = 0.5;        % Disurbance Amplitude
ApplyTime = round(length(t)/2); % Arbitrarly
Disturbance = [zeros(1,ApplyTime) V*ones(1,numel(t) - ApplyTime+1)];

%% Upper and Lower Bounds of the control inputs, assumed to be given.

Lb_u = -1.5;
Ub_u = 1.5;

%% Simulation

u = zeros(size(t));
y = u;
% PA Algorithm parameters
gamma = 0.95;
aalpha = 0.01;
nParams = n+(m+1);                 % Number of Parameters
A_Hat = [1 [-0.1 0]];
B_Hat = [0.5 0.25];
Theta_Hat(:,d+n-1)=[A_Hat(2:end) B_Hat];  

%% Different system delays in order to observe the outcomes

sysD = 2;   % Main System Delay
% sysD = 1;
% sysD = 3;

SuddenIter = 400;   % Sudden Iteration of change of parameters

for i=d+n:length(t)-d
    %% The Effect of Sudden Change of System Parameters
            
%     if i >=SuddenIter
%         
%         A = poly([-0.4 0.1]);
%         B = [-0.1 -0.05];
%         
%     end

    %% PA Algorithm
    
    y(i) = -y(i-1:-1:i-n)*A(2:end)' + [u(i-sysD) u(i-1-sysD)]*B'; % + e(i)
    phi=[-y(i-1:-1:i-n),u(i-d:-1:i-d-n+1)]';
    pa = (gamma*phi)/(aalpha+phi'*phi);
    Theta_Hat(:,i)=Theta_Hat(:,i-1)+pa*(y(i)-phi'*Theta_Hat(:,i-1));
    
    THETAHAT = Theta_Hat(:,i);
    A_Hat = [1 THETAHAT(1:n)'];
    B_Hat = THETAHAT(n+1:end);
    
    %% Controller Parameters Estimation

    [F, G] = OneStepAhead_Dioph(A_Hat,d,n);
    Bprime = B_Hat';
    alpha = G;
    beta = conv(F,Bprime);
    beta0 = beta(1)
    
    Predicted_y = Uc(i+d);  % In order for tracking
        
        if d == 1
            
            u(i) = (-u(i-1)*beta(end) + Predicted_y -...
                [y(i) y(i-1)]*alpha')/beta0; % + Disturbance(i)
            
            % Saturation
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end

        elseif d == 2

            u(i) = (-[u(i-1) u(i-2)]*beta(2:end)' +...
                Predicted_y - [y(i) y(i-1)]*alpha')/beta0; % + Disturbance(i)
            
            % Saturation
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end

        elseif d == 3

            u(i) = (-[u(i-1) u(i-2) u(i-3)]*beta(2:end)'...
                + Predicted_y - [y(i) y(i-1)]*alpha')/beta0; % + Disturbance(i)
            
            % Saturation
            if u(i) <= Ub_u && u(i) >= Lb_u
                u(i) = u(i);
            
            elseif u(i) < Lb_u
                u(i) = Lb_u;
                
            else
                u(i) = Ub_u;
                
            end

        else

            u(i) = (-[u(i-1) u(i-2) u(i-3) u(i-4)]*beta(2:end)'...
                + Predicted_y - [y(i) y(i-1)]*alpha')/beta0; % + Disturbance(i)
            
            % Saturation
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

f1 = figure(1);
plot(y,'b','LineWidth',2)
hold on
plot(Uc,'k--','LineWidth',2.15)
xlabel('Sample','Interpreter','Latex','FontWeight','Bold')
title('Output Vs Reference Pulse','Interpreter','Latex','FontWeight','Bold')
grid on
xlim([0 length(t) - 1]);

f2 = figure(2);
plot(u,'m','LineWidth',2)
xlabel('Sample','Interpreter','Latex','FontWeight','Bold')
title('Control Signal','Interpreter','Latex','Fontweight','Bold')
grid on
xlim([0 length(t) - 1]);

f3 = figure(3);
A_Main = [1 -0.8-0.25 0.25*0.8];
B_Main = [1 -0.75];

Theta_Main = [A_Main(2:end)' 
                B_Main'];
Name = {'a_{1}','a_{2}','b_{0}','b_{1}'};

for i=1:size(Theta_Hat,1)
   
    subplot(size(Theta_Hat,1),1,i)
    plot(Theta_Hat(i,:),'r','LineWidth',2)
    hold on
    grid on
    plot(Theta_Main(i)*ones(size(t)),'k--','LineWidth',2)
    xlabel('Sample','FontWeight','Bold')
    ylabel(Name{i},'FontWeight','Bold')
    legend('\theta_{hat}','\theta')
    if i>1
        legend off
    end
end

movegui(f1,'east')
movegui(f2,'center')
movegui(f3,'west')


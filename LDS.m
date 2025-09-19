clc;
clear all;
close all;


function dydt = forcesdof(t, y, m, c, k, F0, omega_f)
    F_t = F0 * sin(omega_f * t);   % Excitation force
    dydt = zeros(2,1);
    dydt(1) = y(2);
    dydt(2) = (F_t - c*y(2) - k*y(1)) / m;
end


% system parameter
m = 10;    % kg
k = 10000; % Stiffness
c = 20;    % Damping

% not required
omega_n = sqrt(k/m); 
zeta    = c / (2 * m * omega_n);

dt     = 0.001;
t_span = 0:dt:10;

% Initial condition
x0 = 0.001;
v0 = 0;
y0 = [x0; v0];

F0      = 3;
omega_f = 15; % Excitation Freq

[t,y] = ode45(@(t,y) forcesdof(t, y, m, c, k, F0, omega_f), t_span, y0);

figure;
plot(t,y(:,1),'b','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Displacement (m)');
title('Forced Vibration Response');
grid on;

figure;
plot(t,y(:,2),'r','LineWidth',1.5);
xlabel('Time (s)');
ylabel('Velocity (m/s)');
title('Velocity Response');
grid on;


figure;
plot(y(:,1),y(:,2),'r','LineWidth',1.5);
xlabel('Velocity (m/s)');
ylabel('Displacement (m)');
title('State Space');
grid on;
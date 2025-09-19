clc
clear
close all

%% Function to calculate y1 dot  |  y2 dot  |   y3 dot  |  y4 dot
function dydt = forcesdof(t, y, m1, m2, l1, l2, g)
    dydt = zeros(4,1);
    dydt(1) = y(2);         % y1_dot = y2
    dydt(3) = y(4);         % y3_dot = y4
    
    Delta = (m1 + m2)*l1 - m2*l1*(cos(y(1)-y(3)))^2;
    
    dydt(2) = (1/Delta) * ( ...
        -m2*l2*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)) + ...
        m2*g*sin(y(3))*cos(y(1)-y(3)) - ...
        m2*l2*y(4)^2*sin(y(1)-y(3)) - ...
        (m1 + m2)*g*sin(y(1)) ...
    );
    
    denom = l2 - (m2*l2/(m1 + m2))*(cos(y(1)-y(3)))^2;
    
    dydt(4) = (1/denom) * ( ...
        l1*y(2)^2*sin(y(1)-y(3)) - g*sin(y(3)) + ...
        (m2*l2/(m1 + m2))*y(4)^2*sin(y(1)-y(3))*cos(y(1)-y(3)) + ...
        g*sin(y(1))*cos(y(1)-y(3)) ...
    );
end


%% Parameters
dt = 0.001;
m1 = 10.0;   % Mass of first pendulum (kg)
m2 = 20.0;   % Mass of second pendulum (kg)
l1 = 90.0;   % Length of first pendulum (m)
l2 = 90.0;   % Length of second pendulum (m)
g = 9.81;   % Acceleration due to gravity (m/s^2)

% %Initial conditions
y_1 = 1; % Initial angle of first pendulum (rad)
y_2 = 0; % Initial angular velocity of first pendulum (rad/s)
y_3 = 1; % Initial angle of second pendulum (rad)
y_4 = 0; % Initial angular velocity of second pendulum (rad/s)

y0 = [y_1; y_2; y_3; y_4];
t_span = 0:dt:100;

%% Pass all parameters to the function
[t,y] = ode45(@(t,y) forcesdof(t, y, m1, m2, l1, l2, g), t_span, y0);


%% PLotting pendulum path
% Extract angles from solution arrays
theta1 = y(:,1);
theta2 = y(:,3);

% Compute XY coordinates of the pendulum bobs
x1 = l1 * sin(theta1);
y1 = -l1 * cos(theta1);

x2 = x1 + l2 * sin(theta2);
y2 = y1 - l2 * cos(theta2);

% Plot trajectory of pendulum bobs
figure;
plot(x1, y1, 'b', 'LineWidth', 2);      % Plot path of first bob
hold on;
plot(x2, y2, 'r', 'LineWidth', 2);      % Plot path of second bob

xlabel('X position (m)');
ylabel('Y position (m)');
title('Double Pendulum Trajectory');
legend('Bob 1', 'Bob 2');
axis equal;
grid on;


%% PLotting state space trajectories
figure;

subplot(2,2,1);
plot(t, y(:,1));
xlabel('Time (s)');
ylabel('\theta_1 (rad)');
title('Trajectory of theta_1');
grid on;

subplot(2,2,2);
plot(t, y(:,2));
xlabel('Time (s)');
ylabel('\dot{\theta}_1 (rad/s)');
title('Trajectory of theta_1 dot');
grid on;

subplot(2,2,3);
plot(t, y(:,3));
xlabel('Time (s)');
ylabel('\theta_2 (rad)');
title('Trajectory of theta_2');
grid on;

subplot(2,2,4);
plot(t, y(:,4));
xlabel('Time (s)');
ylabel('\dot{\theta}_2 (rad/s)');
title('Trajectory of theta_2 dot');
grid on;

function dy = double_pendulum_ode(t, y, m1, m2, l1, l2, g)
    theta1 = y(1);
    dtheta1 = y(2);
    theta2 = y(3);
    dtheta2 = y(4);

    % Precompute useful quantities
    delta = theta1 - theta2;
    A = (m1 + m2) * l1^2;
    B = m2 * l1 * l2 * cos(delta);
    C = l2;
    D = l1 * cos(delta);
    E = m2 * l1 * l2 * dtheta2^2 * sin(delta) + (m1 + m2) * g * l1 * sin(theta1);
    F = l1 * dtheta1^2 * sin(delta) + g * sin(theta2); 

    % Matrix for equations
    denominator = A*C - B*D;

    ddtheta1 = (C*(-E) - B*(-F)) / denominator;
    ddtheta2 = (A*(-F) - D*(-E)) / denominator;

    dy = zeros(4,1);
    dy(1) = dtheta1;
    dy(2) = ddtheta1;
    dy(3) = dtheta2;
    dy(4) = ddtheta2;
end

% Parameters
m1 = 1.0;   % Mass of first pendulum (kg)
m2 = 1.0;   % Mass of second pendulum (kg)
l1 = 1.0;   % Length of first pendulum (m)
l2 = 1.0;   % Length of second pendulum (m)
g = 9.81;   % Acceleration due to gravity (m/s^2)

% Time span for simulation
tspan = [0 10];

% Initial conditions: [theta1, dtheta1, theta2, dtheta2]
y0 = [pi/2; 0; pi/2; 0];


[t, y] = ode45(@(t, y) double_pendulum_ode(t, y, m1, m2, l1, l2, g), tspan, y0);
length(y)

clear; 
clc; 

%% Constant Definitions
M = 5; % Mass of the rod (kg)
m = 10; % Mass of the cart (kg)
L = 2; % Length of the rod (m)
I = 3; % Inertia of the rod (kg*m^2)
g = 9.81; % Gravitational acceleration (m/s^2)
operate = [0 ; 0 ; pi/2 ; 0 ; 0]; % The poing of operation
cartW = 2; % Cart Width (m)
cartH = 1; % Cart Height (m)
wheelR = 0.25; % Wheel radius (m)
Ts = 0.01; % Time step (s)
finalTime = 50; % Final time (s)

%% Variable Definitions
syms x xD xDD a aD aDD b bD bDD c cD cDD u; 

%% Equation Solutions
ax1 = xDD - L*sin(a)*aDD - (L*sin(a + b + c)*aDD)/2 - (L*sin(a + b + c)*bDD)/2 - ...
    (L*sin(a + b + c)*cDD)/2 - L*cos(a + b)*aD^2 - L*cos(a + b)*bD^2 - L*sin(a + b)*aDD ...
    - L*sin(a + b)*bDD - L*cos(a)*aD^2 - (L*cos(a + b + c)*aD^2)/2 - ...
    (L*cos(a + b + c)*bD^2)/2 - (L*cos(a + b + c)*cD^2)/2 - L*cos(a + b + c)*aD*bD ...
    - L*cos(a + b + c)*aD*cD - L*cos(a + b + c)*bD*cD - 2*L*cos(a + b)*aD*bD; 
ay1 = L*sin(a)*aD^2 - L*cos(a)*aDD + (L*sin(a + b + c)*aD^2)/2 + (L*sin(a + b + c)*bD^2)/2 ...
    + (L*sin(a + b + c)*cD^2)/2 - (L*cos(a + b + c)*aDD)/2 - (L*cos(a + b + c)*bDD)/2 - ...
    (L*cos(a + b + c)*cDD)/2 + L*sin(a + b)*aD^2 + L*sin(a + b)*bD^2 - L*cos(a + b)*aDD - ...
    L*cos(a + b)*bDD + L*sin(a + b + c)*aD*bD + L*sin(a + b + c)*aD*cD + ...
    L*sin(a + b + c)*bD*cD + 2*L*sin(a + b)*aD*bD;

Cx = M*xDD-L/2*M*(aD^2*cos(a)+aDD*sin(a)); 
Cy = (-L*M/2)*(aD^2*sin(a)-aDD*cos(a))+M*g; 

eqn1 = (-L/(2*I))*(M*g*cos(a)+Cy*cos(a)-Cx*sin(a))-aDD; 
eqn2 = (u-Cx)/m-xDD; 

sol = solve([eqn1, eqn2], [xDD, aDD]); 

%% State Declarations
state = [x ; xD ; a ; aD];

f = [xD ; simplify(sol.xDD) ; aD ; simplify(sol.aDD)]; 

%% Linearized System Definitions
A = ACalc(f, state, u, operate); 
B = BCalc(f, state, u, operate); 
C = [1 0 0 0 ; 0 0 1 0]; 
D = [1 0 ; 0 1]; 

%% LQR Controller Definition
Q = diag([100, 10, 100, 10]); 
R = 0.1; 
K = lqr(A, B, Q, R, 0); 

%% Time Definition
t = 0:Ts:finalTime; 
N = length(t); 

%% Covariances
Qekf = diag([1e-5, 1e-3, 1e-5, 1e-3]);
Rekf = diag([1e-3, 1e-3]); 

%% Preallocation
xTrue = zeros(4, N);
xHat  = zeros(4, N);
yMeas = zeros(2, N);
uHist = zeros(1, N);

%% Operation State
xe = [0 ; 0 ; pi/2 ; 0]; 

%% Initial States
xTrue(:, 1) = [0 ; 0 ; pi/2+.6 ; 0]; 
xHat(:, 1) = [0 ; 0 ; pi/2 ; 0]; 
P = diag([0.1, 0.1, 0.1, 0.1]);

v0 = mvnrnd([0 0], Rekf).';
yMeas(:,1) = C*xTrue(:,1) + v0;

%% Main Loop
for k=1:N-1
    % Input Calculation
    u = -K*(xHat(:, k) - xe); 
    uHist(k) = u; 

    % Propogation
    w = mvnrnd([0 0 0 0], Qekf).';
    [~, Xtmp] = ode45(@(tt,xx) pendDerr(xx, u), [t(k) t(k+1)], xTrue(:,k));
    xTrue(:,k+1) = Xtmp(end,:).' + w;
    
    % Measurement
    v = mvnrnd([0 0], Rekf).'; 
    yMeas(:, k+1) = C*xTrue(:, k + 1) + v; 

    % EKF Prediction
    xPred = xHat(:,k) + Ts*pendDerr(xHat(:,k), u);
    Ac = numericalJacobian(@(x) pendDerr(x, u), xHat(:,k));
    F = eye(4) + Ts*Ac;

    PPred = F*P*F' + Qekf;

    % EKF Correction
    H = C; 

    yPred = C*xPred; 
    innov = yMeas(:,k+1) - yPred;

    S = H*PPred*H' + Rekf;

    Kgain = PPred*H'/S;

    xHat(:,k+1) = xPred + Kgain*innov;
    P = (eye(4) - Kgain*H)*PPred;
end 

uHist(N) = uHist(N-1);

%% Theta over time plot
figure(1); 
hold on
grid on
plot(t, xTrue(3, :)*180/pi); 
plot(t, xHat(3, :)*180/pi); 
title("Inverted Pendulum"); 
xlabel("Time (s)"); 
ylabel("Angle (degrees)"); 
legend("True State", "Estimated State"); 
hold off

%% Cart Position
figure(2); 
hold on
grid on
plot(t, xTrue(1, :)); 
plot(t, xHat(1, :)); 
xlim([min(t) max(t)]); 
title("Cart Position"); 
xlabel("Time (s)"); 
ylabel("Position (m)"); 
legend("True State", "Estimated State"); 
hold off

%% Cart Angular Velocity
figure(3);
hold on
grid on
plot(t, xTrue(4, :)); 
plot(t, xHat(4, :)); 
xlabel("Time (s)"); 
ylabel("Angular Velocity (rad/s)"); 
legend("True State", "Estimated State"); 
title("Cart Angular Velocity"); 
hold off

%% Cart Velocity
figure(4); 
hold on
grid on
plot(t, xTrue(2, :)); 
plot(t, xHat(2, :)); 
xlabel("Time (s)"); 
ylabel("Velocity (m/s)"); 
legend("True State", "Estimated State"); 
title("Cart Velocity"); 
hold off

%% Input Force
figure(5); 
hold on
grid on
plot(t, uHist); 
xlabel("Time (s)"); 
ylabel("Force (N)"); 
title("Input Force over Time"); 
hold off

%% Inverted Pendulum on cart animation
figure(6);
hold on
grid on
axis equal
title("Inverted pendulum"); 
xlabel("X Location (m)"); 
ylabel("Y Location (m)"); 

% Axis Definitions
xmin = min(xTrue(1, :)) - L - .5; 
xmax = max(xTrue(1, :)) + L + .5; 
ymin = -10; 
ymax = 10; 
axis([xmin xmax ymin ymax]);

% Ground
plot([xmin xmax], [0 0], 'g', 'LineWidth', 2); 

% Initial Cart Location
cartX = xTrue(1, 1) - cartW/2; 
cartY = wheelR; 

% Cart Definitions
cart = rectangle('Position', [cartX, cartY, cartW, cartH], 'FaceColor', [.2 .6 .8], 'EdgeColor', 'k');
wheel1 = rectangle('Position', [xTrue(1, 1)-cartW/4-wheelR, 0, 2*wheelR, 2*wheelR], 'Curvature', [1 1], 'FaceColor', 'k'); 
wheel2 = rectangle('Position', [xTrue(1, 1)+cartW/4-wheelR, 0, 2*wheelR, 2*wheelR], 'Curvature', [1 1], 'FaceColor', 'k'); 

% Pivot Location
pivotX = xTrue(1, 1); 
pivotY = cartY + cartH/2; 

% Pendulum Location
pendX = pivotX + L*cos(xTrue(3, 1)); 
pendY = pivotY + L*sin(xTrue(3, 1)); 

% Rod Definition
rod = plot([pivotX pendX], [pivotY pendY], 'r', 'LineWidth', 3); 

% Animation Loop
for i = 1:length(t)-1
    % Cart update
    cartX = xTrue(1, i) - cartW/2;
    cartY = wheelR;
    cart.Position = [cartX, cartY, cartW, cartH];

    % Wheels update
    wheel1.Position = [xTrue(1, i)-cartW/4-wheelR, 0, 2*wheelR, 2*wheelR];
    wheel2.Position = [xTrue(1, i)+cartW/4-wheelR, 0, 2*wheelR, 2*wheelR];

    % Pivot location
    pivotX = xTrue(1, i);
    pivotY = cartY + cartH/2;

    % Pendulum tip with theta measured from +x axis
    pendX = pivotX + L*cos(xTrue(3, i));
    pendY = pivotY + L*sin(xTrue(3, i));

    % Update graphics
    set(rod, 'XData', [pivotX pendX], 'YData', [pivotY pendY]);

    drawnow;
    pause(Ts); 
end
hold off

%% A calculation
function A = ACalc(f, state, u, operate)
    A = zeros(4); 
    for index1 = 1:4
        for index2 = 1:4
            A(index1, index2) = subs(subs(subs(subs(subs(diff(f(index1), state(index2)), state(1), operate(1)), state(2), operate(2)), state(3), operate(3)), state(4), operate(4)), u, operate(5)); 
        end
    end
end

%% B calculation
function B = BCalc(f, state, u, operate)
    B = zeros(4, 1); 
    for index = 1:4
        B(index) = subs(subs(subs(subs(subs(diff(f(index), u), state(1), operate(1)), state(2), operate(2)), state(3), operate(3)), state(4), operate(4)), u, operate(5));
    end
end

%% Pendulum Derivative Function
function dx = pendDerr(x, u)
    x1 = x(1); 
    x2 = x(2); 
    x3 = x(3); 
    x4 = x(4); 
    dx = [x2 ...
        ; (160*cos(x3)*x4^2+32*u-981*sin(2*x3))/(100*cos(x3)^2+380) ...
        ; x4 ...
        ; (50*cos(x3)*sin(x3)*x4^2 - 2943*cos(x3) + 10*u*sin(x3))/(10*(5*cos(x3)^2 + 19))]; 
end

%% Jacobian Calculation
function J = numericalJacobian(fun, x)
    n = length(x);
    fx = fun(x);
    m = length(fx);
    J = zeros(m, n);

    epsVal = 1e-6;

    for i = 1:n
        dx = zeros(n,1);
        dx(i) = epsVal;

        J(:,i) = (fun(x + dx) - fun(x - dx)) / (2*epsVal);
    end
end
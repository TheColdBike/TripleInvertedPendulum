clear; clc; 

p.M = .969; % Cart Mass (kg)
p.m1 = .204; % Pendulum 1 Mass (kg)
p.m2 = .15; % Pendulm 2 Mass (kg)
p.m3 = .07; % Pendulum 3 Mass (kg)
p.L1 = .25; % Pendulum 1 Length (m)
p.L2 = .2; % Pendulum 2 Length (m)
p.L3 = .15; % Pendulum 3 Length (m)
p.I1 = .001688; % Pendulum 1 Inertia (kg*m^2)
p.I2 = .0011; % Pendulum 2 Inertia (kg*m^2)
p.I3 = .0003; % Pendulum 3 Inertia (kg*m^2)
p.d1 = .5; % Cart Damping Ratio
p.d2 = .01; % Pendulum 1 Damping Ratio
p.d3 = .01; % Pendulum 2 Damping Ratio
p.d4 = .01; % Pendulum 3 Damping Ratio
p.Wr = .04; % Wheel Radius
p.Ch = .2; % Cart Height
p.g = 9.81; % Gravity

timeStep = 0.01; 
endTime = 20; 
dynamicsOnly = 1; 

% Symbolic Variables
syms x a b c xD aD bD cD xDD aDD bDD cDD real; 
syms M m1 m2 m3 L1 L2 L3 I1 I2 I3 d1 d2 d3 d4 g u real

% State Vector
q = [x ; xD ; a ; aD ; b ; bD ; c ; cD]; 

% Symbolic M matrix for the Triple Inverted Pendulum
M_ = [(M + m1 + m2 + m3) (- (m3*(2*L2*cos(a + b) + 2*L1*cos(a) + L3*cos(a + b + c)))/2 - (m2*(L2*cos(a + b) + 2*L1*cos(a)))/2 - (L1*m1*cos(a))/2) (- (m3*(2*L2*cos(a + b) + L3*cos(a + b + c)))/2 - (L2*m2*cos(a + b))/2) (-(L3*m3*cos(a + b + c))/2) ; ...
(- (L3*m3*cos(a + b + c))/2 - (L2*m2*cos(a + b))/2 - L2*m3*cos(a + b) - (L1*m1*cos(a))/2 - L1*m2*cos(a) - L1*m3*cos(a)) (I1 + I2 + I3 + (L1^2*m1)/8 + L1^2*m2 + L1^2*m3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + (L1^2*m1*cos(2*a))/8 + L1*L2*m2*cos(b) + 2*L1*L2*m3*cos(b) + L2*L3*m3*cos(c) + L1*L3*m3*cos(b + c)) (I2 + I3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + (L1*L2*m2*cos(b))/2 + L1*L2*m3*cos(b) + L2*L3*m3*cos(c) + (L1*L3*m3*cos(b + c))/2) (I3 + (L3^2*m3)/4 + (L2*L3*m3*cos(c))/2 + (L1*L3*m3*cos(b + c))/2) ; ...
(- (L3*m3*cos(a + b + c))/2 - (L2*m2*cos(a + b))/2 - L2*m3*cos(a + b)) (I2 + I3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + (L1*L2*m2*cos(b))/2 + L1*L2*m3*cos(b) + L2*L3*m3*cos(c) + (L1*L3*m3*cos(b + c))/2) (I2 + I3 + (L2^2*m2)/4 + L2^2*m3 + (L3^2*m3)/4 + L2*L3*m3*cos(c)) ((m3*L3^2)/4 + (L2*m3*cos(c)*L3)/2 + I3) ; ...
(-(L3*m3*cos(a + b + c))/2) (I3 + (L3^2*m3)/4 + (L2*L3*m3*cos(c))/2 + (L1*L3*m3*cos(b + c))/2) ((m3*L3^2)/4 + (L2*m3*cos(c)*L3)/2 + I3) ((m3*L3^2)/4 + I3)];  

% Symbolic rhs for the Triple Inverted Pendulum
rhs = [u - (m2*(2*L1*aD^2*sin(a) + L2*sin(a + b)*(aD + bD)^2))/2 - d1*xD - (m3*(2*L1*aD^2*sin(a) + L3*sin(a + b + c)*(aD + bD + cD)^2 + 2*L2*sin(a + b)*(aD + bD)^2))/2 - (L1*aD^2*m1*sin(a))/2 ; ...
(L2*g*m2*sin(a + b))/2 - (L1^2*m1*sin(2*a))/8 - aD*d2 + L2*g*m3*sin(a + b) + (L1*g*m1*sin(a))/2 + L1*g*m2*sin(a) + L1*g*m3*sin(a) + (L1^2*aD^2*m1*sin(2*a))/8 + (L3*g*m3*sin(a + b + c))/2 + (L1*L3*bD^2*m3*sin(b + c))/2 + (L1*L3*cD^2*m3*sin(b + c))/2 + (L1*L2*bD^2*m2*sin(b))/2 + L1*L2*bD^2*m3*sin(b) + (L2*L3*cD^2*m3*sin(c))/2 + L2*L3*bD*cD*m3*sin(c) + L1*L3*aD*bD*m3*sin(b + c) + L1*L3*aD*cD*m3*sin(b + c) + L1*L3*bD*cD*m3*sin(b + c) + L1*L2*aD*bD*m2*sin(b) + 2*L1*L2*aD*bD*m3*sin(b) + L2*L3*aD*cD*m3*sin(c) ; ...
(L2*g*m2*sin(a + b))/2 - bD*d3 + L2*g*m3*sin(a + b) + (L3*g*m3*sin(a + b + c))/2 - (L1*L3*aD^2*m3*sin(b + c))/2 - (L1*L2*aD^2*m2*sin(b))/2 - L1*L2*aD^2*m3*sin(b) + (L2*L3*cD^2*m3*sin(c))/2 + L2*L3*bD*cD*m3*sin(c) + L2*L3*aD*cD*m3*sin(c) ; ...
(L3*g*m3*sin(a + b + c))/2 - cD*d4 - (L1*L3*aD^2*m3*sin(b + c))/2 - (L2*L3*aD^2*m3*sin(c))/2 - (L2*L3*bD^2*m3*sin(c))/2 - L2*L3*aD*bD*m3*sin(c)]; 

% Functions for M and rhs
Mfunc = matlabFunction(M_, ...
    'Vars', {x, a, b, c, ...
             M, m1, m2, m3, I1, I2, I3, L1, L2, L3});

rhsfunc = matlabFunction(rhs, ...
    'Vars', {x, xD, a, aD, b, bD, c, cD, u, ...
             M, m1, m2, m3, I1, I2, I3, L1, L2, L3, g, d1, d2, d3, d4});

% Initial Conditions
q0 = [0 ; 0 ; .15 ; 0 ; 0 ; 0 ; 0 ; 0]; 
tspan = 0:timeStep:endTime; 
N = length(tspan); 

% Linearized System Definitions
qop = [0; 0; 0; 0; 0; 0; 0; 0];
uop = 0;

A = ACalcNumerical(qop, uop, p, Mfunc, rhsfunc);
B = BCalcNumerical(qop, uop, p, Mfunc, rhsfunc);
C = [1 0 0 0 0 0 0 0; 0 0 1 0 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 1 0]; 
D = eye(4); 

%% LQR Controller Definition
Q = diag([100, 10, 25, 10, 25, 10, 25, 10]); 
R = 50; 
K = lqr(A, B, Q, R, 0); 

% Preallocation
qTrue = zeros(8, N); 
uHist = zeros(1, N); 

qTrue(:, 1) = q0; 
uHist(1) = 0; 
operationState = [0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0]; 

xHat = [0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0 ; 0];
P = 0.1*eye(8);

Qukf = 1e-4*eye(8);   % process noise covariance
Rukf = diag([1e-3, 1e-3, 1e-3, 1e-3]);  % measurement noise covariance

xHatHist = zeros(8, N); 
xHatHist(:, 1) = xHat; 

for k = 1:N-1
    u_k = (-K * (qTrue(:,k) - qop) + randn)*dynamicsOnly;
    uHist(k) = u_k;

    [~, qtmp] = ode45(@(tt,xx) diffQ(tt, xx, p, Mfunc, rhsfunc, u_k), [tspan(k) tspan(k+1)], qTrue(:,k));
    qTrue(:,k+1) = qtmp(end,:).';

    % Generate noisy measurement
    yMeas = measurementModelUKF(qTrue(:,k+1)) + chol(Rukf,'lower')*randn(4,1);

    % UKF predict/update
    [xHat, P] = ukfStep(xHat, P, u_k, yMeas, timeStep, Qukf, Rukf, p, Mfunc, rhsfunc);

    xHatHist(:,k+1) = xHat;
end
uHist(N) = uHist(N-1); 

% Cart Position
figure; 
hold on
grid on
plot(tspan, qTrue(1, :)); 
plot(tspan, xHatHist(1, :)); 
legend("True State", "UKF"); 
title("x over time"); 
ylabel("x (m)"); 
xlabel("t (s)"); 
hold off

% Cart Velocity
figure; 
hold on
grid on
plot(tspan, qTrue(2, :)); 
plot(tspan, xHatHist(2, :)); 
legend("True State", "UKF"); 
title("Cart Velocity over time"); 
ylabel("Velocity (m/s)"); 
xlabel("t (s)"); 
hold off

% Angle 1
figure; 
hold on
grid on
plot(tspan, qTrue(3, :) * 180 / pi); 
plot(tspan, xHatHist(3, :)*180/pi); 
legend("True State", "UKF"); 
title("Angle 1 over time"); 
ylabel("Angle (rad)"); 
xlabel("t (s)"); 
hold off

% Angular Velocity 1
figure; 
hold on
grid on
plot(tspan, qTrue(4, :) * 180/pi); 
plot(tspan, xHatHist(4, :)*180/pi); 
legend("True State", "UKF"); 
title("Anglular Velocity 1 over time"); 
ylabel("Anglular Velocity (rad/s)"); 
xlabel("t (s)"); 
hold off

% Angle 2
figure; 
hold on
grid on
plot(tspan, qTrue(5, :) * 180/pi); 
plot(tspan, xHatHist(5, :) * 180/pi); 
legend("True State", "UKF"); 
title("Angle 2 over time"); 
ylabel("Angle (rad)"); 
xlabel("t (s)"); 
hold off

% Angular Velocity 2
figure; 
hold on
grid on
plot(tspan, qTrue(6, :) * 180/pi); 
plot(tspan, xHatHist(6, :) * 180/pi); 
legend("True State", "UKF"); 
title("Anglular Velocity 2 over time"); 
ylabel("Anglular Velocity (rad/s)"); 
xlabel("t (s)"); 
hold off

% Angle 3
figure; 
hold on
grid on
plot(tspan, qTrue(7, :) * 180/pi); 
plot(tspan, xHatHist(7, :) * 180/pi); 
legend("True State", "UKF"); 
title("Angle 3 over time"); 
ylabel("Angle (rad)"); 
xlabel("t (s)"); 
hold off

% Angular Velocity 3
figure; 
hold on
grid on
plot(tspan, qTrue(8, :) * 180/pi); 
plot(tspan, xHatHist(8, :) * 180/pi); 
legend("True State", "UKF"); 
title("Anglular Velocity 3 over time"); 
ylabel("Anglular Velocity (rad/s)"); 
xlabel("t (s)"); 
hold off

% Input Force
figure; 
hold on
grid on
plot(tspan, uHist); 
title("Input Force over Time"); 
ylabel("Input Force (N)"); 
xlabel("t (s)"); 
hold off

% Animation call
animateTripleCartPendulum(tspan, qTrue, p);

function dq = diffQ(t, q, p, Mfunc, rhsfunc, ufunc)
    x = q(1); 
    xD = q(2); 
    a = q(3); 
    aD = q(4); 
    b = q(5); 
    bD = q(6); 
    c = q(7); 
    cD = q(8); 

    u = ufunc; 

    % Evaluate mass matrix
    Mnum = Mfunc(x, a, b, c, ...
                   p.M, p.m1, p.m2, p.m3, p.I1, p.I2, p.I3, p.L1, p.L2, p.L3);

    % Evaluate forcing vector
    rhsnum = rhsfunc(x, xD, a, aD, b, bD, c, cD, u, ...
                       p.M, p.m1, p.m2, p.m3, p.I1, p.I2, p.I3, p.L1, p.L2, p.L3, ...
                       p.g, p.d1, p.d2, p.d3, p.d4);

    % Solve for accelerations
    qDD = Mnum \ rhsnum;

    % State Derivatives
    dq = [xD; qDD(1); aD; qDD(2); bD; qDD(3); cD; qDD(4)];
end

function animateTripleCartPendulum(tspan, qTrue, p)

    % States
    x = qTrue(1, :);
    a = qTrue(3, :);
    b = qTrue(5, :);
    c = qTrue(7, :);

    % Link lengths
    L1 = p.L1;
    L2 = p.L2;
    L3 = p.L3;

    % Cart dimensions
    cartW = .3;
    cartH = p.Ch;
    wheelR = p.Wr;

    % Joint positions
    x0 = x;
    y0 = zeros(size(x)) + 0.5*cartH + 2*wheelR;

    x1 = x0 - L1*sin(a);
    y1 = y0 + L1*cos(a);

    x2 = x1 - L2*sin(a + b);
    y2 = y1 + L2*cos(a + b);

    x3 = x2 - L3*sin(a + b + c);
    y3 = y2 + L3*cos(a + b + c);

    % Figure
    figure;
    hold on;
    axis equal;
    grid on;

    % Plot bounds
    xmin = min([x0 x1 x2 x3]) - 1;
    xmax = max([x0 x1 x2 x3]) + 1;
    ymin = min([y0 y1 y2 y3]) - 1;
    ymax = max([y0 y1 y2 y3]) + 1;

    xlim([xmin xmax]);
    ylim([ymin ymax]);

    % Ground line
    plot([xmin xmax], [0 0], 'g', 'LineWidth', 1.5);

    % Initial cart
    cartX = x0(1) - cartW/2;
    cartY = 2*wheelR;
    cartPatch = rectangle('Position', [cartX, cartY, cartW, cartH], ...
                          'Curvature', 0.05, 'FaceColor', [0.2 0.6 0.8]);

    % Wheels
    th = linspace(0, 2*pi, 100);
    wheel1 = plot(x0(1)-cartW/4 + wheelR*cos(th), wheelR + wheelR*sin(th), ...
                  'b', 'LineWidth', 1.5);
    wheel2 = plot(x0(1)+cartW/4 + wheelR*cos(th), wheelR + wheelR*sin(th), ...
                  'b', 'LineWidth', 1.5);

    % Links
    link1 = plot([x0(1) x1(1)], [y0(1) y1(1)], 'LineWidth', 3);
    link2 = plot([x1(1) x2(1)], [y1(1) y2(1)], 'LineWidth', 3);
    link3 = plot([x2(1) x3(1)], [y2(1) y3(1)], 'LineWidth', 3);

    % Joints / masses
    joint0 = plot(x0(1), y0(1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
    joint1 = plot(x1(1), y1(1), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
    joint2 = plot(x2(1), y2(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
    joint3 = plot(x3(1), y3(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');

    titleHandle = title(sprintf('t = %.2f s', tspan(1)));

    % Animate
    for k = 1:length(tspan)

        % Cart
        cartX = x0(k) - cartW/2;
        set(cartPatch, 'Position', [cartX, cartY, cartW, cartH]);

        % Wheels
        set(wheel1, 'XData', x0(k)-cartW/4 + wheelR*cos(th), ...
                    'YData', wheelR + wheelR*sin(th));
        set(wheel2, 'XData', x0(k)+cartW/4 + wheelR*cos(th), ...
                    'YData', wheelR + wheelR*sin(th));

        % Links
        set(link1, 'XData', [x0(k) x1(k)], 'YData', [y0(k) y1(k)]);
        set(link2, 'XData', [x1(k) x2(k)], 'YData', [y1(k) y2(k)]);
        set(link3, 'XData', [x2(k) x3(k)], 'YData', [y2(k) y3(k)]);

        % Joints
        set(joint0, 'XData', x0(k), 'YData', y0(k));
        set(joint1, 'XData', x1(k), 'YData', y1(k));
        set(joint2, 'XData', x2(k), 'YData', y2(k));
        set(joint3, 'XData', x3(k), 'YData', y3(k));

        set(titleHandle, 'String', sprintf('t = %.2f s', tspan(k)));

        drawnow;
    end
end

% A Calculation
function A = ACalcNumerical(Xop, uop, p, Mfunc, rhsfunc)
    fX = @(X) stateDerivative(X, uop, p, Mfunc, rhsfunc);
    A = numericalJacobian(fX, Xop);
end

% B Calculation
function B = BCalcNumerical(Xop, uop, p, Mfunc, rhsfunc)
    epsVal = 1e-6;
    B = (stateDerivative(Xop, uop + epsVal, p, Mfunc, rhsfunc) - ...
         stateDerivative(Xop, uop - epsVal, p, Mfunc, rhsfunc)) / (2*epsVal);
end

% Jacobian Calculation
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

function dx = stateDerivative(X, u, p, Mfunc, rhsfunc)

    x  = X(1);
    xD = X(2);
    a  = X(3);
    aD = X(4);
    b  = X(5);
    bD = X(6);
    c  = X(7);
    cD = X(8);

    Mnum = Mfunc(x, a, b, c, ...
        p.M, p.m1, p.m2, p.m3, p.I1, p.I2, p.I3, p.L1, p.L2, p.L3);

    rhsnum = rhsfunc(x, xD, a, aD, b, bD, c, cD, u, ...
        p.M, p.m1, p.m2, p.m3, p.I1, p.I2, p.I3, p.L1, p.L2, p.L3, ...
        p.g, p.d1, p.d2, p.d3, p.d4);

    qDD = Mnum \ rhsnum;

    dx = [xD;
          qDD(1);
          aD;
          qDD(2);
          bD;
          qDD(3);
          cD;
          qDD(4)];
end

function Xnext = processModelUKF(X, u, Ts, p, Mfunc, rhsfunc)

    [~, Xtmp] = ode45(@(t,xx) stateDerivative(xx, u, p, Mfunc, rhsfunc), [0 Ts], X);
    Xnext = Xtmp(end,:).';

end

function y = measurementModelUKF(X)
    y = [X(1);
         X(3);
         X(5);
         X(7)];
end

function [xHatNew, PNew] = ukfStep(xHat, P, u, y, Ts, Q, R, p, Mfunc, rhsfunc)

    n = length(xHat);

    alpha = 1e-3;
    beta  = 2;
    kappa = 0;

    lambda = alpha^2*(n + kappa) - n;
    gamma = sqrt(n + lambda);

    % Weights
    Wm = zeros(2*n+1,1);
    Wc = zeros(2*n+1,1);

    Wm(1) = lambda/(n+lambda);
    Wc(1) = lambda/(n+lambda) + (1 - alpha^2 + beta);

    for i = 2:2*n+1
        Wm(i) = 1/(2*(n+lambda));
        Wc(i) = 1/(2*(n+lambda));
    end

    % Sigma points
    S = chol(P, 'lower');
    Xsigma = zeros(n, 2*n+1);
    Xsigma(:,1) = xHat;

    for i = 1:n
        Xsigma(:,i+1)   = xHat + gamma*S(:,i);
        Xsigma(:,i+1+n) = xHat - gamma*S(:,i);
    end

    % Predict sigma points through process model
    XsigmaPred = zeros(n, 2*n+1);
    for i = 1:2*n+1
        XsigmaPred(:,i) = processModelUKF(Xsigma(:,i), u, Ts, p, Mfunc, rhsfunc);
    end

    % Predicted mean
    xPred = zeros(n,1);
    for i = 1:2*n+1
        xPred = xPred + Wm(i)*XsigmaPred(:,i);
    end

    % Predicted covariance
    PPred = zeros(n);
    for i = 1:2*n+1
        dx = XsigmaPred(:,i) - xPred;
        PPred = PPred + Wc(i)*(dx*dx.');
    end
    PPred = PPred + Q;

    % Predicted measurements
    m = length(y);
    Ysigma = zeros(m, 2*n+1);
    for i = 1:2*n+1
        Ysigma(:,i) = measurementModelUKF(XsigmaPred(:,i));
    end

    yPred = zeros(m,1);
    for i = 1:2*n+1
        yPred = yPred + Wm(i)*Ysigma(:,i);
    end

    % Measurement covariance
    Pyy = zeros(m);
    for i = 1:2*n+1
        dy = Ysigma(:,i) - yPred;
        Pyy = Pyy + Wc(i)*(dy*dy.');
    end
    Pyy = Pyy + R;

    % Cross covariance
    Pxy = zeros(n,m);
    for i = 1:2*n+1
        dx = XsigmaPred(:,i) - xPred;
        dy = Ysigma(:,i) - yPred;
        Pxy = Pxy + Wc(i)*(dx*dy.');
    end

    % Gain and update
    Kgain = Pxy / Pyy;

    xHatNew = xPred + Kgain*(y - yPred);
    PNew = PPred - Kgain*Pyy*Kgain.';
end
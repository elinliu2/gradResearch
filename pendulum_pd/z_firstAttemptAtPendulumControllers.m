% From ECE 484 Pset 7 

% For this example we are considering a nonlinear pendulum system
% We set the desired location of the pendulum to be pi/4

% For instability, let's consider the case where if the angle is greater
% than pi/4 it is unstable. (This means we start in the stable region)
% We'll set the stable equilibrium point to be 0. 
% If we have time we can also consider if the absolute value of angular 
% velocity is greater than pi/2 we're unstable. And if we're feeling real
% ballsy I guess we could try to combine those somehow

% [x1, x2, time] = PendulumSim(@ConstantU, "ConstantU", "blue");
% StatePlot(x1, x2, time, 1, "ConstantU");

[x1, x2, time] = PendulumSim(@PIDCtrl, "PIDCtrl", [0.5 0 0.5]);
StatePlot(x1, x2, time, 1, "PIDCtrl");

[x1, x2, time] = PendulumSim(@WeightedAvgTheta, "WeightedAvgTheta", "magenta");
StatePlot(x1, x2, time, 2, "WeightedAvgTheta");

[x1, x2, time] = PendulumSim(@WeightedAvgThetaLog, "WeightedAvgThetaLog", "red");
StatePlot(x1, x2, time, 3, "WeightedAvgThetaLog");

[x1, x2, time] = PendulumSim(@WeightedAvgWhetaLog, "WeightedAvgWhetaLog", "black");
StatePlot(x1, x2, time, 4, "WeightedAvgWhetaLog");

[x1, x2, time] = PendulumSim(@minimalIntervention, "minimalIntervention", [0 0.5 0.5]);
StatePlot(x1, x2, time, 5, "minimalIntervention");

[x1, x2, time] = PendulumSim(@ScaledWeightedAvgTheta, "ScaledWeightedAvgTheta", [0.5 0.5 0]);
StatePlot(x1, x2, time, 6, "ScaledWeightedAvgTheta");

[x1, x2, time] = PendulumSim(@bump, "bump", [0.5 0.75 0.75]);
StatePlot(x1, x2, time, 7, "bump");

[x1, x2, time] = PendulumSim(@optimization, "optimization", [1 0.75 0.75]);
StatePlot(x1, x2, time, 8, "optimization");

function [x1, x2, time] = PendulumSim(u, plotTitle, color)
    % constants
    g = 9.8;
    l = g/2;
    m = 2/g;
    D = 0.1;
    ref = pi/4;

    % setup arrays 
    T = 0.02;
    t_final = 20;
    x1 = zeros([1 t_final/T]);
    x2 = zeros([1 t_final/T]);
    e = zeros([1 t_final/T]);
    time = zeros([1 t_final/T]);

    x = 0;
    y = 0;
    figure(1)
    hold on
    scatter(l*sin(ref), -l*cos(ref), 25, "green", "filled");
    x_lim = linspace(0,9.8/2*sin(pi/4));
    y_lim = -x_lim*tan(pi/4);
    plot(x_lim, y_lim, "green");
    sim = scatter(x, y, 20, color, "filled");
    title(plotTitle);
    xlabel("x");
    ylabel("y");
    xlim([-6 6])
    ylim([-6 0])

    % initial conditions
    x1(1) = 0;
    x2(1) = 0;
    for i = 2:t_final/T
        e(i) = ref - x1(i-1);

        ctrlSig = u(x1(i-1), x2(i-1), e(i), (e(i)-e(i-1))/T, sum(e));

        x1(i) = x1(i-1) + x2(i-1)*T;
        x2(i) = x2(i-1) + (-g/l*sin(x1(i-1))-D/(m*l)*x2(i-1)+1/(m*l)*ctrlSig)*T;

        time(i) = i*T;
        sim.XData = l*sin(x1(i));
        sim.YData = -l*cos(x1(i));
        drawnow
        pause(0.001);
    end
end

function StatePlot(x1, x2, time, index, plotTitle)
    figure(2)
    subplot(4, 2, index);
    scatter(x1, x2, 1, time);
    title(plotTitle);
    xlabel("x1");
    ylabel("x2");
end

function u = ConstantU(x1, x2, e, ed, ei)
    u = 1;
end

function u = PIDCtrl(x1, x2, e, ed, ei)
    p = 2;
    i = 0.05;
    d = 2;
    u = p * e + i * ei + d * ed;
end

function u = WeightedAvgTheta(x1, x2, e, ed, ei)
    boundaryTheta = pi/4;
    dist = abs(x1-boundaryTheta); % dist from RoA boundary

     u = (1/(dist + 1)) * stabilityControllerPos(x1, x2)  ...
        + (1 - 1/(dist + 1)) * PIDCtrl(x1, x2, e, ed, ei);
end

function u = ScaledWeightedAvgTheta(x1, x2, e, ed, ei)
    boundaryTheta = pi/4;
    dist = abs(x1-boundaryTheta); % dist from RoA boundary

     u = (0.001/(dist + 0.001)) * stabilityControllerPos(x1, x2)  ...
        + (1 - 0.001/(dist + 0.001)) * PIDCtrl(x1, x2, e, ed, ei);
end

function u = WeightedAvgThetaLog(x1, x2, e, ed, ei)
    boundaryTheta=pi/4;
    dist = abs(x1-boundaryTheta); % dist from RoA boundary

     u = (1 - log(dist+1)) * stabilityControllerPos(x1, x2)  ...
        + (log(dist+1)) * PIDCtrl(x1, x2, e, ed, ei);
end

function u = stabilityControllerPos(x1, x2)
    k = 1;
    u = - k * x1;
end

function u = WeightedAvgWheta(x1, x2, e, ed, ei)
    boundaryWheta=1;
    dist = abs(x2-boundaryWheta); % dist from RoA boundary

    u = (1/(dist + 1)) * stabilityControllerWheta(x1, x2)  ...
        + (1 - 1/(dist + 1)) * PIDCtrl(x1, x2, e, ed, ei);
end

function u = WeightedAvgWhetaLog(x1, x2, e, ed, ei)
    boundaryWheta=1;
    dist = abs(x2-boundaryWheta); % dist from RoA boundary

    u = (1/(dist + 1)) * stabilityControllerWheta(x1, x2)  ...
        + (1 - 1/(dist + 1)) * PIDCtrl(x1, x2, e, ed, ei);
end

function u = stabilityControllerWheta(x1, x2)
    k = 1;
    u = - k * x2;
end

function u = minimalIntervention(x1, x2, e, ed, ei)
    % constants
    g = 9.8;
    l = g/2;
    m = 2/g;
    boundaryTheta = pi/4;

    % this is the maximum u our controller can apply 
    maxU = 10; 

    % Are we going to hit the border? Check that x1 and x2 have the same
    % sign. If they have different signs, we can assume we're not going to
    % hit the border for now ...
    if (x1/x2 < 0) % They have different signs so we won't hit the border
        u = PIDCtrl(x1, x2, e, ed, ei);
        return
    end 
    
    % If we are going to hit the border, what control effort would we have 
    % to exert at the border for us to stop from going over

    % Conservation of energy ? 
    % ignore friction/damping for now :grimacing: 
    deltaWheta = abs(sqrt(2 * g * (cos(x1)-cos(pi/4))));

    finalWheta = x2 - deltaWheta;
    % What acceleration do we need to have to stay in the boundary region
    netAlpha = -1 * finalWheta^2 / (2 * (boundaryTheta - x1)); 
    % what torque do we need to apply to achieve that alpha?
    minU = m * l * netAlpha + m * g * sin(x1);

    if (abs(minU) >= maxU)
        u = minU ;
    else 
        % okay, we don't need to apply the max stabilizing control effort
        % yet. But then we should make sure we don't cross over the border
        % really slowly either. We want to limit what the PID controller
        % can do. 
        T = 0.02;
        u = PIDCtrl(x1, x2, e, ed, ei);
        alpha = (u - m * g * sin(x1))/(m*l);
        thetaF = x1 + x2 * 5 * T + 1/2 * alpha * T ^2;
        if (thetaF > boundaryTheta)
            u = minU;
        end
    end
end

function u = bump(x1, x2, e, ed, ei)
    boundaryTheta=pi/4;
    dist = abs(x1-boundaryTheta); % dist from RoA boundary

     u = -exp(1/(100*dist)) ...
        + PIDCtrl(x1, x2, e, ed, ei);
end

function u = optimization(x1, x2, e, ed, ei)
    g = 9.8;
    l = g/2;
    m = 2/g;
    D = 0.1;
    boundaryTheta=pi/4;
    T = 0.02;
    
    u_des = PIDCtrl(x1, x2, e, ed, ei);
    f = [1 0];
    % minimize t 
    % subject to 
    % t <= udes - uacc --> t + uacc <= udes
    % -t <= -udes + uacc --> - t - uacc <= - udes
    % 
    % x1 + (x2 + (-g/l*sin(x1)-D/(m*l)*x2+1/(m*l)*uacc)*T)*T <=
    % boundaryTheta
    % 
    b_3 = boundaryTheta - x1 - (x2 + (-g/l*sin(x1)-D/(m*l)*x2)*T)*T;
    A = [ 1 1; -1 -1; -1 0; 0 1/(m*l)*T*T ];
    b = [u_des -u_des  0 b_3];
    x = linprog(f,A,b);
    u = x(2);
end


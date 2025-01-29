addpath('./../utils/')
tjs = trajectorySensitivities;

% Sim Constants
SC.T = 1e-3;
SC.t_final = pi;
SC.start = -pi/2;
SC.stop = pi/2;
SC.step = 0.1;
SC.ref = 0;

% Model
pendulum = PendulumModel(0.2, 0.2);

initialState = [-2, -0.5];

numIterations = 15;
pid_p_tracker = pendulum.kp*ones(1, numIterations);
pid_d_tracker = pendulum.kd*ones(1, numIterations);
G_tracker = zeros(1, numIterations);
for counter=1:numIterations
    pid_p_tracker(counter) = pendulum.kp;
    pid_d_tracker(counter) = pendulum.kd;

    [~, stateProgression, ~] = trajectory(pendulum, SC, initialState, false);
    trajSens = tjs.calculateTrajectorySensitivities(pendulum, SC, stateProgression);
    G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
    trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(pendulum, SC, stateProgression);
    secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(pendulum, SC, stateProgression, trajSens, trajSensParam);

    tp = tjs.findtp(trajSens);
    dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, pendulum);

    scale = 0.1;
    normdGp = sqrt(dGp(1)^2 + dGp(2)^2);
    pendulum.kp = pendulum.kp + dGp(1)/normdGp*scale;
    pendulum.kd = pendulum.kd + dGp(2)/normdGp*scale;
    G_tracker(counter) = G;
end

splotting.uplup(1:counter, pid_p_tracker(1:counter), 2, "pid_p over iteration", "pid_p");
splotting.uplup(1:counter, pid_d_tracker(1:counter), 3, "pid_d over iteration", "pid_d");
splotting.uplup(1:counter, G_tracker(1:counter), 4, "G over iteration", "G");
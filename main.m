addpath('./utils')  
% Sim Constants
SC.T = 1e-2;
SC.t_final = 30;
SC.g = 9.8;
SC.l = SC.g/2;
SC.m = 2/SC.g;
SC.D = 0.1;

SC.start = -2*pi;
SC.stop = 2*pi;
SC.step = 0.1;
SC.ref = 0;

SC.u = 0;
SC.disturbance = 5e5;

tjs = trajectorySensitivities;
norms = normalStuff;
sCtrl = stonkCtrllers;

SC.pid_i = 0;
SC.pid_p = 0.867008;
SC.pid_d = 1.92167;

SC = derivatives_pdctrlmodel.assignDerivatives(SC);
SC.controller = @proposedCtrllers.ParameterizedPIDCtrl;

% plot_all_the_plots(SC);

initialState = [-2, -0.5];
[time, stateProgression, SC] = pendulum_trajectory_v2(SC, initialState, SC.controller, @stonkCtrllers.dummyD, true);
splotting.StatePlot(stateProgression(:,1), stateProgression(:,2), time, 1, "State Progression");
% splotting.uplup(time, stateProgression(:,1), 2, "x1 over time", "x1");
% splotting.uplup(time, stateProgression(:,2), 3, "x2 over time", "x2");


% old_pid_p = 0;
% old_pid_d = 0; 
% initialState = [-2, -0.5];
% 
% pid_p_tracker = SC.pid_p*ones(1, 15);
% pid_d_tracker = SC.pid_d*ones(1, 15);
% G_tracker = zeros(1, 10);
% for counter=1:20
%     pid_p_tracker(counter) = SC.pid_p;
%     pid_d_tracker(counter) = SC.pid_d;
% 
%     old_pid_p = SC.pid_p;
%     old_pid_d = SC.pid_d; 
%     [~, stateProgression, ~] = pendulum_trajectory_v2(SC, initialState, @proposedCtrllers.ParameterizedPIDCtrl, @sCtrl.dummyD, false);
%     trajSens = tjs.calculateTrajectorySensitivities(stateProgression, SC);
%     G = norms.infimumOfReciprocalOfElementOneNorm(trajSens);
%     trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(stateProgression, SC);
%     secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(stateProgression, trajSens, trajSensParam, SC);
% 
%     tp = tjs.findtp(trajSens);
%     dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, SC);
% 
%     scale = 0.1;
%     normdGp = sqrt(dGp(1)^2 + dGp(2)^2);
%     SC.pid_p = SC.pid_p + dGp(1)/normdGp*scale;
%     SC.pid_d = SC.pid_d + dGp(2)/normdGp*scale;
%     G_tracker(counter) = G;
% end
% 
% splotting.uplup(1:counter, pid_p_tracker(1:counter), 2, "pid_p over iteration", "pid_p");
% splotting.uplup(1:counter, pid_d_tracker(1:counter), 3, "pid_d over iteration", "pid_d");
% splotting.uplup(1:counter, G_tracker(1:counter), 4, "G over iteration", "G");


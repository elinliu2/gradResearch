addpath('./utils')  
% Sim Constants
SC.T = 1e-2;
SC.t_final = 20;
SC.g = 9.8;
SC.l = SC.g/2;
SC.m = 2/SC.g;
SC.D = 0.1;

SC.start = -2*pi;
SC.stop = 2*pi;
SC.step = 0.25;
SC.ref = 0;

SC.u = 0;
SC.disturbance = 0.5;

tjs = trajectorySensitivities;
norms = normalStuff;
sCtrl = stonkCtrllers;

SC.pid_i = 0;
SC.pid_p = 0.2;
SC.pid_d = 0.2;

SC = derivatives_pdctrlmodel.assignDerivatives(SC);

% plot_all_the_plots(SC);

SC.controller = @proposedCtrllers.ChangingPIDCtrl;

initialState.x1 = -0.483185;
initialState.x2 = 3.51681;
initialState.e = 0.483185;
initialState.ei = 0;
initialState.ed = -3.51681;

SC.G = 0;

[stateProgression, time, u, SC] = pendulum_trajectory(SC, initialState,SC.controller, @stonkCtrllers.dummyD, false);
[x1, x2] = splatting.extractX1andX2FromState(stateProgression);
splotting.StatePlot(x1, x2, time, 1, "State Progression");
splotting.uplup(time, x1, 2, "x1 over time", "x1");
splotting.uplup(time, x2, 3, "x2 over time", "x2");
splotting.uplup(time, u, 4, "u over time", "u");
splotting.uplup(time, SC.G, 5, "G over time", "G");

% SC.u = 0.5; 
% 
% [stateProgression, ~, ~, ~] = pendulum_trajectory(SC, initialState, @stonkCtrllers.dummyU, @stonkCtrllers.dummyD, false);
% trajSens2 = tjs.calculateTrajectorySensitivities(stateProgression, SC);

% [x1, x2] = splatting.extractX1andX2FromState(state);
% splotting.StatePlot(x1, x2, time, 1, "State Progression");
% splotting.uplup(time, u, 2, "u over time", "u");
% [x1, x2, time, u] = pendulum_trajectory(SC, 2.158, 1.408, @proposedCtrllers.idkWhatImDoingv2, @stonkCtrllers.dummyD, true);
% splotting.StatePlot(x1, x2, time, 3, "idkWhatImDoingv2 controller");
% splotting.uplup(time, x1, 4, "x1 over time", "x1");
% splotting.uplup(time, x2, 5, "x2 over time", "x2");
% splotting.uplup(time, u, 6, "u over time", "u");

% Gx = zeros(size(time,2));
% for i = 1:size(time,2)
%     SC.u = 0;
% 
%     % Gx
%     [x1gx, x2gx] = pendulum_trajectory(SC, x1(i), x2(i), @sCtrl.dummyU, @sCtrl.dummyD, false);
%     trajSens = deriv.calculateTrajectorySensitivities(x1gx, x2gx, SC);
%     Gx(i) = norms.infimumOfReciprocalOfOneNorm(trajSens);
% end
% 
% splotting.uplup(time, Gx, 7, "Gx over time", "Gx");


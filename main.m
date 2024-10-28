addpath('./utils')  
% Sim Constants
SC.T = 1e-2;
SC.t_final = 10;
SC.g = 9.8;
SC.l = SC.g/2;
SC.m = 2/SC.g;
SC.D = 0.1;
SC.u = 0;

SC.start = 1;
SC.stop = 3;
SC.step = 0.05;
SC.ref = 0;

SC.u = 0;
SC.disturbance = 0;

deriv = derivatives;
norms = normalStuff;
sCtrl = stonkCtrllers;

SC.disturbance = 0.00;

plot_all_the_plots(SC);
% [x1, x2, time, u] = pendulum_trajectory(SC, 1.70841, 2.15841, @stonkCtrllers.dummyU, @stonkCtrllers.dummyD, true);
% splotting.StatePlot(x1, x2, time, 1, "State Progression");
% splotting.uplup(time, u, 2, "u over time");
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


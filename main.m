addpath('./utils')  
addpath('./pendulum_pd')

% Sim Constants
SC.T = 1e-2;
SC.t_final = 30;
SC.start = -2*pi;
SC.stop = 2*pi;
SC.step = 0.5;
SC.ref = 0;

% Model
pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

% initialState = [-2, -0.5];
[time, stateProgression, pendulum] = trajectory(pendulum, SC, initialState, false);
% splotting.StatePlot(stateProgression(:,1), stateProgression(:,2), time, 1, "State Progression");
% splotting.uplup(time, stateProgression(:,1), 2, "x1 over time", "x1");
% splotting.uplup(time, stateProgression(:,2), 3, "x2 over time", "x2");


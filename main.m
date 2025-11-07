addpath(genpath('.'))
% warning('off', 'all')

% Sim Constants
SC.T = 1e-3;
SC.t_final = 20;

% Model
% pendulum = PendulumModel(0.2, 0.2);
% plot_all_the_plots(SC, pendulum);

[drone, z0, ref, dist] = initializeDrone(SC);

% tic;
% disp("find closest zbar")
% zbar = findClosestzbar(drone, SC, z0, ref, dist, z0);
% toc;
% 
% tic;
% disp("find optimal p")
% [z, p, p_iteration, dist_zbar] = findOptimalqGradAscent(drone, SC, z0, ref, dist, zbar);
% drone = assignP(drone, p);
% toc;

% t = 10000;
% dist.Fx(1:1+t) = 20*ones(1, 1+t);
% dist.Fy(1:1+t) = 20*ones(1, 1+t);
% dist.Fz(1:1+t) = -20*ones(1, 1+t);

% t = 1;
% dist.Fx(1:1+t) = -6*ones(1, 1+t);
% dist.Fy(1:1+t) = 6*ones(1, 1+t);

[time, stateProgression, drone, stable] = trajectory(drone, SC, z0, false, ref, dist);
% testTrajectorySensitivities(drone, SC, z0, stateProgression, ref, dist);
% trajSens = droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams(drone, SC, stateProgression, ref);
% secondTrajSensParam = droneTrajectorySensitivitiesHybrid.calculateSecondOrderTrajectorySensitivitiesParams(drone, SC, stateProgression, trajSens, ref);
% G = norms.infimumOfReciprocalOfOneNorm(trajSens);
splotting.trajectoryPlot3d(time, stateProgression(:, 1), stateProgression(:, 2), stateProgression(:, 3), 1, "Initial Parameters with No Disturbance State Progression");
% splotting.trajectoryVideo(time, stateProgression);

% figure(2)
% sgtitle('$p^k$ variation over iteration','interpreter','latex')
% for i = 1:18
%     subplot(9, 2, i);
%     p_iterate = p_iteration(i, 1:34);
%     plot(p_iterate)
%     xlabel('iteration','interpreter','latex')
%     ylabel('$p_{' + string(i) + '}$','interpreter','latex')
% 
% end
% 
% figure(3)
% plot(dist_zbar(1, 1:34))
% title('Distance to closest point, $z$, on RoA boundary from $z_0$','interpreter','latex')
% xlabel('iteration','interpreter','latex')
% ylabel('$\frac{1}{2}\|z-z_0\|^2$','interpreter','latex')



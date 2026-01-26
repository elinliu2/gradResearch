function [drone, x0, ref, dist] = initializeDrone(SC)
    % c1 = 10.5;
    % c2 = 2;
    % c3 = 10;
    % c4 = 2;
    % c5 = 2;
    % c6 = 2;
    % c7 = 4;
    % c8 = 2;
    % c9 = 1;
    % c10 = 0.5;
    % c11 = 1;
    % c12 = 0.5;
    % lambda1 = 1e-1;
    % lambda2 = 1e-1;
    % lambda3 = 1e-1;
    % lambda4 = 1e-1;
    % lambda5 = 1e-1;
    % lambda6 = 1e-1;

    c1 = 1;
    c2 = 1;
    c3 = 1;
    c4 = 1;
    c5 = 1;
    c6 = 1;
    c7 = 1;
    c8 = 1;
    c9 = 1;
    c10 = 1;
    c11 = 1;
    c12 = 1;
    lambda1 = 1;
    lambda2 = 1;
    lambda3 = 1;
    lambda4 = 1;
    lambda5 = 1;
    lambda6 = 1;

    m = 1.4; % kg
    g = 9.8;
    drone = droneModelIBHybrid(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g);
    
    ref.x = 1;
    ref.x_dot = 0;
    ref.x_ddot = 0;
    ref.x_dddot = 0;
    ref.x_ddddot = 0;
    ref.y = 1;
    ref.y_dot = 0;
    ref.y_ddot = 0;
    ref.y_dddot = 0;
    ref.y_ddddot = 0;
    ref.z = 100;
    ref.z_dot = 0;
    ref.z_ddot = 0;
    ref.z_dddot = 0;
    ref.z_ddddot = 0;
    ref.psi = 0;
    ref.psi_dot = 0;
    ref.psi_ddot = 0;
    ref.psi_dddot = 0;
    ref.psi_ddddot = 0;

    dist.Fx = zeros(1, floor(SC.t_final/SC.T) + 1);
    dist.Fy = zeros(1, floor(SC.t_final/SC.T) + 1);
    dist.Fz = zeros(1, floor(SC.t_final/SC.T) + 1);
    
    dist.Tx = zeros(1, floor(SC.t_final/SC.T) + 1);
    dist.Ty = zeros(1, floor(SC.t_final/SC.T) + 1);
    dist.Tz = zeros(1, floor(SC.t_final/SC.T) + 1);
    
    x0 = zeros(drone.numStates, 1);
    for i = 1:12
        x0(i) = 0;
    end
    x3 = 100; % Drone starts 100 m in the air
    x0(3) = x3; 
    e7 = ref.z - x0(3);
    e8 = c7*e7 + ref.z_dot + lambda4*x0(16) - x0(6);
    U_z = m/(cos(x0(7))*cos(x0(8))*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x0(16));
    x0(19) = U_z;
    % x0 = load("z0.mat").z0';
end
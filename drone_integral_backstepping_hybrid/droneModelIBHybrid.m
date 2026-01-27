classdef droneModelIBHybrid
    properties
        g;
        m;

        c1;
        c2;
        c3;
        c4;
        c5;
        c6;
        c7;
        c8;
        c9;
        c10;
        c11;
        c12;

        lambda1;
        lambda2;
        lambda3;
        lambda4;
        lambda5;
        lambda6;

        l; % l is the distance between Ob and each propeller; 
        K_psi; %  is the thrust-to-moment gain ;

        Jtheta;
        Jphi;
        Jpsi;

        disturbance;
        modelTrajectory;
        calculateTrajectorySensitivitiesParams;
        nudgeParams;

        numX;
        numY;
        numStates;
        numParams;
        initialdxdxo;

        disturbanceSignal;
        initializeStates;
    end
    
    methods
        function obj = droneModelIBHybrid(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g)
            obj.g = g;
            obj.m = m; % kg

            obj.l = 0.2; % m
            obj.K_psi = 4; % Nm

            obj.Jtheta = 0.03; % kgm^2
            obj.Jphi = 0.03; % kgm^2
            obj.Jpsi = 0.04; % kgm^2

            obj.c1 = c1;
            obj.c2 = c2;
            obj.c3 = c3;
            obj.c4 = c4;
            obj.c5 = c5;
            obj.c6 = c6;
            obj.c7 = c7;
            obj.c8 = c8;
            obj.c9 = c9;
            obj.c10 = c10;
            obj.c11 = c11;
            obj.c12 = c12;

            obj.lambda1 = lambda1;
            obj.lambda2 = lambda2;
            obj.lambda3 = lambda3;
            obj.lambda4 = lambda4;
            obj.lambda5 = lambda5;
            obj.lambda6 = lambda6;

            obj.disturbance = 0;
            
            obj.modelTrajectory = @trajectory;
            obj.calculateTrajectorySensitivitiesParams = @droneTrajectorySensitivitiesHybrid.calculateTrajectorySensitivitiesParams;
            obj.nudgeParams = @derivativesDroneIBHybrid.nudgeParams;

            obj.numX = 19;
            obj.numY = 4;
            obj.numStates = obj.numX + obj.numY;
            obj.numParams = 18;
            obj.initialdxdxo = eye(obj.numStates);
    
            obj.disturbanceSignal = @DroneControllers.dummyD;
        end

        function [t_span, stateProgression, model, stable, SC, index] = trajectory(model, SC, state, shouldPlot, ref, dist)
            % Initialization
            stateProgression = zeros(SC.t_final/SC.T, model.numStates);
            t_span = zeros(1, SC.t_final/SC.T);
            stateProgression(1, :) = state;
            stable = WithinLimits(state);
            % Drone Simulation
            index = 2;
            time = SC.T;
            while(stable && time <= SC.t_final)
                state_prev = stateProgression(index-1, :);
                [stable1, nextState1] = simulateTimestep(state_prev, SC.T, time, model, ref, dist);
                stateProgression(index, 1:19) = nextState1(1:19);
                stateProgression(index, 20) = 1/SC.T * (theta_ref(model, nextState1, ref) - theta_ref(model, state_prev, ref));
                stateProgression(index, 21) = 1/SC.T * (phi_ref(model, nextState1, ref) - phi_ref(model, state_prev, ref));
                stateProgression(index, 22) = 1/SC.T * (stateProgression(index, 20) - state_prev(20));
                stateProgression(index, 23) = 1/SC.T * (stateProgression(index, 21) - state_prev(21));
                t_span(index) = time;
                index = index + 1;
                time = time + SC.T;
                
                state_prev = stateProgression(index-1, :);
                [stable2, nextState2] = simulateTimestep(state_prev, SC.T, time, model, ref, dist);
                stateProgression(index, 1:19) = nextState2(1:19);
                stateProgression(index, 20) = 1/SC.T * (theta_ref(model, nextState2, ref) - theta_ref(model, nextState1, ref));
                stateProgression(index, 21) = 1/SC.T * (phi_ref(model, nextState2, ref) - phi_ref(model, nextState1, ref));
                stateProgression(index, 22) = 1/SC.T * (stateProgression(index, 20) - state_prev(20));
                stateProgression(index, 23) = 1/SC.T * (stateProgression(index, 21) - state_prev(21));
                t_span(index) = time;
                index = index + 1;
                time = time + SC.T;
                
                [stableDouble, nextStateDouble] = simulateTimestep(state_prev, SC.T*2, time, model, ref, dist);

                if (norm(nextState2(1:19) - nextStateDouble(1:19), inf) < 1e-7 && stableDouble)
                    SC.T = SC.T*2;   
                elseif (norm(nextState2(1:19) - nextStateDouble(1:19), inf) > 1e-4 || ~stable1 || ~stable2)
                    if SC.T > 1e-3
                        SC.T = SC.T/2;
                    end 
                end  

                stable = stable1 && stable2;
                
            end

            if stable
                final = stateProgression(index-1, :);
                prev = stateProgression(index-2, :);
                prevprev = stateProgression(index-3, :);
                stable = stable & CheckConvergence(final, prev, prevprev, ref);
            end

            % Plot
            if (shouldPlot)
                PlotDrone(stateProgression);
            end

            index = index-1;
        end 

        function p = extractP(drone)
            p = zeros(drone.numParams, 1);
            p(1) = drone.c1;
            p(2) = drone.c2;
            p(3) = drone.c3;
            p(4) = drone.c4;
            p(5) = drone.c5;
            p(6) = drone.c6;
            p(7) = drone.c7;
            p(8) = drone.c8;
            p(9) = drone.c9;
            p(10) = drone.c10;
            p(11) = drone.c11;
            p(12) = drone.c12;
        
            p(13) = drone.lambda1;
            p(14) = drone.lambda2;
            p(15) = drone.lambda3;
            p(16) = drone.lambda4;
            p(17) = drone.lambda5;
            p(18) = drone.lambda6;
        end
        
        function drone = assignP(drone, p)
            drone.c1      = p(1);
            drone.c2      = p(2);
            drone.c3      = p(3);
            drone.c4      = p(4);
            drone.c5      = p(5);
            drone.c6      = p(6);
            drone.c7      = p(7);
            drone.c8      = p(8);
            drone.c9      = p(9);
            drone.c10     = p(10);
            drone.c11     = p(11);
            drone.c12     = p(12);
            
            drone.lambda1 = p(13);
            drone.lambda2 = p(14);
            drone.lambda3 = p(15);
            drone.lambda4 = p(16);
            drone.lambda5 = p(17);
            drone.lambda6 = p(18);
        end
    end
end

function [stable, state_guess] = simulateTimestep(state_prev, timestep, time, model, ref, dist)
    state_guess = zeros(1, model.numStates);
    x3 = 100; % Drone starts 100 m in the air
    state_guess(3) = x3; 
    state_guess(19) = U_z(model, state_guess, ref);
    state_guess(20:23) = state_prev(20:23);
    count = 0;
    tol = 1e-12; 
    stable = 1;
    while abs(norm(F(state_prev, state_guess, model, timestep, ref, dist, time))) > tol  & count < 20
        if det(DF(state_guess, model, timestep, ref)) ~= 0
            state_guess(1:19) = state_guess(1:19) - transpose(DF(state_guess, model, timestep, ref)\transpose(F(state_prev, state_guess, model, timestep, ref, dist, time)));
            count = count + 1;
        else 
            stable = 0;
        end
    end
    stable = stable & count < 20 ; % & WithinLimits(state_guess);
end

function F = F(state_prev, state_guess, model, timestep, ref, dist, time)
    F = state_prev(1:model.numX) - state_guess(1:model.numX) + timestep/2*(f(state_prev, model, ref, dist, time-timestep) + f(state_guess, model, ref, dist, time));
end

function DF = DF(state_guess, model, timestep, ref)
    DF = -eye(model.numX);
    dfdx_plus = derivativesDroneIBHybrid.dfdx(model, state_guess, ref);
    DF = DF + timestep/2*dfdx_plus;
end
        
function f = f(state, model, ref, dist, time)

    x1 = state(1);
    x2 = state(2);
    x3 = state(3);
    x4 = state(4);
    x5 = state(5);
    x6 = state(6);
    x7 = state(7);
    x8 = state(8);
    x9 = state(9);
    x10 = state(10);
    x11 = state(11);
    x12 = state(12);

    [U_theta, U_phi, U_psi, U_z, U_zdot, theta_ref, phi_ref] = ...
        IntegralBackSteppingController(state, ref, model);

    f1 = x4; 
    f2 = x5;
    f3 = x6;
    
    f4 = 4*(U_z + dist.Fx(time))/model.m*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9)) ;
    f5 = 4*(U_z + dist.Fy(time))/model.m*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9));
    f6 = 4*(U_z + dist.Fz(time))/model.m*(cos(x7)*cos(x8)) - model.g;

    % f4 = 4*(U_z)/model.m*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9)) + dist.Fx(index);
    % f5 = 4*(U_z)/model.m*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9)) + dist.Fy(index);
    % f6 = 4*(U_z)/model.m*(cos(x7)*cos(x8)) - model.g + dist.Fz(index);
    
    f7 = x10; 
    f8 = x11;
    f9 = x12;
    
    f10 = 1/(model.Jtheta) * model.l * (2*U_theta) + dist.Tx(time); 
    f11 = 1/(model.Jphi) * model.l * (2*U_phi) + dist.Ty(time);
    f12 = 1/(model.Jpsi) * model.K_psi * (4*U_psi) + dist.Tz(time);
    
    f13 = theta_ref - x7;
    f14 = phi_ref - x8;
    f15 = ref.psi - x9;
    f16 = ref.z - x3;
    f17 = ref.x - x1;
    f18 = ref.y - x2;

    f19 = U_zdot;

    f = [f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 f13 f14 f15 f16 f17 f18 f19];
end  

function U_z = U_z(model, state, ref)
    x3 = state(3);
    x6 = state(6);
    x7 = state(7);
    x8 = state(8);
    x16 = state(16);

    c7 = model.c7;
    c8 = model.c8;
    lambda4 = model.lambda4;
    m = model.m; 
    g = model.g; 
    
    e7 = ref.z - x3;
    e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
    U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
end 

function theta_ref = theta_ref(model, state, ref)
    refx = ref.x;
    refxdot = ref.x_dot;
    refxddot = ref.x_ddot;

    m = model.m;
    c9 = model.c9;
    c10 = model.c10;
    lambda5 = model.lambda5; 

    x1 = state(1);
    x4 = state(4);
    x17 = state(17);
    x19 = state(19);

    e9 = refx - x1;
    e10 = c9*e9 + refxdot + lambda5*x17 - x4;
    U_z = x19;

    theta_ref = m/U_z*((1-c9^2  + lambda5)*e9 + (c9 +c10)*e10 + refxddot - c9*lambda5*x17);
end

function phi_ref = phi_ref(model, state, ref)
    refy = ref.y;
    refydot = ref.y_dot;
    refyddot = ref.y_ddot;

    m = model.m;
    c11 = model.c11;
    c12 = model.c12;
    lambda6 = model.lambda6; 

    x2 = state(2);
    x5 = state(5);
    x18 = state(18);
    x19 = state(19);

    e11 = refy - x2;
    e12 = c11*e11 + refydot + lambda6*x18 - x5;
    U_z = x19;

    phi_ref =   m/U_z*( -(1-c11^2 + lambda6)*e11 - (c11+c12)*e12 - refyddot + c11*lambda6*x18);

end

function [U_theta, U_phi, U_psi, U_z, U_zdot, theta_reference, phi_reference] = IntegralBackSteppingController(state, ref, model)

    refz = ref.z;
    refzdot = ref.z_dot;
    refzddot = ref.z_ddot;
    refpsi = ref.psi;
    refpsidot = ref.psi_dot;
    refpsiddot = ref.psi_ddot;

    x3 = state(3);
    x6 = state(6);
    x7 = state(7);
    x8 = state(8);
    x9 = state(9);
    x10 = state(10);
    x11 = state(11);
    x12 = state(12);
    x13 = state(13);
    x14 = state(14);
    x15 = state(15);
    x16 = state(16);
    x19 = state(19);
    x20 = state(20);
    x21 = state(21);
    x22 = state(22);
    x23 = state(23);

    c1 = model.c1;
    c2 = model.c2;
    c3 = model.c3;
    c4 = model.c4;
    c5 = model.c5;
    c6 = model.c6;
    c7 = model.c7; 
    c8 = model.c8;
    lambda1 = model.lambda1;
    lambda2 = model.lambda2;
    lambda3 = model.lambda3;
    lambda4 = model.lambda4;
    m = model.m;
    g = model.g;
    l = model.l;
    K_psi = model.K_psi;
    Jtheta = model.Jtheta;
    Jphi = model.Jphi;
    Jpsi = model.Jpsi;
    
    e7 = refz - x3;
    e8 = c7*e7 + refzdot + lambda4*x16 - x6;
    U_z = x19;

    theta_reference = theta_ref(model, state, ref);
    phi_reference =   phi_ref(model, state, ref);
    
    theta_ref_dot = x20;
    phi_ref_dot = x21;
    theta_ref_ddot = x22;
    phi_ref_ddot = x23;

    e1 = theta_reference - x7;
    e2 = c1*e1 + theta_ref_dot + lambda1*x13 - x10;
    e3 = phi_reference - x8;
    e4 = c3*e3 + phi_ref_dot + lambda2*x14 - x11;
    e5 = refpsi - x9;
    e6 = c5*e5 + refpsidot + lambda3*x15 - x12;

    z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
    e8_dot = c7*(refzdot - x6)+refzddot+lambda4*e7-z_ddot;
    ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
    ez_dot = ((1-c7^2+lambda4)*(refzdot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7);
    U_zdot = m/4*( cos(x7)*cos(x8)*ez_dot + ez*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11) )/(cos(x7)*cos(x8))^2;

    U_theta = Jtheta/(2*l)*((1-c1^2+lambda1)*e1 + (c1+c2)*e2 + theta_ref_ddot - c1*lambda1*x13);
    U_phi = Jphi/(2*l)*((1-c3^2+lambda2)*e3 + (c3+c4)*e4 + phi_ref_ddot - c3*lambda2*x14);
    U_psi = Jpsi/(4*K_psi)*((1-c5^2+lambda3)*e5 + (c5+c6)*e6 + refpsiddot - c5*lambda3*x15);
end

function [] = PlotDrone(stateProgression) 
    figure(100)
    sim = scatter3(0, 0, 0, 20, "blue", "filled");
    h_annot = annotation('textbox',[.9 .5 .1 .2],'String',0,'EdgeColor','none');
    title("Drone Sim");
    xlabel("x");
    ylabel("y");
    zlabel("z");
    zlim([0 15])

    for i = 1 : length(t_span)
        t = t_span(i);
        % update plot 
        sim.XData = stateProgression(i, 1);
        sim.YData = stateProgression(i, 2);
        sim.ZData = stateProgression(i, 3);
        set(h_annot,'String',t)
        drawnow
        pause(0);
    end
end

function converging = CheckConvergence(final, prev, prevprev, ref)

    finalTol = 1e-1;
    convergingTol = 1e-3;
    speedTol = 5e-2;
    numericalTol = 1e-8;

    converging = abs(final(1) - ref.x) <= finalTol;
    converging = converging & abs(final(2) - ref.y) <= finalTol;
    converging = converging & abs(final(3) - ref.z) <= finalTol;

    % converging = converging & abs(final(1) - ref.x) <= abs(prev(1) - ref.x) + numericalTol;
    % converging = converging & abs(final(2) - ref.y) <= abs(prev(2) - ref.y) + numericalTol;
    % converging = converging & abs(final(3) - ref.z) <= abs(prev(3) - ref.z) + numericalTol;

    converging = converging & abs(final(4)) <= speedTol;
    converging = converging & abs(final(5)) <= speedTol;
    converging = converging & abs(final(6)) <= speedTol;

    converging = converging & abs(final(1) - prev(1)) <= convergingTol;
    converging = converging & abs(final(2) - prev(2)) <= convergingTol;
    converging = converging & abs(final(3) - prev(3)) <= convergingTol;
    converging = converging & abs(final(1) - prev(1)) <= abs(prev(1)-prevprev(1))  + numericalTol;
    converging = converging & abs(final(2) - prev(2)) <= abs(prev(2)-prevprev(2))  + numericalTol;
    converging = converging & abs(final(3) - prev(3)) <= abs(prev(3)-prevprev(3))  + numericalTol;
end

function stable = WithinLimits(state)
    stable = true;
    % translational velocities
    stable = stable & (abs(state(4)) < 1e6);
    stable = stable & (abs(state(5)) < 1e6);
    stable = stable & (abs(state(6)) < 1e6);
    % rotational velocities
    stable = stable & (abs(state(10)) < 1e6);
    stable = stable & (abs(state(11)) < 1e6);
    stable = stable & (abs(state(12)) < 1e6);
    % Uz and discrete derivatives
    stable = stable & (abs(state(19)) < 1e6);
    stable = stable & (abs(state(20)) < 1e6);
    stable = stable & (abs(state(21)) < 1e6);
    stable = stable & (abs(state(22)) < 1e6);
    stable = stable & (abs(state(23)) < 1e6);
end

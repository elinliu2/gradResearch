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

        function [t_span, stateProgression, model, stable] = trajectory(model, SC, state, shouldPlot, ref)
            % Initialization
            t_span = 0:SC.T:SC.t_final; % amount of time we will simulate for; fixed step sizes
            tol = 1e-12; 
            stateProgression = zeros(length(t_span), model.numStates);
            stateProgression(1, :) = state;
            stable = WithinLimits(state);
            % Drone Simulation
            i = 2;
            while(stable & (i <= length(t_span)))
                state_guess = zeros(1, model.numStates);
                x3 = 100; % Drone starts 100 m in the air
                state_guess(3) = x3; 
                state_guess(19) = U_z(model, state_guess, ref);
                state_guess(20:23) = stateProgression(i-1, 20:23);

                state_prev = stateProgression(i-1, :);
                count = 0;
                while abs(norm(F(state_prev, state_guess, model, SC, ref))) > tol & count < 5000
                    state_guess(1:19) = state_guess(1:19) - transpose(DF(state_guess, model, SC, ref)\transpose(F(state_prev, state_guess, model, SC, ref)));
                    count = count + 1;
                end
                stable = stable & (count < 5000) & WithinLimits(state_guess);
                stateProgression(i, 1:19) = state_guess(1:19);
                stateProgression(i, 20) = 1/SC.T * (theta_ref(model, state_guess, ref) - theta_ref(model, state_prev, ref));
                stateProgression(i, 21) = 1/SC.T * (phi_ref(model, state_guess, ref) - phi_ref(model, state_prev, ref));
                stateProgression(i, 22) = 1/SC.T * (stateProgression(i, 20) - state_prev(20));
                stateProgression(i, 23) = 1/SC.T * (stateProgression(i, 21) - state_prev(21));
                i = i + 1;
            end

            % Plot
            if (shouldPlot)
                PlotDrone(stateProgression);
            end

        end 
    end
end

function F = F(state_prev, state_guess, model, SC, ref)
    F = state_prev(1:model.numX) - state_guess(1:model.numX) + SC.T/2*(f(state_prev, model, ref) + f(state_guess, model, ref));
end

function DF = DF(state_guess, model, SC, ref)
    DF = -eye(model.numX);
    dfdx_plus = derivativesDroneIBHybrid.dfdx(model, state_guess, ref);
    DF = DF + SC.T/2*dfdx_plus;
end
        
function f = f(state, model, ref)

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
    
    f4 = 4*U_z/model.m*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9));
    f5 = 4*U_z/model.m*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9));
    f6 = 4*U_z/model.m*(cos(x7)*cos(x8)) - model.g;
    
    f7 = x10; 
    f8 = x11;
    f9 = x12;
    
    f10 = 1/(model.Jtheta) * model.l * (2*U_theta); 
    f11 = 1/(model.Jphi) * model.l * (2*U_phi);
    f12 = 1/(model.Jpsi) * model.K_psi * (4*U_psi);
    
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

function stable = WithinLimits(state)
    stable = true;
    % translational velocities
    stable = stable & (abs(state(4)) < 10000);
    stable = stable & (abs(state(5)) < 10000);
    stable = stable & (abs(state(6)) < 10000);
    % rotational velocities
    stable = stable & (abs(state(10)) < 10000);
    stable = stable & (abs(state(11)) < 10000);
    stable = stable & (abs(state(12)) < 10000);
end

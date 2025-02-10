classdef DroneModelIB
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

        dfdx;
        dfdp;
        d2fdx2;
        d2fdxdp;
        nudgeParams;

        numStates;
        numParams;
        initialdxdxo;

        disturbanceSignal;
        initializeStates;
    end
    
    methods
        function obj = DroneModelIB(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6)
            obj.g = 9.8;
            obj.m = 1.4; % kg

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
            
            obj.dfdx = @derivatives_drone_pidmodel.dfdx;
            obj.dfdp = @derivatives_drone_pidmodel.dfdp; 
            obj.d2fdx2 = @derivatives_drone_pidmodel.d2fdx2;
            obj.d2fdxdp = @derivatives_drone_pidmodel.d2fdxdp;
            obj.nudgeParams = @derivatives_drone_pidmodel.nudgeParams;
    
            obj.numStates = 16;
            obj.numParams = 12;
            obj.initialdxdxo = eye(16);
    
            obj.disturbanceSignal = @DroneControllers.dummyD;
            obj.initializeStates = @InitializeStatesDrone;
        end

        function [time, stateProgression, model] = trajectory(model, SC, initialState, shouldPlot)
            t_span = 0:SC.T:SC.t_final; % amount of time we will simulate for; fixed step sizes
            options = odeset('RelTol',1e-12, 'abstol', 1e-12); % will be used to set options for ODE solver
        
            % initialize plot
            if (shouldPlot)
                figure(100)
                sim = scatter3(0, 0, 0, 20, "blue", "filled");
                h_annot = annotation('textbox',[.9 .5 .1 .2],'String',0,'EdgeColor','none');
                title("Drone Sim");
                xlabel("x");
                ylabel("y");
                zlabel("z");
                zlim([0 15])
            end
            
            [time, stateProgression] = ode45(@droneDynamics, t_span, initialState, options, model);
            global MODEL_COPY;
            model = MODEL_COPY;
            
            for i = 1 : length(time)
                t = time(i);
                if (shouldPlot)
                    % update plot 
                    sim.XData = stateProgression(i, 1);
                    sim.YData = stateProgression(i, 2);
                    sim.ZData = stateProgression(i, 3);
                    set(h_annot,'String',t)
                    drawnow
                    pause(0);
                end 
            end
        end

        function [dState] = droneDynamics(t, state, model)
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

            x = state(1);
            y = state(2);
            z = state(3); 

            dx1 = state(4); 
            dx2 = state(5);
            dx3 = state(6);

            theta = state(7); % theta (pitch)
            phi = state(8); % phi (roll)
            psi = state(9); % psi (yaw)

            dx7 = state(10);
            dx8 = state(11);
            dx9 = state(12);

            [U_theta, U_phi, U_psi, U_z, theta_ref, phi_ref] = IntegralBackSteppingController(state, ref, model);
            
            F1 = U_theta + U_psi + U_z;
            F2 = -U_theta + U_psi + U_z;
            F3 = U_phi - U_psi + U_z;
            F4 = -U_phi - U_psi + U_z;

            % sat_limit = 25;
            % F1 = min(max(F1, sat_limit) -sat_limit);
            % F2 = min(max(F2, sat_limit) -sat_limit);
            % F3 = min(max(F3, sat_limit) -sat_limit);
            % F4 = min(max(F4, sat_limit) -sat_limit);

            dx4 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*cos(psi) + sin(phi)*sin(psi));
            dx5 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*sin(psi) - sin(phi)*cos(psi));
            dx6 = 1/(model.m) * (F1+F2+F3+F4) * (cos(theta)*cos(phi)) - model.g;

            dx10 = 1/(model.Jtheta) * model.l * (2*U_theta); 
            dx11 = 1/(model.Jphi) * model.l * (2*U_phi);
            dx12 = 1/(model.Jpsi) * model.K_psi * (4*U_psi);
      
            dx13 = theta_ref - theta;
            dx14 = phi_ref - phi;
            dx15 = ref.psi - psi;
            dx16 = ref.z - z;
            dx17 = ref.x - x;
            dx18 = ref.y - y;

            dState = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10; dx11; dx12; dx13; dx14; dx15; dx16; dx17; dx18];

            % This allows the controller to update the parameter values
            global MODEL_COPY;
            MODEL_COPY = model;
        end

        function [U_theta, U_phi, U_psi, U_z, theta_ref, phi_ref] = IntegralBackSteppingController(state, ref, model)
            x = state(1);
            y = state(2);
            z = state(3); 

            x_dot = state(4);
            y_dot = state(5);
            z_dot = state(6);

            theta = state(7); % theta (pitch)
            phi = state(8); % phi (roll)
            psi = state(9); % psi (yaw)

            theta_dot = state(10);
            phi_dot = state(11);
            psi_dot = state(12);

            theta_integral_error = state(13);
            phi_integral_error = state(14);
            psi_integral_error = state(15);
            z_integral_error = state(16);
            x_integral_error = state(17);
            y_integral_error = state(18);

            e7 = ref.z - z;
            e8 = model.c7*e7 + ref.z_dot + model.lambda4*z_integral_error - z_dot;
            U_z = model.m/(cos(theta)*cos(phi)*4)*(model.g + (1-model.c7^2+model.lambda4)*e7 + (model.c7+model.c8)*e8 - model.c7*model.lambda4*z_integral_error);

            e9 = ref.x - x;
            e10 = model.c9*e9 + ref.x_dot + model.lambda5*x_integral_error - x_dot;
            theta_ref = model.m/U_z*((1-model.c9^2  + model.lambda5)*e9 + (model.c9 +model.c10)*e10 + ref.x_ddot - model.c9*model.lambda5*x_integral_error);

            e11 = ref.y - y;
            e12 = model.c11*e11 + ref.y_dot + model.lambda6*y_integral_error - y_dot;
            phi_ref =   model.m/U_z*( -(1-model.c11^2 + model.lambda6)*e11 - (model.c11+model.c12)*e12 - ref.y_ddot + model.c11*model.lambda6*y_integral_error);

            z_ddot = 4*U_z*cos(theta)*cos(phi)/model.m-model.g;
            e8_dot = model.c7*(ref.z_dot - z_dot)*ref.z_ddot+model.lambda4*e7-z_ddot;
            ez_dot = (1-model.c7^2+model.lambda4)*(ref.z_dot-z_dot)+(model.c7+model.c8)*e8_dot-model.c7*model.lambda4*e7;
            ez     = (1-model.c7^2+model.lambda4)*e7+(model.c7+model.c8)*e8-model.c7*model.lambda4*z_integral_error+model.g;
            U_zdot = model.m/4*(cos(theta)*cos(phi)*ez_dot - ez*(-sin(theta)*cos(phi)-cos(theta)*sin(phi)))/(cos(theta)*cos(phi))^2;
            
            x_ddot = 4*U_z*(sin(theta)*cos(phi)*cos(psi) + sin(phi)*sin(psi))/(model.m);
            e10_dot= model.c9*(ref.x_dot - x_dot) - ref.x_ddot + model.lambda5*e9 - x_ddot;
            ex_dot = (1-model.c9^2+model.lambda5)*(ref.x_dot-x_dot)+(model.c9+model.c10)*e10_dot+ref.x_dddot-model.c9*model.lambda5*e9;
            ex     = (1-model.c9^2+model.lambda5)*e9 + (model.c9+model.c10)*e10 + ref.x_ddot - model.c9*model.lambda5*x_integral_error;

            y_ddot = 4*U_z*(sin(theta)*cos(phi)*sin(psi) - sin(phi)*cos(psi));
            e12_dot= model.c11*(ref.y_dot - y_dot) - ref.y_ddot + model.lambda6*e11 - y_ddot;
            ey_dot = -(1-model.c11^2+model.lambda6)*(ref.y_dot-y_dot)-(model.c11+model.c12)*e12_dot-ref.y_dddot+model.c11*model.lambda6*e11;
            ey     = -(1-model.c11^2+model.lambda6)*e11 - (model.c11+model.c12)*e12 - ref.y_ddot + model.c11*model.lambda6*y_integral_error;
            
            theta_ref_dot = model.m * (U_z*ex_dot - ex*U_zdot)/U_z^2;
            phi_ref_dot = model.m * (U_z*ey_dot - ey*U_zdot)/U_z^2;

            x_dddot = 4/model.m*( ...
                U_zdot*(sin(phi)*sin(psi) + sin(theta)*cos(phi)*cos(psi)) + ...
                U_z*(cos(phi)*sin(psi)*phi_dot + sin(phi)*cos(psi)*psi_dot + ...
                     cos(theta)*theta_dot*cos(phi)*cos(psi) - ...
                     sin(theta)*sin(phi)*phi_dot*cos(psi) - ...
                     sin(theta)*cos(phi)*sin(psi)*psi_dot));

            e10_ddot = model.c9*(ref.x_ddot - x_ddot) - ref.x_dddot + model.lambda5*(ref.x_dot - x_dot) - x_dddot;
            ex_ddot = -(1-model.c9^2+model.lambda5)*(ref.x_ddot - x_ddot) - (model.c9+model.c10)*e10_ddot - ref.x_ddddot + model.c9*model.lambda5*(ref.x_dot - x_dot);
            
            z_dddot = 4/(model.m)*(U_zdot*cos(phi)*cos(theta)-U_z*sin(phi)*phi_dot*cos(theta)-U_z*cos(theta)*sin(theta)*theta_dot);
            e8_ddot = model.c7*(ref.z_ddot-z_ddot)+ref.z_dddot+model.lambda4*(ref.z_dot-z_dot)-z_dddot;
            ez_ddot = (1-model.c7^2+model.lambda4)*(ref.z_ddot-z_ddot)+(model.c7+model.c8)*e8_ddot-model.c7*model.lambda4*(ref.z_dot-z_dot);
            U_zdot_top_dot = -sin(theta)*cos(phi)*ez_dot - cos(theta)*sin(phi)*ez_dot + cos(theta)*cos(phi)*ez_ddot + ...
                cos(theta)*cos(phi)*ez - sin(theta)*sin(phi)*ez + sin(theta)*cos(phi)*ez_dot - ...
                sin(theta)*sin(phi)*ez + cos(theta)*cos(phi)*ez + cos(theta)*sin(phi)*ez_dot;
            U_zdot_top = cos(theta)*cos(phi)*ez_dot - ez*(-sin(theta)*cos(phi)-cos(theta)*sin(phi));
            U_zddot = model.m/4*(cos(theta)*cos(phi)^2*U_zdot_top_dot - U_zdot_top*2*cos(theta)*cos(phi)*(-sin(theta*cos(phi)-cos(theta)*sin(phi))))/(cos(theta)*cos(phi))^4;

            komh_dot = (U_z*ex_ddot)-(ex*U_zddot);
            komh = U_z*ex_dot - ex*U_zdot;
            theta_ref_ddot = model.m*(U_z^2*komh_dot-komh*2*U_z*U_zdot)/U_z^4;

            y_dddot = 4/model.m*(U_zdot*(-sin(phi)*cos(psi)+sin(theta)*cos(phi)*sin(psi)) ...
                    + U_z *(-cos(phi)*phi_dot*cos(psi)+sin(phi)*sin(psi)*psi_dot ...
                            +cos(theta)*theta_dot*cos(psi)*sin(psi) ...
                            -sin(theta)*sin(phi)*phi_dot*sin(psi) ...
                            +sin(theta)*cos(phi)*cos(psi)*psi_dot));
            e12_ddot = model.c11*(ref.y_ddot-y_ddot) - ref.y_dddot + model.lambda6*(ref.y_dot-y_dot)-y_dddot;
            ey_ddot = (1-model.c11^2+model.lambda6)*(ref.y_ddot-y_ddot)+(model.c11+model.c12)*e12_ddot+ref.y_ddddot-model.c11*model.lambda6*(ref.y_dot-y_dot);
            ciwyw_dot = (U_zdot*ey_dot+U_z*ey_ddot)-(U_zdot*ey_dot+ey*U_zddot);
            ciwyw = U_z*ey_dot - ey*U_zdot;
            phi_ref_ddot = model.m*(U_z^2*ciwyw_dot-ciwyw*2*U_z*U_zdot)/U_z^4;

            e1 = theta_ref - theta;
            e2 = model.c1*e1 + theta_ref_dot + model.lambda1*theta_integral_error - theta_dot;
            U_theta = model.Jtheta/(2*model.l)*((1-model.c1^2+model.lambda1)*e1 + (model.c1+model.c2)*e2 + theta_ref_ddot - model.c1*model.lambda1*theta_integral_error);

            e3 = phi_ref - phi;
            e4 = model.c3*e3 + phi_ref_dot + model.lambda2*phi_integral_error - phi_dot;
            U_phi = model.Jphi/(2*model.l)*((1-model.c3^2+model.lambda2)*e3 + (model.c3+model.c4)*e4 + phi_ref_ddot - model.c3*model.lambda2*phi_integral_error);

            e5 = ref.psi - psi;
            e6 = model.c5*e5 + ref.psi_dot + model.lambda3*psi_integral_error - psi_dot;
            U_psi = model.Jpsi/(4*model.K_psi)*((1-model.c5^2+model.lambda3)*e5 + (model.c5+model.c5)*e6 + ref.psi_ddot - model.c5*model.lambda3*psi_integral_error);
        end

        function initialState = InitializeStatesDrone(~, i, j)
            initialState = zeros(16, 1);
            initialState(3) = 100; % Drone starts 100 m in the air
            initialState(7) = i;
            initialState(8) = j;
        end
    end
end


classdef DroneModel
    properties
        g;
        m;

        kp1;
        ki1;
        kd1;

        kp2;
        ki2;
        kd2;

        kp3;
        ki3;
        kd3;

        kp4;
        ki4;
        kd4;

        l; % l is the distance between Ob and each propeller; 
        K_psi; %  is the thrust-to-moment gain ;

        Ja;
        Jb;
        Jc;

        disturbance;

        dfdx;
        dfdp;
        d2fdx2;
        d2fdxdp;
        nudgeParams;

        numStates;
        numParams;
        initialdxdxo;

        controller;
        disturbanceSignal;
        initializeStates;
    end
    
    methods
        function obj = DroneModel(kp1, ki1, kd1, kp2, ki2, kd2, kp3, ki3, kd3, kp4, ki4, kd4)
            obj.g = 9.8;
            obj.m = 1.4; % kg

            obj.l = 0.2; % m
            obj.K_psi = 4; % Nm

            obj.Ja = 0.03; % kgm^2
            obj.Jb = 0.03; % kgm^2
            obj.Jc = 0.04; % kgm^2

            obj.kp1 = kp1;
            obj.ki1 = ki1;
            obj.kd1 = kd1;

            obj.kp2 = kp2;
            obj.ki2 = ki2;
            obj.kd2 = kd2;

            obj.kp3 = kp3;
            obj.ki3 = ki3;
            obj.kd3 = kd3;

            obj.kp4 = kp4;
            obj.ki4 = ki4;
            obj.kd4 = kd4;

            obj.disturbance = 0;
            
            obj.dfdx = @derivatives_drone_pidmodel.dfdx;
            obj.dfdp = @derivatives_drone_pidmodel.dfdp; 
            obj.d2fdx2 = @derivatives_drone_pidmodel.d2fdx2;
            obj.d2fdxdp = @derivatives_drone_pidmodel.d2fdxdp;
            obj.nudgeParams = @derivatives_drone_pidmodel.nudgeParams;
    
            obj.numStates = 16;
            obj.numParams = 12;
            obj.initialdxdxo = eye(16);
    
            obj.controller = DroneControllers();
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
            pdot_x = state(4);
            pdot_y = state(5);
            pdot_z = state(6);

            theta = state(7); % theta (pitch)
            phi = state(8); % phi (roll)
            psi = state(9); % psi (yaw)

            dot_a = state(10);
            dot_b = state(11);
            dot_c = state(12);

            dx1 = pdot_x; 
            dx2 = pdot_y;
            dx3 = pdot_z;

            [F1, F2, F3, F4, model] = DroneControllers.dummyU(state, t, model);

            dx4 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*cos(psi) + sin(phi)*sin(psi));
            dx5 = 1/(model.m) * (F1+F2+F3+F4) * (sin(theta)*cos(phi)*sin(psi) - sin(phi)*cos(psi));
            dx6 = 1/(model.m) * (F1+F2+F3+F4) * (cos(theta)*cos(phi)) - model.g;

            dx7 = dot_a;
            dx8 = dot_b;
            dx9 = dot_c;

            dx10 = 1/(model.Ja) * model.l * (F1 - F2); 
            dx11 = 1/(model.Jb) * model.l * (F3 - F4);
            dx12 = 1/(model.Jc) * model.K_psi * (F1 + F2 - F3 - F4);

            ref(3) = 101;
            delta_z = ref(3) - state(3);

            % theta_ref = atan2(delta_x, 1);
            % phi_ref = atan2(-delta_y, 1);
            % psi_ref = 0; 

            theta_ref = pi/6;
            phi_ref = 0;
            psi_ref = t;
      
            dx13 = theta_ref - state(7);
            dx14 = phi_ref   - state(8);
            dx15 = psi_ref   - state(9);
            % dx16 = norm(ref-state(1:3));
            dx16 = delta_z;

            dState = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10; dx11; dx12; dx13; dx14; dx15; dx16];

            % This allows the controller to update the parameter values
            global MODEL_COPY;
            MODEL_COPY = model;
        end

        function initialState = InitializeStatesDrone(~, i, j)
            initialState = zeros(16, 1);
            initialState(3) = 100; % Drone starts 100 m in the air
            initialState(7) = i;
            initialState(8) = j;
        end
    end
end


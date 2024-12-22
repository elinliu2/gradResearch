classdef DroneModel
    properties
        g;
        m;

        kp;
        ki;
        kd;

        d_x1; % offset of propeller from CoM
        d_x2;
        d_x3;
        d_x4;

        d_y1;
        d_y2;
        d_y3;
        d_y4;

        d_tau1; % drag
        d_tau2;
        d_tau3;
        d_tau4;

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
    end
    
    methods
        function obj = DroneModel(kp, ki, kd)
            obj.g = 9.8;
            obj.m = 1; % 1 kg

            obj.d_x1 = 0.2;
            obj.d_x2 = 0.2;
            obj.d_x3 = 0.2;
            obj.d_x4 = 0.2;
    
            obj.d_y1 = 0;
            obj.d_y2 = 0;
            obj.d_y3 = 0;
            obj.d_y4 = 0;
    
            obj.d_tau1 = 0.002;
            obj.d_tau2 = 0.002;
            obj.d_tau3 = 0.002;
            obj.d_tau4 = 0.002;

            obj.Ja = 0.4;
            obj.Jb = 0.4;
            obj.Jc = 0.4;


            obj.kp = kp;
            obj.ki = ki;
            obj.kd = kd;

            obj.disturbance = 0;
            
            obj.dfdx = @derivatives_drone_pidmodel.dfdx;
            obj.dfdp = @derivatives_drone_pidmodel.dfdp;
            obj.d2fdx2 = @derivatives_drone_pidmodel.d2fdx2;
            obj.d2fdxdp = @derivatives_drone_pidmodel.d2fdxdp;
            obj.nudgeParams = @derivatives_drone_pidmodel.nudgeParams;
    
            obj.numStates = 12;
            obj.numParams = 3;
            obj.initialdxdxo = eye(12);
    
            obj.controller = @DroneControllers.ParameterizedPIDCtrl;
            obj.disturbanceSignal = @DroneControllers.dummyD;
        end

        function [time, stateProgression, model] = trajectory(model, SC, initialState, shouldPlot)

            t_span = 0:SC.T:SC.t_final; % amount of time we will simulate for; fixed step sizes
            options = odeset('RelTol',1e-12, 'abstol', 1e-12); % will be used to set options for ODE solver
        
            % initialize plot
            if (shouldPlot)
                figure(100)
                hold on
                sim = scatter3d(0, 0, 0, 20, [0.5 0 0.5], "filled");
                h_annot = annotation('textbox',[.9 .5 .1 .2],'String',0,'EdgeColor','none');
                title("Pendulum Sim");
                xlabel("x");
                ylabel("y");
                xlim([-6 6])
                ylim([-6 6])
            else
                sim = 0;
                h_annot = 0;
            end

            
            [time, stateProgression] = ode45(@pendulumDynamics, t_span, initialState, options, model);
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

            a = state(7); % alpha
            b = state(8); % beta
            c = state(9); % gamma

            w_a = state(10);
            w_b = state(11);
            w_c = state(12);

            dx1 = pdot_x + model.disturbanceSignal(model.disturbance, t);
            dx2 = pdot_y;
            dx3 = pdot_z;

            [F1, F2, F3, F4] = model.controller;

            dx4 = 1/(model.m) * (F1+F2+F3+F4) * (cos(a)*sin(b)*cos(c) + sin(a)*sin(c));
            dx5 = 1/(model.m) * (F1+F2+F3+F4) * (sin(a)*sin(b)*cos(c) - cos(a)*sin(c));
            dx6 = 1/(model.m) * (F1+F2+F3+F4) * (cos(b)*cos(c)) - model.g;

            dx7 = w_a;
            dx8 = w_b;
            dx9 = w_c;

            w = [w_a; w_b; w_c];
            Jw = [model.Ja*w_a; model.Jb*w_b; model.Jc*w_c];
            wJw = cross(w, Jw);
            dx10 = 1/(model.Ja) * (F1*model.d_y1 + F2*model.d_y2 + F3*model.d_y3 + F4*model.d_y4) - wJw(1) ; 
            dx11 = 1/(model.Jb) * (-F1*model.d_x1 - F2*model.d_x2 - F3*model.d_x3 - F4*model.d_x4) - wJw(1) ; 
            dx12 = 1/(model.Jc) * (F1*model.d_tau1 + F2*model.d_tau2 + F3*model.d_tau3 + F4*model.d_tau4) - wJw(1) ; 

            dState = [dx1; dx2; dx3; dx4; dx5; dx6; dx7; dx8; dx9; dx10; dx11; dx12];

            % This allows the controller to update the parameter values
            global MODEL_COPY;
            MODEL_COPY = model;
        end
        
        
    end
end


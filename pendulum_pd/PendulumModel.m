classdef PendulumModel
    properties
        g;
        l;
        m;
        D;

        u;
        disturbance;

        kp;
        kd;

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
        % Constructor to initialize the properties
        function obj = PendulumModel(kp, kd)
            obj.g = 9.8;
            obj.l = obj.g/2;
            obj.m = 2/obj.g;
            obj.D = 0.1;
    
            obj.u = 0;
            obj.disturbance = 0;
    
            obj.kp = kp;
            obj.kd = kd;
    
            obj.dfdx = @derivatives_pdctrlmodel.dfdx;
            obj.dfdp = @derivatives_pdctrlmodel.dfdp;
            obj.d2fdx2 = @derivatives_pdctrlmodel.d2fdx2;
            obj.d2fdxdp = @derivatives_pdctrlmodel.d2fdxdp;
            obj.nudgeParams = @derivatives_pdctrlmodel.nudgeParams;
    
            obj.numStates = 2;
            obj.numParams = 2;
            obj.initialdxdxo = eye(2);
    
            obj.controller = @PendulumControllers.ParameterizedPDCtrl;
            obj.disturbanceSignal = @PendulumControllers.dummyD;
        end

        function [time, stateProgression, model] = trajectory(model, SC, initialState, shouldPlot)

            t_span = 0:SC.T:SC.t_final; % amount of time we will simulate for; fixed step sizes
            options = odeset('RelTol',1e-12, 'abstol', 1e-12); % will be used to set options for ODE solver
        
            % initialize plot
            if (shouldPlot)
                figure(100)
                hold on
                sim = scatter(0, 0, 20, [0.5 0 0.5], "filled");
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
                    sim.XData = model.l*sin(stateProgression(i, 1));
                    sim.YData = -model.l*cos(stateProgression(i, 1));
                    set(h_annot,'String',t)
                    drawnow
                    pause(0);
                end 
            end
        end

        function [dState] = pendulumDynamics(t, state, model)
            x1 = state(1);
            x2 = state(2);
        
            dx1 = x2 + model.disturbanceSignal(model.disturbance, t);
        
            [uVal, model] = model.controller(x1, x2, model);

            dx2 = -model.g/model.l*sin(x1) - model.D/(model.m*model.l)*x2 + 1/(model.m*model.l)*(uVal);

            dState = [dx1;dx2];

            % This allows the controller to update the parameter values
            global MODEL_COPY;
            MODEL_COPY = model;
        end
    end
end


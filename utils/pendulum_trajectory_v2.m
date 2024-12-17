function [time, stateProgression, SC] = pendulum_trajectory_v2(SC, initialState, u, d, shouldPlot)
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
    global SC_copy;
    SC_copy = SC;
    [time, stateProgression] = ode45(@traj, t_span, initialState, options, u, d);
    
    SC = SC_copy;

    for i = 1 : length(time)
        t = time(i);
        if (shouldPlot)
            % update plot 
            sim.XData = SC.l*sin(stateProgression(i, 1));
            sim.YData = -SC.l*cos(stateProgression(i, 1));
            set(h_annot,'String',t)
            drawnow
            pause(0);
        end 
    end

end

function [dState] = traj(t, state, u, d)
    global SC_copy;
    SC = SC_copy;
    x1 = state(1);
    x2 = state(2);

    dx1 = x2 + d(SC, t);

    [uVal, SC] = u(x1, x2, SC);

    dx2 = -SC.g/SC.l*sin(x1) - SC.D/(SC.m*SC.l)*x2 + 1/(SC.m*SC.l)*(uVal);
    dState = [dx1;dx2];
end


    


    
 
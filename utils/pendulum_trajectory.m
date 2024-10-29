function [stateProgression, time, uPlot] = pendulum_trajectory(SC, x1_o, x2_o, u, d, plot)
    % setup arrays 
    t_final = SC.t_final;
    size = t_final/SC.T;
    x1 = zeros([1 size]);
    x2 = zeros([1 size]);
    time = zeros([1 size]);
    ei = zeros([1 size]);
    e = zeros([1 size]);
    ed = zeros([1 size]);
    
    uPlot = zeros([1 size]);
    
    % initial conditions
    x1(1) = x1_o;
    x2(1) = x2_o;
    tolerance = 1e-5;

    state.x1 = x1(1);
    state.x2 = x2(1);
    state.e = e(1);
    state.ei = ei(1);
    state.ed = ed(1);

    stateProgression =  repmat(state, 1, size );

    
    % initialize plot
    if (plot)
        figure(100)
        hold on
        sim = scatter(0, 0, 20, [0.5 0 0.5], "filled");
        title("Pendulum Sim");
        xlabel("x");
        ylabel("y");
        xlim([-6 6])
        ylim([-6 6])
    end 
    
    for i = 2:t_final/SC.T
        SC.u = u(x1(i-1), x2(i-1), e(i), (e(i)-e(i-1))/SC.T, sum(e), SC);

        x1_guess = 0;
        x2_guess = 0;
        while 1
            dhdxh_term = derivativeh(x1_guess, SC)\h(x1(i-1), x1_guess, x2(i-1), x2_guess, SC);
            x1_guess = x1_guess - dhdxh_term(1);
            x2_guess = x2_guess - dhdxh_term(2);
            
            habs = max(abs( h(x1(i-1), x1_guess, x2(i-1), x2_guess, SC)));
            if ( habs < tolerance)
                break;
            end
        end
        x1(i) = x1_guess + d(SC);
        x2(i) = x2_guess;
        
        e(i) = SC.ref - x1(i-1);
        ei(i) = ei(i-1) + e(i);
        ed(i) = (e(i)-e(i-1))/SC.T;

        state.x1 = x1(i);
        state.x2 = x2(i);
        state.e = e(i);
        state.ei = ei(i);
        state.ed = ed(i);
        stateProgression(i) = state;

        time(i) = i*SC.T; 
        uPlot(i) = SC.u;
        if (plot)
            % update plot 
            sim.XData = SC.l*sin(x1(i));
            sim.YData = -SC.l*cos(x1(i));
            drawnow
            pause(0);
        end 
    end
end

function dhdx = derivativeh(x1k1, SC)
    dhdx = [-1 SC.T/2; (-SC.T/2*SC.g/SC.l*cos(x1k1)) (-1-SC.T/2*SC.D/(SC.m*SC.l))];
end

function h = h(x1k, x1k1, x2k, x2k1, SC)
    h = [x1k - x1k1 + SC.T/2*x2k + SC.T/2*x2k1; 
         x2k - x2k1 - SC.T/2*SC.g/SC.l*sin(x1k) - SC.T/2*SC.D/(SC.m*SC.l)*x2k + SC.T/(2*SC.m*SC.l) * SC.u - SC.T/2*SC.g/SC.l*sin(x1k1) - SC.T/2*SC.D/(SC.m*SC.l)*x2k1 + SC.T/(2*SC.m*SC.l) * SC.u];
end

classdef DroneModelIB_syms
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
        function obj = DroneModelIB_syms(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, lambda1, lambda2, lambda3, lambda4, lambda5, lambda6, m, g)
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
            
            obj.dfdx = @derivatives_drone_IB_syms_dfdx.dfdx;
            obj.dfdp = @derivatives_drone_IB_syms_dfdp.dfdp; 
            obj.d2fdx2 = @derivatives_drone_IB_syms.d2fdx2;
            obj.d2fdxdp = @derivatives_drone_IB_syms.d2fdxdp;
            obj.nudgeParams = @derivatives_drone_IB_syms.nudgeParams;
    
            obj.numStates = 24;
            obj.numParams = 18;
            obj.initialdxdxo = eye(obj.numStates);
    
            obj.disturbanceSignal = @DroneControllers.dummyD;
            obj.initializeStates = @InitializeStatesDrone;
        end

        function [time, stateProgression, model] = trajectory(model, SC, initialState, shouldPlot, ref)
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
            
            [time, stateProgression] = ode45(@droneDynamics, t_span, initialState, options, model, ref);
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

        function [dState] = droneDynamics(t, state, model, ref)

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

            [U_theta, U_phi, U_psi, U_z, theta_ref, phi_ref, ...
                U_zdot, theta_ref_dot, phi_ref_dot, theta_ref_ddot, phi_ref_ddot, U_zddot] = ...
                IntegralBackSteppingController(state, ref, model);
            
            % F1 = U_theta + U_psi + U_z;
            % F2 = -U_theta + U_psi + U_z;
            % F3 = U_phi - U_psi + U_z;
            % F4 = -U_phi - U_psi + U_z;

            % sat_limit = 25;
            % F1 = min(max(F1, sat_limit) -sat_limit);
            % F2 = min(max(F2, sat_limit) -sat_limit);
            % F3 = min(max(F3, sat_limit) -sat_limit);
            % F4 = min(max(F4, sat_limit) -sat_limit);

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
            f20 = U_zddot;

            f21 = theta_ref_dot;
            f22 = phi_ref_dot;
            f23 = theta_ref_ddot;
            f24 = phi_ref_ddot;

            dState = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13; f14; f15; f16; f17; f18; f19; f20; f21; f22; f23; f24];

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

function [U_theta, U_phi, U_psi, U_z, theta_ref, phi_ref, U_zdot, theta_ref_dot, phi_ref_dot, theta_ref_ddot, phi_ref_ddot, U_zddot] = IntegralBackSteppingController(state, ref, model)
    refx = ref.x;
    refxdot = ref.x_dot;
    refxddot = ref.x_ddot;
    refxdddot = ref.x_dddot;
    refxddddot = ref.x_ddddot;
    refy = ref.y;
    refydot = ref.y_dot;
    refyddot = ref.y_ddot;
    refydddot = ref.y_dddot;
    refyddddot = ref.y_ddddot;
    refz = ref.z;
    refzdot = ref.z_dot;
    refzddot = ref.z_ddot;
    refzdddot = ref.z_dddot;
    refpsi = ref.psi;
    refpsidot = ref.psi_dot;
    refpsiddot = ref.psi_ddot;

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
    x13 = state(13);
    x14 = state(14);
    x15 = state(15);
    x16 = state(16);
    x17 = state(17);
    x18 = state(18);
    x19 = state(19);
    x20 = state(20);
    x21 = state(21);
    x22 = state(22);
    x23 = state(23);
    x24 = state(24);

    c1 = model.c1;
    c2 = model.c2;
    c3 = model.c3;
    c4 = model.c4;
    c5 = model.c5;
    c6 = model.c6;
    c7 = model.c7; 
    c8 = model.c8;
    c9 = model.c9;
    c10 = model.c10;
    c11 = model.c11;
    c12 = model.c12;
    lambda1 = model.lambda1;
    lambda2 = model.lambda2;
    lambda3 = model.lambda3;
    lambda4 = model.lambda4;
    lambda5 = model.lambda5;
    lambda6 = model.lambda6;
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

    e9 = refx - x1;
    e10 = c9*e9 + refxdot + lambda5*x17 - x4;
    theta_ref = x21;
    
    e11 = refy - x2;
    e12 = c11*e11 + refydot + lambda6*x18 - x5;
    phi_ref =   x22;
    
    z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
    e8_dot = c7*(refzdot - x6)+refzddot+lambda4*e7-z_ddot;
    ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
    ez_dot = ((1-c7^2+lambda4)*(refzdot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7);
    U_zdot = x20;

    x_ddot = 4*U_z*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(m);
    e10_dot= c9*(refxdot - x4) - refxddot + lambda5*e9 - x_ddot;
    ex_dot = (1-c9^2+lambda5)*(refxdot-x4)+(c9+c10)*e10_dot+refxdddot-c9*lambda5*e9;
    ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + refxddot - c9*lambda5*x17;
    
    y_ddot = 4*U_z*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9))/(m);
    e12_dot= c11*(refydot - x5) - refyddot + lambda6*e11 - y_ddot;
    ey_dot = -(1-c11^2+lambda6)*(refydot-x5)-(c11+c12)*e12_dot-refydddot+c11*lambda6*e11;
    ey     = -(1-c11^2+lambda6)*e11 - (c11+c12)*e12 - refyddot + c11*lambda6*x18;
    
    theta_ref_dot = x23;
    phi_ref_dot = x24;

    z_dddot = 4/m*(U_z*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11) + U_zdot*(cos(x7)*cos(x8)));
    e8_ddot = c7*(refzddot-z_ddot) + refzdddot + lambda4*(refzdot - x6) - z_dddot;
    ez_ddot = (1-c7^2+lambda4)*(refzddot - z_ddot) + (c7+c8)*e8_ddot - c7*lambda4*(refzdddot - z_dddot);

    e1 = theta_ref - x7;
    e2 = c1*e1 + theta_ref_dot + lambda1*x13 - x10;
    e3 = phi_ref - x8;
    e4 = c3*e3 + phi_ref_dot + lambda2*x14 - x11;
    e5 = refpsi - x9;
    e6 = c5*e5 + refpsidot + lambda3*x15 - x12;
    
    x_dddot = 4/m*( ...
        U_zdot*(sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9)) ...
        + U_z*(cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 ...
             + cos(x7)*x10*cos(x8)*cos(x9) ...
             - sin(x7)*sin(x8)*x11*cos(x9) ...
             - sin(x7)*cos(x8)*sin(x9)*x12));
    e10_ddot = c9*(refxddot - x_ddot) - refxdddot + lambda5*(refxdot - x4) - x_dddot;
    ex_ddot = (1-c9^2+lambda5)*(refxddot - x_ddot) + (c9+c10)*e10_ddot + refxddddot - c9*lambda5*(refxdot - x4);

    y_dddot = 4/m*( ...
        U_zdot*(-sin(x8)*cos(x9) + sin(x7)*cos(x8)*sin(x9)) ...
        + U_z *(-cos(x8)*cos(x9)*x11 + sin(x8)*sin(x9)*x12 ...
             + cos(x7)*x10*cos(x9)*sin(x9) ...
             - sin(x7)*sin(x8)*x11*sin(x9) ...
             + sin(x7)*cos(x8)*cos(x9)*x12));

    e12_ddot = c11*(refyddot-y_ddot) - refydddot + lambda6*(refydot-x5)-y_dddot;
    ey_ddot = -(1-c11^2+lambda6)*(refyddot-y_ddot) - (c11+c12)*e12_ddot - refyddddot + c11*lambda6*(refydot-x5);

    U_zdottop = m/4*(cos(x7)*cos(x8)*ez_dot + ez*(-sin(x7)*cos(x8)*x10 - cos(x7)*sin(x8)*x11));
    
    betty = m*( (ex/(U_z*cos(x7)*cos(x8))^2*( ...
      + sin(x7)*cos(x8)*x10*ez_dot ...
      + cos(x7)*sin(x8)*x11*ez_dot ...
      - cos(x7)*cos(x8)*ez_ddot ...
      + cos(x7)*cos(x8)*x10^2*ez - sin(x7)*sin(x8)*x10*x11*ez ...
      + cos(x7)*cos(x8)*x11^2*ez - sin(x7)*sin(x8)*x10*x11*ez ...
      + sin(x7)*cos(x8)*x10*ez_dot + cos(x7)*sin(x8)*x11*ez_dot )) ...
      + ex/U_z^2*(U_zdottop*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11)*2*cos(x7)*cos(x8))/(cos(x7)*cos(x8))^3 ...
      + ex_ddot/U_z - (U_z*ex_dot-ex*U_zdot)*2*U_zdot/U_z^3 );

    august = m*( (ey/(U_z*cos(x7)*cos(x8))^2*( ...
      + sin(x7)*cos(x8)*x10*ez_dot ...
      + cos(x7)*sin(x8)*x11*ez_dot ...
      - cos(x7)*cos(x8)*ez_ddot ...
      + cos(x7)*cos(x8)*x10^2*ez - sin(x7)*sin(x8)*x10*x11*ez ...
      + cos(x7)*cos(x8)*x11^2*ez - sin(x7)*sin(x8)*x10*x11*ez ...
      + sin(x7)*cos(x8)*x10*ez_dot + cos(x7)*sin(x8)*x11*ez_dot )) ...
      + ey/U_z^2*(U_zdottop*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11)*2*cos(x7)*cos(x8))/(cos(x7)*cos(x8))^3 ...
      + ey_ddot/U_z - (U_z*ey_dot-ey*U_zdot)*2*U_zdot/U_z^3 );

    ivy =   (1-c1^2+lambda1)*e1 + (c1+c2)*e2 - lambda1*x13;
    rwylm = (1-c3^2+lambda2)*e3 + (c3+c4)*e4 - lambda2*x14;

    peace = m*ey*ez/(U_z*cos(x7)*cos(x8))^2;
    hoax =  m*ex*ez/(U_z*cos(x7)*cos(x8))^2;

    waolom = sin(x7)*cos(x8)*ivy + cos(x7)*sin(x8)*rwylm;
    bdilh  = 1/(1-hoax*sin(x7)*cos(x8));

    phi_ref_ddot = (peace*sin(x7)*cos(x8)*bdilh*(hoax*waolom+betty) + peace*waolom + august)...
    /(1-peace*cos(x7)*sin(x8)-peace*sin(x7)*cos(x8)*bdilh*hoax*cos(x7)*sin(x8));

    theta_ref_ddot = (hoax*(waolom + cos(x7)*sin(x8)*phi_ref_ddot) + betty) ...
        /(1-hoax*sin(x7)*cos(x8)); 

    U_theta = Jtheta/(2*l)*((1-c1^2+lambda1)*e1 + (c1+c2)*e2 + theta_ref_ddot - c1*lambda1*x13);
    U_phi = Jphi/(2*l)*((1-c3^2+lambda2)*e3 + (c3+c4)*e4 + phi_ref_ddot - c3*lambda2*x14);
    U_psi = Jpsi/(4*K_psi)*((1-c5^2+lambda3)*e5 + (c5+c6)*e6 + refpsiddot - c5*lambda3*x15);

    theta_ddot = ((1-c1^2+lambda1)*e1 + (c1+c2)*e2 + theta_ref_ddot - c1*lambda1*x13);
    phi_ddot = ((1-c3^2+lambda2)*e3 + (c3+c4)*e4 + phi_ref_ddot - c3*lambda2*x14);

    U_zdottopdot = -sin(x7)*cos(x8)*x10*ez_dot - cos(x7)*sin(x8)*x11*ez_dot + cos(x7)*cos(x8)*ez_ddot ...
        - cos(x7)*cos(x8)*x10^2*ez + sin(x7)*sin(x8)*x10*x11*ez ...
        - ez_dot*sin(x7)*cos(x8)*x10 - ez*sin(x7)*cos(x8)*theta_ddot ...
        + sin(x7)*sin(x8)*x10*x11*ez - cos(x7)*cos(x8)*x11^2*ez ...
        - ez_dot*cos(x7)*sin(x8)*x11 - ez*cos(x7)*sin(x8)*phi_ddot;

    U_zddot = ((cos(x7)*cos(x8))^2*U_zdottopdot - U_zdottop*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11))/(cos(x7)*cos(x8))^4;

end


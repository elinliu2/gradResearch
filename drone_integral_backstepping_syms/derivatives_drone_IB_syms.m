classdef derivatives_drone_IB_syms
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(model.numStates, model.numStates);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = derivatives_drone_IB_syms.df4dx(model, state, ref);
            dfdx(5, :) = derivatives_drone_IB_syms.df5dx(model, state, ref);
            dfdx(6, :) = derivatives_drone_IB_syms.df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = derivatives_drone_IB_syms.df10dx(model, state, ref);
            dfdx(11, :) = derivatives_drone_IB_syms.df11dx(model, state, ref);
            dfdx(12, :) = derivatives_drone_IB_syms.df12dx(model, state, ref);

            dfdx(13, :) = derivatives_drone_IB_syms.df13dx(model, state, ref);
            dfdx(14, :) = derivatives_drone_IB_syms.df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;
            dfdx(19, :) = derivatives_drone_IB_syms.df19dx(model, state, ref);
        end

        function df4dx = df4dx(model, state, ~)
            df4dx = zeros(1, model.numStates);
            
            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);

            df4dx(7)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
            df4dx(8)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
            df4dx(9)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
            df4dx(19)=(4*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
        end

        function df5dx = df5dx(model, state, ~)
            df5dx = zeros(1, model.numStates);
            
            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);

            df5dx(7)=(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
            df5dx(8)=-(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
            df5dx(9)=(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
            df5dx(19)=-(4*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
        end

        function df6dx = df6dx(model, state, ~)
            df6dx = zeros(1, model.numStates);

            m = model.m;

            x7 = state(7);
            x8 = state(8);
            x19 = state(19);

            df6dx(7)=-(4*x19*cos(x8)*sin(x7))/m;
            df6dx(8)=-(4*x19*cos(x7)*sin(x8))/m;
            df6dx(19)=(4*cos(x7)*cos(x8))/m;
        end

        function df10dx = df10dx(model, state, ref)
            df10dx  = zeros(1, model.numStates);
            
            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            refxdddot = ref.x_dddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c1 = model.c1;
            c2 = model.c2;
            c7 = model.c7; 
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            lambda1 = model.lambda1;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;
            m = model.m;
            g = model.g;

            x1 = state(1);
            x3 = state(3);
            x4 = state(4);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x17 = state(17);
            x19 = state(19);

           df10dx(1)=(c1 + c2)*((m*(x19*(c9*lambda5 - lambda5*(c9 + c10)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19) - (m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           df10dx(3)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dx(4)=- (c1 + c2)*((m*(x19*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - (m*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c1*m*(c9 + c10))/x19) - (m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19;
           df10dx(6)=(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dx(7)=c1^2 - (c1 + (m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/m))/x19^2)*(c1 + c2) - lambda1 - 1;
           df10dx(8)=-(m*(c1 + c2)*((4*x19^2*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           df10dx(9)=-4*(c1 + c2)*(c9 + c10)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9));
           df10dx(10)=-(c1 + c2)*((m^2*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)) + 1);
           df10dx(11)=-(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           df10dx(13)=lambda1*(c1 + c2) - c1*lambda1;
           df10dx(16)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dx(17)=- (c1 + c2)*((c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19 - (m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) - (m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19;
           df10dx(19)=- (c1 + c2)*((2*m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1) - (4*x19*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^2 + (c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2) - (m*(lambda1 - c1^2 + 1)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
        end

        function df11dx = df11dx(model, state, ref)
            df11dx  = zeros(1, model.numStates);
            
            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;
            refydddot = ref.y_dddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c3 = model.c3;
            c4 = model.c4;
            c7 = model.c7; 
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            lambda2 = model.lambda2;
            lambda4 = model.lambda4;
            lambda6 = model.lambda6;
            m = model.m;
            g = model.g;

            x2 = state(2);
            x3 = state(3);
            x5 = state(5);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x18 = state(18);
            x19 = state(19);
            
            df11dx(2)=(m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 - (c3 + c4)*((m*(x19*(c11*lambda6 - lambda6*(c11 + c12)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19);
           df11dx(3)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dx(5)=(c3 + c4)*((m*(x19*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - (m*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c3*m*(c11 + c12))/x19) + (m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19;
           df11dx(6)=-(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dx(7)=(m*(c3 + c4)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + 4*x19^2*cos(x7)*cos(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           df11dx(8)=c3^2 - lambda2 - (c3 + c4)*(c3 - (m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2) - 1;
           df11dx(9)=4*m*(c3 + c4)*(c11 + c12)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7));
           df11dx(10)=(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           df11dx(11)=(c3 + c4)*((m^2*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2) - 1);
           df11dx(14)=lambda2*(c3 + c4) - c3*lambda2;
           df11dx(16)=-(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dx(18)=(c3 + c4)*((c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19 - (m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) + (m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19;
           df11dx(19)=(c3 + c4)*((2*m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18) + (refydot - x5)*(lambda6 - c11^2 + 1) + x19*(c11 + c12)*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) - c11*lambda6*(refy - x2)))/x19^2 + (c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2) + (m*(lambda2 - c3^2 + 1)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
        end

        function df12dx = df12dx(model, ~, ~)
            df12dx  = zeros(1, model.numStates);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

            df12dx(9)=c5^2 - lambda3 - c5*(c5 + c6) - 1;
            df12dx(12)=- c5 - c6;
            df12dx(15)=lambda3*(c5 + c6) - c5*lambda3;
        end

        function df13dx = df13dx(model, state, ref)
            df13dx = zeros(1, model.numStates);
            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;
            m = model.m;

            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

            df13dx(1)=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           df13dx(4)=-(m*(c9 + c10))/x19;
           df13dx(7)=-1;
           df13dx(17)=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           df13dx(19)=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
        end

        function df14dx = df14dx(model, state, ref)
            df14dx  = zeros(1, model.numStates);
            
            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;

            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;
            m = model.m;

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

            df14dx(2)=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           df14dx(5)=(m*(c11 + c12))/x19;
           df14dx(8)=-1;
           df14dx(18)=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           df14dx(19)=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
        end

        function df19dx = df19dx(model, state, ref)
            df19dx  = zeros(1, model.numStates);
            
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;
            
            m = model.m;
            g = model.g;
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

            df19dx(3)=-(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(6)=-(m*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(7)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
           df19dx(8)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
           df19dx(10)=(m*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8));
           df19dx(11)=(m*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)^2);
           df19dx(16)=-(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(19)=- c7 - c8;

        end

        function dfdp = dfdp(model, state, ~, ref)
            dfdp = zeros(model.numStates, model.numParams);
            
            dfdp(10, :) = derivatives_drone_IB_syms.df10dp(model, state, ref);
            dfdp(11, :) = derivatives_drone_IB_syms.df11dp(model, state, ref);
            dfdp(12, :) = derivatives_drone_IB_syms.df12dp(model, state, ref);

            dfdp(13, :) = derivatives_drone_IB_syms.df13dp(model, state, ref);
            dfdp(14, :) = derivatives_drone_IB_syms.df14dp(model, state, ref);
            dfdp(19, :) = derivatives_drone_IB_syms.df19dp(model, state, ref);
        end
        
        function df10dp = df10dp(model, state, ref)
            df10dp = zeros(1, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            refxdddot = ref.x_dddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c1 = model.c1;
            c2 = model.c2;
            c7 = model.c7; 
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            lambda1 = model.lambda1;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;
            m = model.m;
            g = model.g;

            x1 = state(1);
            x3 = state(3);
            x4 = state(4);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
            x13 = state(13);
            x16 = state(16);
            x17 = state(17);
            x19 = state(19);

           df10dp(1)=c1*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19) - (c1 + c2)*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19) - x10 + (m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^2;
           df10dp(2)=lambda1*x13 - x10 - c1*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19) + (m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^2;
           df10dp(7)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dp(8)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dp(9)=(m*(lambda1 - c1^2 + 1)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19 - (c1 + c2)*((m*(x19*(refxddot - (c9 + c10)*(refxdot - x4) + c9*(refxdot - x4) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c1*m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19);
           df10dp(10)=(m*(lambda1 - c1^2 + 1)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19 - (c1 + c2)*((m*(x19*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c1*m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19);
           df10dp(13)=x13*(c1 + c2) - c1*x13 - x7 + (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19;
           df10dp(16)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df10dp(17)=(c1 + c2)*((m*(x19*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c1*m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19) + (m*(lambda1 - c1^2 + 1)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19;
        end
        
        function df11dp = df11dp(model, state, ref)
            df11dp = zeros(1, model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;
            refydddot = ref.y_dddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c3 = model.c3;
            c4 = model.c4;
            c7 = model.c7; 
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            lambda2 = model.lambda2;
            lambda4 = model.lambda4;
            lambda6 = model.lambda6;
            m = model.m;
            g = model.g;

            x2 = state(2);
            x3 = state(3);
            x5 = state(5);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
            x14 = state(14);
            x16 = state(16);
            x18 = state(18);
            x19 = state(19);

            df11dp(3)=c3*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19) - (c3 + c4)*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19) - x11 - (m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^2;
           df11dp(4)=lambda2*x14 - x11 - c3*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19) - (m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^2;
           df11dp(7)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dp(8)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dp(11)=(c3 + c4)*((m*(x19*(refyddot - (c11 + c12)*(refydot - x5) + c11*(refydot - x5) - 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c3*m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19) - (m*(lambda2 - c3^2 + 1)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19;
           df11dp(12)=- (c3 + c4)*((m*(x19*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c3*m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19) - (m*(lambda2 - c3^2 + 1)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19;
           df11dp(14)=x14*(c3 + c4) - c3*x14 - x8 - (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19;
           df11dp(16)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           df11dp(18)=- (c3 + c4)*((m*(x19*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c3*m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19) - (m*(lambda2 - c3^2 + 1)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19;

        end
        
        function df12dp = df12dp(model, state, ref)
            df12dp = zeros(1, model.numParams);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;
            refpsi = ref.psi; 
            refpsidot = ref.psi_dot;

            x9 = state(9);
            x12 = state(12);
            x15 = state(15); 

            df12dp(5)=refpsidot - x12 + (c5 + c6)*(refpsi - x9) - c5*(refpsi - x9);
           df12dp(6)=refpsidot - x12 + lambda3*x15 + c5*(refpsi - x9);
           df12dp(15)=refpsi - x9 - c5*x15 + x15*(c5 + c6);
        end
        
        function df13dp = df13dp(model, state, ref)
            df13dp = zeros(1, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;
            m = model.m;

            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

            df13dp(9)=(m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19;
           df13dp(10)=(m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19;
           df13dp(17)=(m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19;
        end
        
        function df14dp = df14dp(model, state, ref)
            df14dp = zeros(1, model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;
            m = model.m;

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

            df14dp(11)=-(m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19;
           df14dp(12)=-(m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19;
           df14dp(18)=-(m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19;
        end

        function df19dp = df19dp(model, state, ref)
            df19dp = zeros(1, model.numParams);

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;
            
            m = model.m;
            g = model.g;
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);
            
            df19dp(7)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*cos(x7)^2*cos(x8)^2);
           df19dp(8)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(4*cos(x7)^2*cos(x8)^2);
           df19dp(16)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2);
        end

        function d2fdx2 = d2fdx2(model, state, ref)
            d2fdx2 = zeros(model.numStates, model.numStates, model.numStates);

            d2fdx2(4, :, :) = derivatives_drone_IB_syms.d2f4dx2(model, state, ref);
            d2fdx2(5, :, :) = derivatives_drone_IB_syms.d2f5dx2(model, state, ref);
            d2fdx2(6, :, :) = derivatives_drone_IB_syms.d2f6dx2(model, state, ref);

            d2fdx2(10, :, :) = derivatives_drone_IB_syms.d2f10dx2(model, state, ref);
            d2fdx2(11, :, :) = derivatives_drone_IB_syms.d2f11dx2(model, state, ref);
            d2fdx2(13, :, :) = derivatives_drone_IB_syms.d2f13dx2(model, state, ref);
            d2fdx2(14, :, :) = derivatives_drone_IB_syms.d2f14dx2(model, state, ref);
            d2fdx2(19, :, :) = derivatives_drone_IB_syms.d2f19dx2(model, state, ref);
        end

        function d2f4dx2 = d2f4dx2(model, state, ~)
            d2f4dx2 = zeros(model.numStates, model.numStates);

            
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);
            
            
            m = model.m;

            d2f4dx2(7,7)=-(4*x19*cos(x8)*cos(x9)*sin(x7))/m;
           d2f4dx2(7,8)=-(4*x19*cos(x7)*cos(x9)*sin(x8))/m;
           d2f4dx2(7,9)=-(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
           d2f4dx2(7,19)=(4*cos(x7)*cos(x8)*cos(x9))/m;
           d2f4dx2(8,7)=-(4*x19*cos(x7)*cos(x9)*sin(x8))/m;
           d2f4dx2(8,8)=-(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f4dx2(8,9)=(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f4dx2(8,19)=(4*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f4dx2(9,7)=-(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
           d2f4dx2(9,8)=(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f4dx2(9,9)=-(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f4dx2(9,19)=(4*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f4dx2(19,7)=(4*cos(x7)*cos(x8)*cos(x9))/m;
           d2f4dx2(19,8)=(4*cos(x8)*sin(x9) - 4*cos(x9)*sin(x7)*sin(x8))/m;
           d2f4dx2(19,9)=(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9))/m;
        end
                
        function d2f5dx2 = d2f5dx2(model, state, ~)
            d2f5dx2 = zeros(model.numStates, model.numStates);

           
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);
            
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;

            d2f5dx2(3,7)=(sin(x7)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(cos(x7)^2*cos(x8)) - sin(x9)*(lambda4 - c7^2 + c7*(c7 + c8) + 1);
            d2f5dx2(7,7)=-(4*x19*cos(x8)*sin(x7)*sin(x9))/m;
           d2f5dx2(7,8)=-(4*x19*cos(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(7,9)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
           d2f5dx2(7,19)=(4*cos(x7)*cos(x8)*sin(x9))/m;
           d2f5dx2(8,7)=-(4*x19*cos(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(8,8)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f5dx2(8,9)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f5dx2(8,19)=-(4*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
           d2f5dx2(9,7)=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
           d2f5dx2(9,8)=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
           d2f5dx2(9,9)=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
           d2f5dx2(9,19)=(4*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
           d2f5dx2(19,7)=(4*cos(x7)*cos(x8)*sin(x9))/m;
           d2f5dx2(19,8)=-(4*cos(x8)*cos(x9) + 4*sin(x7)*sin(x8)*sin(x9))/m;
           d2f5dx2(19,9)=(4*sin(x8)*sin(x9) + 4*cos(x8)*cos(x9)*sin(x7))/m;
        end

        function d2f6dx2 = d2f6dx2(model, state, ~)
            d2f6dx2 = zeros(model.numStates, model.numStates);

           
            x7 = state(7);
            x8 = state(8);
            x19 = state(19);

            m = model.m;

            d2f6dx2(7,7)=-(4*x19*cos(x7)*cos(x8))/m;
           d2f6dx2(7,8)=(4*x19*sin(x7)*sin(x8))/m;
           d2f6dx2(7,19)=-(4*cos(x8)*sin(x7))/m;
           d2f6dx2(8,7)=(4*x19*sin(x7)*sin(x8))/m;
           d2f6dx2(8,8)=-(4*x19*cos(x7)*cos(x8))/m;
           d2f6dx2(8,19)=-(4*cos(x7)*sin(x8))/m;
           d2f6dx2(19,7)=-(4*cos(x8)*sin(x7))/m;
           d2f6dx2(19,8)=-(4*cos(x7)*sin(x8))/m;
        end

        function d2f10dx2 = d2f10dx2(model, state, ref)
            d2f10dx2 = zeros(model.numStates, model.numStates);
            
            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            refxdddot = ref.x_dddot;
            
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            x1 = state(1);
            x3 = state(3);
            x4 = state(4);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x17 = state(17);
            x19 = state(19);

            c1 = model.c1;
            c2 = model.c2;
            
            c7 = model.c7; 
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            
            lambda1 = model.lambda1;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;
            m = model.m;
            g = model.g;
            
            d2f10dx2(1,3)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(1,6)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(1,7)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(2*cos(x7)^3*cos(x8)^2))*(c1 + c2))/x19^2;
           d2f10dx2(1,8)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(2*cos(x7)^2*cos(x8)^3))*(c1 + c2))/x19^2;
           d2f10dx2(1,10)=(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(1,11)=(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(1,16)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(1,19)=(m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2 - (c1 + c2)*((2*m*(x19*(c9*lambda5 - lambda5*(c9 + c10)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*((c7 + c8)*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - c9*lambda5 + lambda5*(c9 + c10)))/x19^2 - (c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2);
           d2f10dx2(3,1)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(3,4)=-(m^2*(c1 + c2)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(3,7)=(m^2*(c1 + c2)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x7)*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f10dx2(3,8)=(m^2*(c1 + c2)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x8)*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f10dx2(3,10)=(m^2*sin(x7)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(3,11)=(m^2*sin(x8)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(3,17)=-(m^2*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(3,19)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(4,3)=-(m^2*(c1 + c2)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(4,6)=-(m^2*(c1 + c2)*(c9 + c10)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(4,7)=(m*(c1 + c2)*((m*(c9 + c10)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dx2(4,8)=(m*(c1 + c2)*((m*(c9 + c10)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dx2(4,10)=(m^2*sin(x7)*(c1 + c2)*(c9 + c10)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(4,11)=(m^2*sin(x8)*(c1 + c2)*(c9 + c10)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(4,16)=-(m^2*(c1 + c2)*(c9 + c10)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(4,19)=(c1 + c2)*((2*m*(x19*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - (m*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(lambda5 + (c7 + c8)*(c9 + c10) - c9^2 + c9*(c9 + c10) + 1))/x19^2 + (c1*m*(c9 + c10))/x19^2) + (m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19^2;
           d2f10dx2(6,1)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(6,4)=-(m^2*(c1 + c2)*(c9 + c10)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(6,7)=(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x7)*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f10dx2(6,8)=(m^2*(c1 + c2)*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x8)*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f10dx2(6,10)=(m^2*sin(x7)*(c1 + c2)*(c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(6,11)=(m^2*sin(x8)*(c1 + c2)*(c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(6,17)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(c9*lambda5 - lambda5*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(6,19)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(7,1)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(2*cos(x7)^3*cos(x8)^2))*(c1 + c2))/x19^2;
           d2f10dx2(7,3)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2))*(c1 + c2))/x19^2;
           d2f10dx2(7,4)=(m*(c1 + c2)*((m*(c9 + c10)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dx2(7,6)=(m*(c1 + c2)*((m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dx2(7,7)=-(m*(c1 + c2)*((m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^2) - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x8)^2*sin(x7)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(cos(x7)^3*cos(x8)^2) + (3*m*sin(x7)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^4*cos(x8)^2) - (4*x19^2*cos(x8)*cos(x9)*sin(x7)*(c9 + c10))/m))/x19^2;
           d2f10dx2(7,8)=-(m*(c1 + c2)*((m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) - (4*x19^2*cos(x7)*cos(x9)*sin(x8)*(c9 + c10))/m + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(cos(x7)^3*cos(x8)^3)))/x19^2;
           d2f10dx2(7,9)=4*cos(x7)*cos(x8)*sin(x9)*(c1 + c2)*(c9 + c10);
           d2f10dx2(7,10)=-(m*(c1 + c2)*((m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)*cos(x8)) + (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8))))/x19^2;
           d2f10dx2(7,11)=-(m^2*sin(x7)*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(7,16)=(m*(c1 + c2)*((m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dx2(7,17)=(m*(c1 + c2)*((m*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dx2(7,19)=(c1 + c2)*((2*m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/m))/x19^3 - (8*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/x19);
           d2f10dx2(8,1)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(2*cos(x7)^2*cos(x8)^3))*(c1 + c2))/x19^2;
           d2f10dx2(8,3)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3))*(c1 + c2))/x19^2;
           d2f10dx2(8,4)=(m*(c1 + c2)*((m*(c9 + c10)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dx2(8,6)=(m*(c1 + c2)*((m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dx2(8,7)=-(m*(c1 + c2)*((m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) - (4*x19^2*cos(x7)*cos(x9)*sin(x8)*(c9 + c10))/m + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(cos(x7)^3*cos(x8)^3)))/x19^2;
           d2f10dx2(8,8)=-(m*(c1 + c2)*((m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^2) - (4*x19^2*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x7)^2*sin(x8)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(cos(x7)^2*cos(x8)^3) + (3*m*sin(x8)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^4)))/x19^2;
           d2f10dx2(8,9)=-4*(c1 + c2)*(c9 + c10)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9));
           d2f10dx2(8,10)=-(m^2*sin(x7)*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(8,11)=-(m*(c1 + c2)*((m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)*cos(x8)) + (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)*cos(x8)^3)))/x19^2;
           d2f10dx2(8,16)=(m*(c1 + c2)*((m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dx2(8,17)=(m*(c1 + c2)*((m*(c9*lambda5 - lambda5*(c9 + c10))*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dx2(8,19)=(2*m*(c1 + c2)*((4*x19^2*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^3 - (8*(c1 + c2)*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/x19;
           d2f10dx2(9,7)=cos(x7)*cos(x8)*sin(x9)*(c9 + c10)*(4*c1 + 4*c2);
           d2f10dx2(9,8)=-(c9 + c10)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9))*(4*c1 + 4*c2);
           d2f10dx2(9,9)=(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7))*(4*c1 + 4*c2);
           d2f10dx2(10,1)=(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,3)=(m^2*sin(x7)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,4)=(m^2*sin(x7)*(c1 + c2)*(c9 + c10)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,6)=(m^2*sin(x7)*(c1 + c2)*(c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,7)=-(c1 + c2)*((m^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)) + (m^2*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^3*cos(x8)));
           d2f10dx2(10,8)=-(m^2*sin(x7)*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(10,16)=(m^2*sin(x7)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,17)=(m^2*sin(x7)*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(10,19)=(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8));
           d2f10dx2(11,1)=(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,3)=(m^2*sin(x8)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,4)=(m^2*sin(x8)*(c1 + c2)*(c9 + c10)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,6)=(m^2*sin(x8)*(c1 + c2)*(c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,7)=-(m^2*sin(x7)*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(11,8)=- (m^2*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)) - (m^2*sin(x8)^2*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)*cos(x8)^3);
           d2f10dx2(11,16)=(m^2*sin(x8)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,17)=(m^2*sin(x8)*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(11,19)=(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)*cos(x8)^2);
           d2f10dx2(16,1)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(16,4)=-(m^2*(c1 + c2)*(c9 + c10)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(16,7)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x7)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f10dx2(16,8)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x8)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f10dx2(16,10)=(m^2*sin(x7)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(16,11)=(m^2*sin(x8)*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(16,17)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(c9*lambda5 - lambda5*(c9 + c10))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(16,19)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(17,3)=-(m^2*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(17,6)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(c9*lambda5 - lambda5*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(17,7)=(c1 + c2)*((m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x7)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^2*cos(x7)^3*cos(x8)^2));
           d2f10dx2(17,8)=(c1 + c2)*((m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x8)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^2*cos(x7)^2*cos(x8)^3));
           d2f10dx2(17,10)=(m^2*sin(x7)*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dx2(17,11)=(m^2*sin(x8)*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dx2(17,16)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(c9*lambda5 - lambda5*(c9 + c10))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dx2(17,19)=(m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19^2 - (c1 + c2)*((m*(c7 + c8)*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2 - (c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2 + (m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f10dx2(19,1)=(m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2 - (c1 + c2)*((2*m*(x19*(c9*lambda5 - lambda5*(c9 + c10)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*((c7 + c8)*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - c9*lambda5 + lambda5*(c9 + c10)))/x19^2 - (c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2);
           d2f10dx2(19,3)=-(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(19,4)=(c1 + c2)*((2*m*(x19*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - (m*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(lambda5 + (c7 + c8)*(c9 + c10) - c9^2 + c9*(c9 + c10) + 1))/x19^2 + (c1*m*(c9 + c10))/x19^2) + (m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19^2;
           d2f10dx2(19,6)=-(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(19,7)=(c1 + c2)*((2*m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/m))/x19^3 - (8*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/x19);
           d2f10dx2(19,8)=(c1 + c2)*((2*m*((4*x19^2*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^3 - (8*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/x19);
           d2f10dx2(19,10)=(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8));
           d2f10dx2(19,11)=(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)*cos(x8)^2);
           d2f10dx2(19,16)=-(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dx2(19,17)=(m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19^2 - (c1 + c2)*((m*(c7 + c8)*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2 - (c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2 + (m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f10dx2(19,19)=(c1 + c2)*((6*m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^4 - (4*m*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1) - (4*x19*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^3 - (8*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/x19^2 + (2*c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3) + (2*m*(lambda1 - c1^2 + 1)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3;
        end

        function d2f11dx2 = d2f11dx2(model, state, ref)
            d2f11dx2 = zeros(model.numStates, model.numStates);

            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;
            refydddot = ref.y_dddot;
            refyddddot = ref.y_ddddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            x2 = state(2);
            x3 = state(3);
            x5 = state(5);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
            x16 = state(16);
            x18 = state(18);
            x19 = state(19);

            
            c3 = model.c3;
            c4 = model.c4;
            c7 = model.c7; 
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            
            lambda2 = model.lambda2;
            
            lambda4 = model.lambda4;
            
            lambda6 = model.lambda6;
            m = model.m;
            g = model.g;
           
            d2f11dx2(2,3)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(2,6)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(2,7)=-(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(2*cos(x7)^3*cos(x8)^2))*(c3 + c4))/x19^2;
           d2f11dx2(2,8)=-(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(2*cos(x7)^2*cos(x8)^3))*(c3 + c4))/x19^2;
           d2f11dx2(2,10)=-(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(2,11)=-(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(2,16)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(2,19)=(c3 + c4)*((2*m*(x19*(c11*lambda6 - lambda6*(c11 + c12)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*((c7 + c8)*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - c11*lambda6 + lambda6*(c11 + c12)))/x19^2 - (c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2) - (m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f11dx2(3,2)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(3,5)=(m^2*(c3 + c4)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(3,7)=- (m^2*(c3 + c4)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x7)*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f11dx2(3,8)=- (m^2*(c3 + c4)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x8)*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f11dx2(3,10)=-(m^2*sin(x7)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(3,11)=-(m^2*sin(x8)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(3,18)=(m^2*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(3,19)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(5,3)=(m^2*(c3 + c4)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(5,6)=(m^2*(c3 + c4)*(c11 + c12)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(5,7)=-(m*(c3 + c4)*((m*(c11 + c12)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dx2(5,8)=-(m*(c3 + c4)*((m*(c11 + c12)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dx2(5,10)=-(m^2*sin(x7)*(c3 + c4)*(c11 + c12)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(5,11)=-(m^2*sin(x8)*(c3 + c4)*(c11 + c12)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(5,16)=(m^2*(c3 + c4)*(c11 + c12)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(5,19)=- (c3 + c4)*((2*m*(x19*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - (m*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(lambda6 + (c7 + c8)*(c11 + c12) - c11^2 + c11*(c11 + c12) + 1))/x19^2 + (c3*m*(c11 + c12))/x19^2) - (m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19^2;
           d2f11dx2(6,2)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(6,5)=(m^2*(c3 + c4)*(c11 + c12)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(6,7)=- (m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x7)*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f11dx2(6,8)=- (m^2*(c3 + c4)*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x8)*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f11dx2(6,10)=-(m^2*sin(x7)*(c3 + c4)*(c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(6,11)=-(m^2*sin(x8)*(c3 + c4)*(c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(6,18)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(c11*lambda6 - lambda6*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(6,19)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(7,2)=-(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(2*cos(x7)^3*cos(x8)^2))*(c3 + c4))/x19^2;
           d2f11dx2(7,3)=-(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2))*(c3 + c4))/x19^2;
           d2f11dx2(7,5)=-(m*(c3 + c4)*((m*(c11 + c12)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dx2(7,6)=-(m*(c3 + c4)*((m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dx2(7,7)=(m*(c3 + c4)*((m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^2) - 4*x19^2*cos(x8)*sin(x7)*sin(x9)*(c11 + c12) - (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x8)^2*sin(x7)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(cos(x7)^3*cos(x8)^2) + (3*m*sin(x7)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^4*cos(x8)^2)))/x19^2;
           d2f11dx2(7,8)=(m*(c3 + c4)*((m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*cos(x7)*sin(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2) + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(cos(x7)^3*cos(x8)^3)))/x19^2;
           d2f11dx2(7,9)=4*m*cos(x7)*cos(x8)*cos(x9)*(c3 + c4)*(c11 + c12);
           d2f11dx2(7,10)=(m*(c3 + c4)*((m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)*cos(x8)) + (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8))))/x19^2;
           d2f11dx2(7,11)=(m^2*sin(x7)*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(7,16)=-(m*(c3 + c4)*((m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dx2(7,18)=-(m*(c3 + c4)*((m*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dx2(7,19)=(8*m*cos(x7)*cos(x8)*sin(x9)*(c3 + c4)*(c11 + c12))/x19 - (2*m*(c3 + c4)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + 4*x19^2*cos(x7)*cos(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^3;
           d2f11dx2(8,2)=-(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(2*cos(x7)^2*cos(x8)^3))*(c3 + c4))/x19^2;
           d2f11dx2(8,3)=-(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3))*(c3 + c4))/x19^2;
           d2f11dx2(8,5)=-(m*(c3 + c4)*((m*(c11 + c12)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dx2(8,6)=-(m*(c3 + c4)*((m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dx2(8,7)=(m*(c3 + c4)*((m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*cos(x7)*sin(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2) + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(cos(x7)^3*cos(x8)^3)))/x19^2;
           d2f11dx2(8,8)=(m*(c3 + c4)*(4*x19^2*(c11 + c12)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^2) - (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x7)^2*sin(x8)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(cos(x7)^2*cos(x8)^3) + (3*m*sin(x8)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^4)))/x19^2;
           d2f11dx2(8,9)=4*m*(c3 + c4)*(c11 + c12)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8));
           d2f11dx2(8,10)=(m^2*sin(x7)*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(8,11)=(m*(c3 + c4)*((m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)*cos(x8)) + (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)*cos(x8)^3)))/x19^2;
           d2f11dx2(8,16)=-(m*(c3 + c4)*((m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dx2(8,18)=-(m*(c3 + c4)*((m*(c11*lambda6 - lambda6*(c11 + c12))*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dx2(8,19)=-(c3 + c4)*((2*m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^3 + (8*m*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/x19);
           d2f11dx2(9,7)=4*m*cos(x7)*cos(x8)*cos(x9)*(c3 + c4)*(c11 + c12);
           d2f11dx2(9,8)=4*m*(c3 + c4)*(c11 + c12)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8));
           d2f11dx2(9,9)=4*m*(c3 + c4)*(c11 + c12)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9));
           d2f11dx2(10,2)=-(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,3)=-(m^2*sin(x7)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,5)=-(m^2*sin(x7)*(c3 + c4)*(c11 + c12)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,6)=-(m^2*sin(x7)*(c3 + c4)*(c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,7)=(m^2*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)) + (m^2*sin(x7)^2*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^3*cos(x8));
           d2f11dx2(10,8)=(m^2*sin(x7)*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(10,16)=-(m^2*sin(x7)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,18)=-(m^2*sin(x7)*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(10,19)=-(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8));
           d2f11dx2(11,2)=-(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,3)=-(m^2*sin(x8)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,5)=-(m^2*sin(x8)*(c3 + c4)*(c11 + c12)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,6)=-(m^2*sin(x8)*(c3 + c4)*(c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,7)=(m^2*sin(x7)*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(11,8)=(c3 + c4)*((m^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)) + (m^2*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)*cos(x8)^3));
           d2f11dx2(11,16)=-(m^2*sin(x8)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,18)=-(m^2*sin(x8)*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(11,19)=-(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)*cos(x8)^2);
           d2f11dx2(16,2)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(16,5)=(m^2*(c3 + c4)*(c11 + c12)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(16,7)=- (m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x7)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^3*cos(x8)^2);
           d2f11dx2(16,8)=- (m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (m^2*sin(x8)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^2*cos(x7)^2*cos(x8)^3);
           d2f11dx2(16,10)=-(m^2*sin(x7)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(16,11)=-(m^2*sin(x8)*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(16,18)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(c11*lambda6 - lambda6*(c11 + c12))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(16,19)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(18,3)=(m^2*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(18,6)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(c11*lambda6 - lambda6*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(18,7)=-(c3 + c4)*((m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x7)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^2*cos(x7)^3*cos(x8)^2));
           d2f11dx2(18,8)=-(c3 + c4)*((m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*x19^2*cos(x7)^2*cos(x8)^2) + (m^2*sin(x8)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^2*cos(x7)^2*cos(x8)^3));
           d2f11dx2(18,10)=-(m^2*sin(x7)*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dx2(18,11)=-(m^2*sin(x8)*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dx2(18,16)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(c11*lambda6 - lambda6*(c11 + c12))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dx2(18,19)=(c3 + c4)*((m*(c7 + c8)*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2 - (c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2 + (m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^3*cos(x7)^2*cos(x8)^2)) - (m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19^2;
           d2f11dx2(19,2)=(c3 + c4)*((2*m*(x19*(c11*lambda6 - lambda6*(c11 + c12)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*((c7 + c8)*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - c11*lambda6 + lambda6*(c11 + c12)))/x19^2 - (c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2) - (m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f11dx2(19,3)=(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(19,5)=- (c3 + c4)*((2*m*(x19*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - (m*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(lambda6 + (c7 + c8)*(c11 + c12) - c11^2 + c11*(c11 + c12) + 1))/x19^2 + (c3*m*(c11 + c12))/x19^2) - (m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19^2;
           d2f11dx2(19,6)=(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(19,7)=-(c3 + c4)*((2*m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + 4*x19^2*cos(x7)*cos(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^3 - (8*m*cos(x7)*cos(x8)*sin(x9)*(c11 + c12))/x19);
           d2f11dx2(19,8)=-((m*(4*x19*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + x19*(c11 + c12)*(4*cos(x8)*cos(x9) + 4*sin(x7)*sin(x8)*sin(x9))))/x19^2 + (2*m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^3)*(c3 + c4);
           d2f11dx2(19,9)=(c3 + c4)*((m*(4*x19*(c11 + c12)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)) + x19*(c11 + c12)*(4*sin(x8)*sin(x9) + 4*cos(x8)*cos(x9)*sin(x7))))/x19^2 - (8*m*(c11 + c12)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/x19);
           d2f11dx2(19,10)=-(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8));
           d2f11dx2(19,11)=-(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)*cos(x8)^2);
           d2f11dx2(19,16)=(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dx2(19,18)=(c3 + c4)*((m*(c7 + c8)*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2 - (c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2 + (m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*x19^3*cos(x7)^2*cos(x8)^2)) - (m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19^2;
           d2f11dx2(19,19)=- (c3 + c4)*((6*m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^4 - (4*m*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18) + (refydot - x5)*(lambda6 - c11^2 + 1) + x19*(c11 + c12)*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) - c11*lambda6*(refy - x2)))/x19^3 + (2*c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3 + (2*m*(c11 + c12)*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)))/x19^2) - (2*m*(lambda2 - c3^2 + 1)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3;           
        end

        function d2f13dx2 = d2f13dx2(model, state, ref)
            d2f13dx2 = zeros(model.numStates, model.numStates);

            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            
            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;
            m = model.m;
            
            d2f13dx2(1,19)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f13dx2(4,19)=(m*(c9 + c10))/x19^2;
           d2f13dx2(17,19)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f13dx2(19,1)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f13dx2(19,4)=(m*(c9 + c10))/x19^2;
           d2f13dx2(19,17)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f13dx2(19,19)=(2*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3;

         end

        function d2f14dx2 = d2f14dx2(model, state, ref)
            d2f14dx2 = zeros(model.numStates, model.numStates);

            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;
           
            m = model.m;

            d2f14dx2(2,19)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f14dx2(5,19)=-(m*(c11 + c12))/x19^2;
           d2f14dx2(18,19)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f14dx2(19,2)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f14dx2(19,5)=-(m*(c11 + c12))/x19^2;
           d2f14dx2(19,18)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f14dx2(19,19)=-(2*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3;

        end

        function d2f19dx2 = d2f19dx2(model, state, ref)
            d2f19dx2 = zeros(model.numStates, model.numStates);

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;
            
            m = model.m;
            g = model.g;
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

            d2f19dx2(3,7)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(3,8)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(3,10)=-(m*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(3,11)=-(m*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(6,7)=- (m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(6,8)=- (m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(6,10)=-(m*sin(x7)*(c7 + c8))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(6,11)=-(m*sin(x8)*(c7 + c8))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(7,3)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,6)=- (m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) - cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,7)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^2) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x8)^2*sin(x7)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (3*m*sin(x7)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^4*cos(x8)^2) + (m*sin(x7)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,8)=(m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(2*cos(x7)^3*cos(x8)^2) + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(cos(x7)^3*cos(x8)^3);
           d2f19dx2(7,10)=(m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) + (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)^3*cos(x8));
           d2f19dx2(7,11)=(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(7,16)=- (m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(8,3)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,6)=- (m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) - cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,7)=(m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(2*cos(x7)^2*cos(x8)^3) - (m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(2*cos(x7)^3*cos(x8)^2) + (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(cos(x7)^3*cos(x8)^3);
           d2f19dx2(8,8)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^2) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m + (8*x19*cos(x7)^2*sin(x8)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) + (3*m*sin(x8)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^4) + (m*sin(x8)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,10)=(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(8,11)=(m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) + (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)*cos(x8)^3);
           d2f19dx2(8,16)=- (m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(10,3)=-(m*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(10,6)=-(m*sin(x7)*(c7 + c8))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(10,7)=(m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) + (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)^3*cos(x8));
           d2f19dx2(10,8)=(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(10,16)=-(m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(11,3)=-(m*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(11,6)=-(m*sin(x8)*(c7 + c8))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(11,7)=(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(11,8)=(m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) + (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)*cos(x8)^3);
           d2f19dx2(11,16)=-(m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(16,7)=- (m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(16,8)=- (m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(16,10)=-(m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(16,11)=-(m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);

        end
        
        function d2fdxdp = d2fdxdp(model, state, ref)
            d2fdxdp = zeros(model.numStates, model.numStates, model.numParams);

            
            d2fdxdp(10, :, :) = derivatives_drone_IB_syms.d2f10dxdp(model, state, ref);
            d2fdxdp(11, :, :) = derivatives_drone_IB_syms.d2f11dxdp(model, state, ref);
            d2fdxdp(12, :, :) = derivatives_drone_IB_syms.d2f12dxdp(model, state, ref);

            d2fdxdp(13, :, :) = derivatives_drone_IB_syms.d2f13dxdp(model, state, ref);
            d2fdxdp(14, :, :) = derivatives_drone_IB_syms.d2f14dxdp(model, state, ref);
            d2fdxdp(19, :, :) = derivatives_drone_IB_syms.d2f19dxdp(model, state, ref);
        end

        function d2f10dxdp = d2f10dxdp(model, state, ref)
            d2f10dxdp = zeros(model.numStates, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;
            refxdddot = ref.x_dddot;
            
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            x1 = state(1);
            x3 = state(3);
            x4 = state(4);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
            x16 = state(16);
            x17 = state(17);
            x19 = state(19);

            c1 = model.c1;
            c2 = model.c2;
            
            c7 = model.c7; 
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            
            lambda1 = model.lambda1;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;
            m = model.m;
            g = model.g;

            d2f10dxdp(1,1)=(m*(x19*(c9*lambda5 - lambda5*(c9 + c10)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19 + (c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,2)=(m*(x19*(c9*lambda5 - lambda5*(c9 + c10)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,7)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(1,8)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(1,9)=- (c1 + c2)*((c1*c10*m)/x19 - (c10*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) - (c10*m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(1,10)=- (c1 + c2)*((m*(lambda5*x19 - (c9*m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c1*c9*m)/x19) - (c9*m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(1,13)=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,16)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(1,17)=- (c1 + c2)*((c1*m)/x19 + (m*(c10*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) - (m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(3,1)=(m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,2)=(m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,7)=(c8*m^2*(c1 + c2)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,8)=(m^2*(c1 + c2)*(c7*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + lambda4*cos(x7)*cos(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,9)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,10)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,16)=(m^2*(c1 + c2)*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(3,17)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(4,1)=(c1*m*(c9 + c10))/x19 - (m*(x19*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - (m*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (m*(c1 + c2)*(c9 + c10))/x19;
           d2f10dxdp(4,2)=- (m*(x19*(lambda5 - c9^2 + c9*(c9 + c10) + 1) - (m*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c1*m*(c9 + c10))/x19;
           d2f10dxdp(4,7)=(m^2*(c1 + c2)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(4,8)=(m^2*(c1 + c2)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(4,9)=- (c1 + c2)*((c1*m)/x19 + (m*(c10*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) - (m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(4,10)=- (c1 + c2)*((c1*m)/x19 + (m*(c9*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) - (m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(4,13)=-(m*(c9 + c10))/x19;
           d2f10dxdp(4,16)=(m^2*(c1 + c2)*(c9 + c10)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(4,17)=-(m*(c1 + c2))/x19;
           d2f10dxdp(6,1)=(m^2*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,2)=(m^2*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,7)=(m^2*(c1 + c2)*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,8)=(m^2*(c1 + c2)*(c7*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,9)=(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,10)=(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(6,16)=(m^2*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8));
           d2f10dxdp(6,17)=(m^2*(c1 + c2)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(7,1)=- c2 - (m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/m))/x19^2;
           d2f10dxdp(7,2)=- c1 - (m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9)*(c9 + c10))/m))/x19^2;
           d2f10dxdp(7,7)=-(m*(c1 + c2)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dxdp(7,8)=-(m*(c1 + c2)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dxdp(7,9)=-(m*(c1 + c2)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9))/m + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dxdp(7,10)=-(m*(c1 + c2)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2) + (4*x19^2*cos(x7)*cos(x8)*cos(x9))/m + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dxdp(7,13)=-1;
           d2f10dxdp(7,16)=-(m*(c1 + c2)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x8)*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f10dxdp(7,17)=-(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(2*cos(x7)^3*cos(x8)^2))*(c1 + c2))/x19^2;
           d2f10dxdp(8,1)=-(m*((4*x19^2*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,2)=-(m*((4*x19^2*(c9 + c10)*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,7)=-(m*(c1 + c2)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,8)=-(m*(c1 + c2)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,9)=-(m*(c1 + c2)*((4*x19^2*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,10)=-(m*(c1 + c2)*((4*x19^2*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m + (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,16)=-(m*(c1 + c2)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f10dxdp(8,17)=-(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(2*cos(x7)^2*cos(x8)^3))*(c1 + c2))/x19^2;
           d2f10dxdp(9,1)=-4*(c9 + c10)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9));
           d2f10dxdp(9,2)=-4*(c9 + c10)*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9));
           d2f10dxdp(9,9)=-(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))*(4*c1 + 4*c2);
           d2f10dxdp(9,10)=-(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))*(4*c1 + 4*c2);
           d2f10dxdp(10,1)=- (m^2*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)) - 1;
           d2f10dxdp(10,2)=- (m^2*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)) - 1;
           d2f10dxdp(10,7)=-(m^2*sin(x7)*(c1 + c2)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(10,8)=-(m^2*sin(x7)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(10,9)=-(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(10,10)=-(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(10,16)=-(m^2*sin(x7)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(10,17)=-(m^2*sin(x7)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f10dxdp(11,1)=-(m^2*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,2)=-(m^2*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,7)=-(m^2*sin(x8)*(c1 + c2)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,8)=-(m^2*sin(x8)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,9)=-(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,10)=-(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,16)=-(m^2*sin(x8)*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(11,17)=-(m^2*sin(x8)*(c1 + c2)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f10dxdp(13,2)=lambda1;
           d2f10dxdp(13,13)=c2;
           d2f10dxdp(16,1)=(m^2*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,2)=(m^2*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,8)=-(lambda4*m^2*(c1 + c2)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,9)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,10)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,16)=-(c8*m^2*(c1 + c2)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(16,17)=(m^2*(c1 + c2)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(17,1)=(c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19 - (m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19 + (m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(17,2)=(m^2*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2) - (c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           d2f10dxdp(17,7)=(m^2*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(17,8)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(c9*lambda5 - lambda5*(c9 + c10)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(17,10)=(c1 + c2)*((c1*lambda5*m)/x19 - (lambda5*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) + (lambda5*m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(17,13)=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           d2f10dxdp(17,16)=(m^2*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(17,17)=(c1 + c2)*((c1*c10*m)/x19 - (c10*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) + (c10*m*(lambda1 - c1^2 + 1))/x19;
           d2f10dxdp(19,1)=(m*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1) - (4*x19*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^2 - (2*m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 + (c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,2)=(m*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (c7 + c8)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1) - (4*x19*(c9 + c10)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^2 - (2*m*(x19*(refxdddot - (c9 + c10)*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (refxdot - x4)*(lambda5 - c9^2 + 1) - c9*lambda5*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,7)=(c1 + c2)*((m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 + (m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f10dxdp(19,8)=(c1 + c2)*((m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 + (m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f10dxdp(19,9)=- (c1 + c2)*((m*(refxddot - (c9 + c10)*(refxdot - x4) + c9*(refxdot - x4) - (c7 + c8)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)) + (8*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^2 - (2*m*(x19*(refxddot - (c9 + c10)*(refxdot - x4) + c9*(refxdot - x4) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (c1*m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2) - (m*(lambda1 - c1^2 + 1)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2;
           d2f10dxdp(19,10)=(c1 + c2)*((2*m*(x19*(refxddot - c9*(refxdot - x4) - lambda5*(refx - x1) + (4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*(c9*(refxdot - x4) - refxddot + lambda5*(refx - x1) + (c7 + c8)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) - (8*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m))/x19^2 - (c1*m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2) - (m*(lambda1 - c1^2 + 1)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2;
           d2f10dxdp(19,13)=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,16)=(m^2*(c1 + c2)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f10dxdp(19,17)=- (c1 + c2)*((2*m*(x19*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refx - x1 - c9*x17 + x17*(c9 + c10)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refxdot - x4 + (c9 + c10)*(refx - x1) + (c7 + c8)*(refx - x1 - c9*x17 + x17*(c9 + c10)) - c9*(refx - x1)))/x19^2 + (c1*m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2) - (m*(lambda1 - c1^2 + 1)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2;
        end

        function d2f11dxdp = d2f11dxdp(model, state, ref)
            d2f11dxdp = zeros(model.numStates, model.numParams);

            
            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;
            refydddot = ref.y_dddot;
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;
            x2 = state(2);
            x3 = state(3);
            x5 = state(5);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x18 = state(18);
            x19 = state(19);
            
            c3 = model.c3;
            c4 = model.c4;
            c7 = model.c7; 
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            
            lambda2 = model.lambda2;
            
            lambda4 = model.lambda4;
            
            lambda6 = model.lambda6;
            m = model.m;
            g = model.g;
           
            d2f11dxdp(2,3)=(m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 - (m*(x19*(c11*lambda6 - lambda6*(c11 + c12)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           d2f11dxdp(2,4)=(c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 - (m*(x19*(c11*lambda6 - lambda6*(c11 + c12)) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*cos(x7)^2*cos(x8)^2)))/x19^2;
           d2f11dxdp(2,7)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(2,8)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(2,11)=(c3 + c4)*((c3*c12*m)/x19 - (c12*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) + (c12*m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(2,12)=(c3 + c4)*((m*(lambda6*x19 - (c11*m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c3*c11*m)/x19) + (c11*m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(2,14)=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           d2f11dxdp(2,16)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(2,18)=(c3 + c4)*((c3*m)/x19 + (m*(c12*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) + (m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(3,3)=-(m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,4)=-(m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,7)=-(c8*m^2*(c3 + c4)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,8)=-(m^2*(c3 + c4)*(c7*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + lambda4*cos(x7)*cos(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,11)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,12)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,16)=-(m^2*(c3 + c4)*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(3,18)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8)))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(5,3)=(m*(x19*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - (m*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 - (c3*m*(c11 + c12))/x19 + (m*(c3 + c4)*(c11 + c12))/x19;
           d2f11dxdp(5,4)=(m*(x19*(lambda6 - c11^2 + c11*(c11 + c12) + 1) - (m*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2 + (c3*m*(c11 + c12))/x19;
           d2f11dxdp(5,7)=-(m^2*(c3 + c4)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(5,8)=-(m^2*(c3 + c4)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(5,11)=(c3 + c4)*((c3*m)/x19 + (m*(c12*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) + (m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(5,12)=(c3 + c4)*((c3*m)/x19 + (m*(c11*x19 - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2)))/x19^2) + (m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(5,14)=(m*(c11 + c12))/x19;
           d2f11dxdp(5,16)=-(m^2*(c3 + c4)*(c11 + c12)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(5,18)=(m*(c3 + c4))/x19;
           d2f11dxdp(6,3)=-(m^2*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,4)=-(m^2*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,7)=-(m^2*(c3 + c4)*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,8)=-(m^2*(c3 + c4)*(c7*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,11)=-(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,12)=-(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(6,16)=-(m^2*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8));
           d2f11dxdp(6,18)=-(m^2*(c3 + c4)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(7,3)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + 4*x19^2*cos(x7)*cos(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,4)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + 4*x19^2*cos(x7)*cos(x8)*sin(x9)*(c11 + c12) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,7)=(m*(c3 + c4)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,8)=(m*(c3 + c4)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,11)=(m*(c3 + c4)*(4*x19^2*cos(x7)*cos(x8)*sin(x9) + (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,12)=(m*(c3 + c4)*(4*x19^2*cos(x7)*cos(x8)*sin(x9) + (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,16)=(m*(c3 + c4)*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x8)*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^3*cos(x8)^2)))/x19^2;
           d2f11dxdp(7,18)=(m*((m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(2*cos(x7)^3*cos(x8)^2))*(c3 + c4))/x19^2;
           d2f11dxdp(8,3)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2 - c4;
           d2f11dxdp(8,4)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(c11 + c12)*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2 - c3;
           d2f11dxdp(8,7)=(m*(c3 + c4)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dxdp(8,8)=(m*(c3 + c4)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dxdp(8,11)=(m*(c3 + c4)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dxdp(8,12)=(m*(c3 + c4)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2) - 4*x19^2*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dxdp(8,14)=-1;
           d2f11dxdp(8,16)=(m*(c3 + c4)*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*cos(x7)^2*cos(x8)^3)))/x19^2;
           d2f11dxdp(8,18)=(m*((m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(2*cos(x7)^2*cos(x8)^3))*(c3 + c4))/x19^2;
           d2f11dxdp(9,3)=4*m*(c11 + c12)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7));
           d2f11dxdp(9,4)=4*m*(c11 + c12)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7));
           d2f11dxdp(9,11)=4*m*(c3 + c4)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7));
           d2f11dxdp(9,12)=4*m*(c3 + c4)*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7));
           d2f11dxdp(10,3)=(m^2*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,4)=(m^2*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,7)=(m^2*sin(x7)*(c3 + c4)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,8)=(m^2*sin(x7)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,11)=(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,12)=(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,16)=(m^2*sin(x7)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(10,18)=(m^2*sin(x7)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8));
           d2f11dxdp(11,3)=(m^2*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2) - 1;
           d2f11dxdp(11,4)=(m^2*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2) - 1;
           d2f11dxdp(11,7)=(m^2*sin(x8)*(c3 + c4)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(11,8)=(m^2*sin(x8)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(11,11)=(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(11,12)=(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(11,16)=(m^2*sin(x8)*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(11,18)=(m^2*sin(x8)*(c3 + c4)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*x19^2*cos(x7)*cos(x8)^2);
           d2f11dxdp(14,4)=lambda2;
           d2f11dxdp(14,14)=c4;
           d2f11dxdp(16,3)=-(m^2*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,4)=-(m^2*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,8)=(lambda4*m^2*(c3 + c4)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,11)=-(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,12)=-(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,16)=(c8*m^2*(c3 + c4)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(16,18)=-(m^2*(c3 + c4)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,3)=(m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19 - (c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19 - (m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,4)=(c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19 - (m^2*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,7)=-(m^2*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,8)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(c11*lambda6 - lambda6*(c11 + c12)))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,12)=- (c3 + c4)*((c3*lambda6*m)/x19 - (lambda6*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) - (lambda6*m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(18,14)=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           d2f11dxdp(18,16)=-(m^2*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12))*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(18,18)=- (c3 + c4)*((c3*c12*m)/x19 - (c12*m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(4*x19^2*cos(x7)^2*cos(x8)^2)) - (c12*m*(lambda2 - c3^2 + 1))/x19;
           d2f11dxdp(19,3)=(2*m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18) + (refydot - x5)*(lambda6 - c11^2 + 1) + x19*(c11 + c12)*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) - c11*lambda6*(refy - x2)))/x19^2 + (m*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 - (c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,4)=(2*m*(x19*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (refydot - x5)*(lambda6 - c11^2 + 1) - c11*lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refydddot + (c11 + c12)*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) + (c7 + c8)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18) + (refydot - x5)*(lambda6 - c11^2 + 1) + x19*(c11 + c12)*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) - c11*lambda6*(refy - x2)))/x19^2 + (c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,7)=-(c3 + c4)*((m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 + (m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f11dxdp(19,8)=-(c3 + c4)*((m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 + (m^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2));
           d2f11dxdp(19,11)=(m*(lambda2 - c3^2 + 1)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2 - (c3 + c4)*((2*m*(x19*(refyddot - (c11 + c12)*(refydot - x5) + c11*(refydot - x5) - 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9))) + (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 + (m*((c11 + c12)*(refydot - x5) - refyddot - c11*(refydot - x5) + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + x19*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) + (c7 + c8)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2))))/x19^2 - (c3*m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2);
           d2f11dxdp(19,12)=(c3 + c4)*((2*m*(x19*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(c11*(refydot - x5) - refyddot + 4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)) + x19*(4*cos(x9)*sin(x8) - 4*cos(x8)*sin(x7)*sin(x9)) + lambda6*(refy - x2) + (c7 + c8)*(refydot - x5 + lambda6*x18 + c11*(refy - x2))))/x19^2 + (c3*m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2) + (m*(lambda2 - c3^2 + 1)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2;
           d2f11dxdp(19,14)=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,16)=-(m^2*(c3 + c4)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/(2*x19^3*cos(x7)^2*cos(x8)^2);
           d2f11dxdp(19,18)=(c3 + c4)*((2*m*(x19*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)))*(refy - x2 - c11*x18 + x18*(c11 + c12)))/(4*cos(x7)^2*cos(x8)^2)))/x19^3 - (m*(refydot - x5 + (c11 + c12)*(refy - x2) + (c7 + c8)*(refy - x2 - c11*x18 + x18*(c11 + c12)) - c11*(refy - x2)))/x19^2 + (c3*m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2) + (m*(lambda2 - c3^2 + 1)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2;
        end

        function d2f12dxdp = d2f12dxdp(model, ~, ~)
            d2f12dxdp = zeros(model.numStates, model.numParams);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

            d2f12dxdp(9,5)=- 4*c5 - c6;
           d2f12dxdp(9,6)=-c5;
           d2f12dxdp(9,15)=-1;
           d2f12dxdp(12,5)=-1;
           d2f12dxdp(12,6)=-1;
           d2f12dxdp(15,6)=lambda3;
           d2f12dxdp(15,15)=c6;
        end
        
        function d2f13dxdp = d2f13dxdp(model, state, ref)
            d2f13dxdp = zeros(model.numStates, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;

            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;
            m = model.m;

            x1 = state(1);
            x4 = state(4);
            
            x17 = state(17);
            x19 = state(19);

            d2f13dxdp(1,9)=-(c10*m)/x19;
           d2f13dxdp(1,10)=-(c9*m)/x19;
           d2f13dxdp(1,17)=-m/x19;
           d2f13dxdp(4,9)=-m/x19;
           d2f13dxdp(4,10)=-m/x19;
           d2f13dxdp(17,10)=(lambda5*m)/x19;
           d2f13dxdp(17,17)=(c10*m)/x19;
           d2f13dxdp(19,9)=-(m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2;
           d2f13dxdp(19,10)=-(m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2;
           d2f13dxdp(19,17)=-(m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2;

        end

        function d2f14dxdp = d2f14dxdp(model, state, ref)
            d2f14dxdp = zeros(model.numStates, model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;

            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;
            m = model.m;

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);
            
            d2f14dxdp(2,11)=(c12*m)/x19;
           d2f14dxdp(2,12)=(c11*m)/x19;
           d2f14dxdp(2,18)=m/x19;
           d2f14dxdp(5,11)=m/x19;
           d2f14dxdp(5,12)=m/x19;
           d2f14dxdp(18,12)=-(lambda6*m)/x19;
           d2f14dxdp(18,18)=-(c12*m)/x19;
           d2f14dxdp(19,11)=(m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2;
           d2f14dxdp(19,12)=(m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2;
           d2f14dxdp(19,18)=(m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2;
        end
        
        function d2f19dxdp = d2f19dxdp(model, state, ref)
            d2f19dxdp = zeros(model.numStates, model.numParams);
            
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;
            
            m = model.m;
            g = model.g;
            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);
            
            d2f19dxdp(3,7)=-(c8*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(3,8)=-(m*(c7*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) + lambda4*cos(x7)*cos(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(3,16)=-(m*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,7)=-(m*(c8*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,8)=-(m*(c7*cos(x7)*cos(x8) + x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,16)=-m/(4*cos(x7)*cos(x8));
           d2f19dxdp(7,7)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(7,8)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x8)*sin(x7)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(7,16)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x8)*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(8,7)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(8,8)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*sin(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) + (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(8,16)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(10,7)=(m*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(10,8)=(m*sin(x7)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(10,16)=(m*sin(x7)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(11,7)=(m*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(11,8)=(m*sin(x8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(11,16)=(m*sin(x8)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(16,8)=(lambda4*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(16,16)=(c8*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(19,7)=-1;
           d2f19dxdp(19,8)=-1;
        end
        
        function model = nudgeParams(model, delta, index)
            if index == 1
                model.c1 = model.c1 + delta;
            elseif index == 2
                model.c2 = model.c2+ delta;
            elseif index == 3
                model.c3 = model.c3+ delta;
            elseif index == 4
                model.c4 = model.c4+ delta;
            elseif index == 5
                model.c5 = model.c5+ delta;
            elseif index == 6
                model.c6 = model.c6+ delta;
            elseif index == 7
                model.c7 = model.c7+ delta;
            elseif index == 8
                model.c8 = model.c8+ delta;
            elseif index == 9
                model.c9 = model.c9+ delta;
            elseif index == 10
                model.c10 = model.c10+ delta;
            elseif index == 11
                model.c11 = model.c11+ delta;
            elseif index == 12
                model.c12 = model.c12+ delta;
            elseif index == 13
                model.lambda1 = model.lambda1 + delta;
            elseif index == 14
                model.lambda2 = model.lambda2 + delta;
            elseif index == 15
                model.lambda3 = model.lambda3 + delta;
            elseif index == 16
                model.lambda4 = model.lambda4 + delta;
            elseif index == 17
                model.lambda5 = model.lambda5 + delta;
            elseif index == 18
                model.lambda6 = model.lambda6 + delta;
            end
        end
    end
end


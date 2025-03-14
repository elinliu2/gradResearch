classdef derivatives_drone_IB_hybrid
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(19, model.numStates);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = derivatives_drone_IB_hybrid.df4dx(model, state, ref);
            dfdx(5, :) = derivatives_drone_IB_hybrid.df5dx(model, state, ref);
            dfdx(6, :) = derivatives_drone_IB_hybrid.df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = derivatives_drone_IB_hybrid.df10dx(model, state, ref);
            dfdx(11, :) = derivatives_drone_IB_hybrid.df11dx(model, state, ref);
            dfdx(12, :) = derivatives_drone_IB_hybrid.df12dx(model, state, ref);

            dfdx(13, :) = derivatives_drone_IB_hybrid.df13dx(model, state, ref);
            dfdx(14, :) = derivatives_drone_IB_hybrid.df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;

            dfdx(19, :) = derivatives_drone_IB_hybrid.df19dx(model, state, ref);
            
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
            df10dx = zeros(1, model.numStates);

            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;

            m = model.m;
            c1 = model.c1;
            c2 = model.c2;
            c9 = model.c9;
            c10 = model.c10;
            lambda1 = model.lambda1;
            lambda5 = model.lambda5;

            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

            df10dx(1)=- (m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19 - (c1*m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           df10dx(4)=- (m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2)*(c9 + c10))/x19;
           df10dx(7)=c1^2 - lambda1 - c1*(c1 + c2) - 1;
           df10dx(10)=- c1 - c2;
           df10dx(13)=lambda1*(c1 + c2) - c1*lambda1;
           df10dx(17)=- (m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           df10dx(19)=- (m*(lambda1 - c1^2 + 1)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 - (c1*m*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           df10dx(20)=c1 + c2;
           df10dx(22)=1;
        end

        function df11dx = df11dx(model, state, ref)
            df11dx = zeros(1, model.numStates);

            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;

            m = model.m;
            c3 = model.c3;
            c4 = model.c4;
            c11 = model.c11;
            c12 = model.c12;
            lambda2 = model.lambda2;
            lambda6 = model.lambda6;

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

            df11dx(2)=(m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 + (c3*m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           df11dx(5)=(m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4)*(c11 + c12))/x19;
           df11dx(8)=c3^2 - lambda2 - c3*(c3 + c4) - 1;
           df11dx(11)=- c3 - c4;
           df11dx(14)=lambda2*(c3 + c4) - c3*lambda2;
           df11dx(18)=(m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           df11dx(19)=(m*(lambda2 - c3^2 + 1)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 + (c3*m*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           df11dx(21)=c3 + c4;
           df11dx(23)=1;
            
        end
        function df12dx = df12dx(model, ~, ~)
            df12dx = zeros(1, model.numStates);

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

            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
        
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
            df14dx = zeros(1, model.numStates);

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

            df14dx(2)=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           df14dx(5)=(m*(c11 + c12))/x19;
           df14dx(8)=-1;
           df14dx(18)=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           df14dx(19)=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
        end
        function df19dx = df19dx(model, state, ref)
            df19dx = zeros(1, model.numStates);

            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;
            g = model.g;

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

           df19dx(3)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(6)=(m*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(7)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
           df19dx(8)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
           df19dx(10)=-(m*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8));
           df19dx(11)=-(m*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)^2);
           df19dx(16)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(19)=- c7 - c8;
        end

        function dfdp = dfdp(model, state, t)
            
        end
        function d2fdx2 = d2fdx2(model, state)
            
        end
        function d2fdxdp = d2fdxdp(model, ~)
          
        end
        function model = nudgeParams(model, delta, index)
            
        end
        
    end
end


classdef derivativesDroneIBHybrid
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(model.numX, model.numX);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = derivativesDroneIBHybrid.df4dx(model, state, ref);
            dfdx(5, :) = derivativesDroneIBHybrid.df5dx(model, state, ref);
            dfdx(6, :) = derivativesDroneIBHybrid.df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = derivativesDroneIBHybrid.df10dx(model, state, ref);
            dfdx(11, :) = derivativesDroneIBHybrid.df11dx(model, state, ref);
            dfdx(12, :) = derivativesDroneIBHybrid.df12dx(model, state, ref);

            dfdx(13, :) = derivativesDroneIBHybrid.df13dx(model, state, ref);
            dfdx(14, :) = derivativesDroneIBHybrid.df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;

            dfdx(19, :) = derivativesDroneIBHybrid.df19dx(model, state, ref);
        end

        function df4dx = df4dx(model, state, ~)
            df4dx = zeros(1, model.numX);

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
            df5dx = zeros(1, model.numX);

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
            df6dx = zeros(1, model.numX);

            m = model.m;
            x7 = state(7);
            x8 = state(8);
            x19 = state(19);

            df6dx(7)=-(4*x19*cos(x8)*sin(x7))/m;
            df6dx(8)=-(4*x19*cos(x7)*sin(x8))/m;
            df6dx(19)=(4*cos(x7)*cos(x8))/m;
        end

        function df10dx = df10dx(model, state, ref)
            df10dx = zeros(1, model.numX);

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
           
        end

        function df11dx = df11dx(model, state, ref)
            df11dx = zeros(1, model.numX);

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
           
            
        end
        
        function df12dx = df12dx(model, ~, ~)
            df12dx = zeros(1, model.numX);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

            df12dx(9)=c5^2 - lambda3 - c5*(c5 + c6) - 1;
           df12dx(12)=- c5 - c6;
           df12dx(15)=lambda3*(c5 + c6) - c5*lambda3;
        end

        function df13dx = df13dx(model, state, ref)
            df13dx = zeros(1, model.numX);

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
            df14dx = zeros(1, model.numX);

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
            df19dx = zeros(1, model.numX);

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;
            g = model.g;
         
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

           df19dx(3)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(6)=(m*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(7)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
           df19dx(8)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
           df19dx(10)=-(m*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8));
           df19dx(11)=-(m*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)^2);
           df19dx(16)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           df19dx(19)=- c7 - c8;
        end
        
        function dfdy = dfdy(model, ~, ~)
            dfdy = zeros(model.numX, model.numY);
            dfdy(10, 1)=model.c1 + model.c2;
            dfdy(10, 3)=1;

            dfdy(11, 2)=model.c3 + model.c4;
            dfdy(11, 4)=1;
        end
        
        function dfdq = dfdq(model, state, ref)
            dfdq = zeros(model.numX, model.numStates + model.numParams);

            dfdq(10, :) = derivativesDroneIBHybrid.df10dq(model, state, ref);
            dfdq(11, :) = derivativesDroneIBHybrid.df11dq(model, state, ref);
            dfdq(12, :) = derivativesDroneIBHybrid.df12dq(model, state, ref);

            dfdq(13, :) = derivativesDroneIBHybrid.df13dq(model, state, ref);
            dfdq(14, :) = derivativesDroneIBHybrid.df14dq(model, state, ref);

            dfdq(19, :) = derivativesDroneIBHybrid.df19dq(model, state, ref);
        end

        function df10dq = df10dq(model, state, ref)
            df10dq = zeros(1, model.numStates + model.numParams);

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
            x7 = state(7);
            x10 = state(10);
            x13 = state(13);
            x17 = state(17);
            x19 = state(19);
            x20 = state(20);

           df10dq(24)=x20 - x10 - (c1 + c2)*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19) + c1*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19);
           df10dq(25)=x20 - x10 + lambda1*x13 - c1*(x7 - (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19);
           df10dq(32)=(m*(lambda1 - c1^2 + 1)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19 + (c1*m*(c1 + c2)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19;
           df10dq(33)=(m*(lambda1 - c1^2 + 1)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19 + (c1*m*(c1 + c2)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19;
           df10dq(36)=x13*(c1 + c2) - c1*x13 - x7 + (m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19;
           df10dq(40)=(m*(lambda1 - c1^2 + 1)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19 + (c1*m*(c1 + c2)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19;
        end

        function df11dq = df11dq(model, state, ref)
            df11dq = zeros(1, model.numStates + model.numParams);

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
            x8 = state(8);
            x11 = state(11);
            x14 = state(14);
            x18 = state(18);
            x19 = state(19);
            x21 = state(21);

           df11dq(26)=x21 - x11 - (c3 + c4)*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19) + c3*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19);
           df11dq(27)=x21 - x11 + lambda2*x14 - c3*(x8 + (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19);
           df11dq(34)=- (m*(lambda2 - c3^2 + 1)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19 - (c3*m*(c3 + c4)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19;
           df11dq(35)=- (m*(lambda2 - c3^2 + 1)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19 - (c3*m*(c3 + c4)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19;
           df11dq(37)=x14*(c3 + c4) - c3*x14 - x8 - (m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19;
           df11dq(41)=- (m*(lambda2 - c3^2 + 1)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19 - (c3*m*(c3 + c4)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19;
        end

        function df12dq = df12dq(model, state, ref)
            df12dq = zeros(1, model.numStates + model.numParams);

            refpsi = ref.psi;
            refpsidot = ref.psi_dot;

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;
                   
            x9 = state(9);
            x12 = state(12);
            x15 = state(15);

           df12dq(28)=refpsidot - x12 + (c5 + c6)*(refpsi - x9) - c5*(refpsi - x9);
           df12dq(29)=refpsidot - x12 + lambda3*x15 + c5*(refpsi - x9);
           df12dq(38)=refpsi - x9 - c5*x15 + x15*(c5 + c6);
        end

        function df13dq = df13dq(model, state, ref)
            df13dq = zeros(1, model.numStates + model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;

            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
        
            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

           df13dq(32)=(m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19;
           df13dq(33)=(m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19;
           df13dq(40)=(m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19;
        end

        function df14dq = df14dq(model, state, ref)
            df14dq = zeros(1, model.numStates + model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;

            m = model.m;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 
        
            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

           df14dq(34)=-(m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19;
           df14dq(35)=-(m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19;
           df14dq(41)=-(m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19;

        end
        
        function df19dq = df19dq(model, state, ref)
            df19dq = zeros(1, model.numStates + model.numParams);
           
            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;
            g = model.g;

            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

           df19dq(30)=-(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(4*cos(x7)^2*cos(x8)^2);
           df19dq(31)=-(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(4*cos(x7)^2*cos(x8)^2);
           df19dq(39)=-(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2);
        end
        
        function dgdx = dgdx(model, state, ref)
            dgdx = zeros(2, model.numX);

            refx = ref.x;
            refxdot = ref.x_dot;
            refxddot = ref.x_ddot;

            refy = ref.y;
            refydot = ref.y_dot;
            refyddot = ref.y_ddot;
        
            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;
        
            x1 = state(1);
            x4 = state(4);
            x17 = state(17);

            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

           dgdx(1,1)=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           dgdx(1,4)=-(m*(c9 + c10))/x19;
           dgdx(1,17)=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           dgdx(1,19)=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;

           dgdx(2,2)=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           dgdx(2,5)=(m*(c11 + c12))/x19;
           dgdx(2,18)=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           dgdx(2,19)=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;

        end

        function dgdq = dgdq(model, state, ref)
            dgdq = zeros(2, model.numStates + model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
            refy = ref.y;
            refydot = ref.y_dot;

            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 
        
            x1 = state(1);
            x2 = state(2);
            x4 = state(4);
            x5 = state(5);
            x17 = state(17);
            x18 = state(18);
            x19 = state(19);
            
           dgdq(1,32)=(m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19;
           dgdq(1,33)=(m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19;
           dgdq(1,40)=(m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19;
           
           dgdq(2,34)=-(m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19;
           dgdq(2,35)=-(m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19;
           dgdq(2,41)=-(m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19;
        end

        function d2fdx2 = d2fdx2(model, state, ref)
            d2fdx2 = zeros(model.numX, model.numX, model.numX);

            d2fdx2(4, :, :) = derivativesDroneIBHybrid.d2f4dx2(model, state, ref);
            d2fdx2(5, :, :) = derivativesDroneIBHybrid.d2f5dx2(model, state, ref);
            d2fdx2(6, :, :) = derivativesDroneIBHybrid.d2f6dx2(model, state, ref);

            d2fdx2(10, :, :) = derivativesDroneIBHybrid.d2f10dx2(model, state, ref);
            d2fdx2(11, :, :) = derivativesDroneIBHybrid.d2f11dx2(model, state, ref);
            
            d2fdx2(13, :, :) = derivativesDroneIBHybrid.d2f13dx2(model, state, ref);
            d2fdx2(14, :, :) = derivativesDroneIBHybrid.d2f14dx2(model, state, ref);
            d2fdx2(19, :, :) = derivativesDroneIBHybrid.d2f19dx2(model, state, ref);

        end

        function d2f4dx2 = d2f4dx2(model, state, ~)
            d2f4dx2 = zeros(model.numX, model.numX);

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
            d2f5dx2 = zeros(model.numX, model.numX);

            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x19 = state(19);

            m = model.m;

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
            d2f6dx2 = zeros(model.numX, model.numX);

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
            d2f10dx2 = zeros(model.numX, model.numX);

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

           d2f10dx2(1,19)=(m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2 + (c1*m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f10dx2(4,19)=(m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19^2 + (c1*m*(c1 + c2)*(c9 + c10))/x19^2;
           d2f10dx2(17,19)=(m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19^2 + (c1*m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f10dx2(19,1)=(m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2 + (c1*m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f10dx2(19,4)=(m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19^2 + (c1*m*(c1 + c2)*(c9 + c10))/x19^2;
           d2f10dx2(19,17)=(m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19^2 + (c1*m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f10dx2(19,19)=(2*m*(lambda1 - c1^2 + 1)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3 + (2*c1*m*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3;

        end

        function d2f11dx2 = d2f11dx2(model, state, ref)
            d2f11dx2 = zeros(model.numX, model.numX);

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

           d2f11dx2(2,19)=- (m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2 - (c3*m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f11dx2(5,19)=- (m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19^2 - (c3*m*(c3 + c4)*(c11 + c12))/x19^2;
           d2f11dx2(18,19)=- (m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19^2 - (c3*m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f11dx2(19,2)=- (m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2 - (c3*m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f11dx2(19,5)=- (m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19^2 - (c3*m*(c3 + c4)*(c11 + c12))/x19^2;
           d2f11dx2(19,18)=- (m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19^2 - (c3*m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f11dx2(19,19)=- (2*m*(lambda2 - c3^2 + 1)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3 - (2*c3*m*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3;

        end

        function d2f13dx2 = d2f13dx2(model, state, ref)
            d2f13dx2 = zeros(model.numX, model.numX);

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

            d2f13dx2(1,19)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f13dx2(4,19)=(m*(c9 + c10))/x19^2;
           d2f13dx2(17,19)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f13dx2(19,1)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2f13dx2(19,4)=(m*(c9 + c10))/x19^2;
           d2f13dx2(19,17)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2f13dx2(19,19)=(2*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3;

        end

        function d2f14dx2 = d2f14dx2(model, state, ref)
            d2f14dx2 = zeros(model.numX, model.numX);

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

            d2f14dx2(2,19)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f14dx2(5,19)=-(m*(c11 + c12))/x19^2;
           d2f14dx2(18,19)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f14dx2(19,2)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2f14dx2(19,5)=-(m*(c11 + c12))/x19^2;
           d2f14dx2(19,18)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2f14dx2(19,19)=-(2*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3;

        end

        function d2f19dx2 = d2f19dx2(model, state, ref)
            d2f19dx2 = zeros(model.numX, model.numX);

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;
            g = model.g;
         
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

            d2f19dx2(3,7)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(3,8)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(3,10)=(m*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(3,11)=(m*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(6,7)=(m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) + cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(6,8)=(m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) + cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(6,10)=(m*sin(x7)*(c7 + c8))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(6,11)=(m*sin(x8)*(c7 + c8))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(7,3)=(m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x8)*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,6)=(m*((c7 + c8)*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)) + cos(x8)*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,7)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m - (8*x19*cos(x8)^2*sin(x7)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^2) - (3*m*sin(x7)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^4*cos(x8)^2) - (m*sin(x7)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(cos(x7)^3*cos(x8)^2);
           d2f19dx2(7,8)=(m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(2*cos(x7)^2*cos(x8)^3) - (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(2*cos(x7)^3*cos(x8)^2) - (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(cos(x7)^3*cos(x8)^3);
           d2f19dx2(7,10)=- (m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) - (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)^3*cos(x8));
           d2f19dx2(7,11)=-(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(7,16)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(8,3)=(m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) - cos(x7)*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,6)=(m*((c7 + c8)*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)) + cos(x7)*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,7)=(m*((x10*cos(x7)*sin(x8) + x11*cos(x8)*sin(x7))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + sin(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (12*x19*cos(x7)*cos(x8)*sin(x7)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(2*cos(x7)^2*cos(x8)^3) - (m*sin(x7)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(2*cos(x7)^3*cos(x8)^2) - (m*sin(x7)*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(cos(x7)^3*cos(x8)^3);
           d2f19dx2(8,8)=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) + (4*x19*cos(x7)^2*cos(x8)^2*(c7 + c8))/m - (8*x19*cos(x7)^2*sin(x8)^2*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^2) - (3*m*sin(x8)^2*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^4) - (m*sin(x8)*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(cos(x7)^2*cos(x8)^3);
           d2f19dx2(8,10)=-(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(8,11)=- (m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) - (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)*cos(x8)^3);
           d2f19dx2(8,16)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(10,3)=(m*sin(x7)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(10,6)=(m*sin(x7)*(c7 + c8))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(10,7)=- (m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) - (m*sin(x7)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)^3*cos(x8));
           d2f19dx2(10,8)=-(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(10,16)=(m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(11,3)=(m*sin(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(11,6)=(m*sin(x8)*(c7 + c8))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(11,7)=-(m*sin(x7)*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dx2(11,8)=- (m*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)) - (m*sin(x8)^2*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(2*cos(x7)*cos(x8)^3);
           d2f19dx2(11,16)=(m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);
           d2f19dx2(16,7)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dx2(16,8)=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2) + (m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dx2(16,10)=(m*sin(x7)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dx2(16,11)=(m*sin(x8)*(c7*lambda4 - lambda4*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);

        end

        function d2thetarefdx2 = d2thetarefdx2(model, state, ref)
            d2thetarefdx2 = zeros(model.numX, model.numX);

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

           d2thetarefdx2(1,19)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2thetarefdx2(4,19)=(m*(c9 + c10))/x19^2;
           d2thetarefdx2(17,19)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2thetarefdx2(19,1)=(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19^2;
           d2thetarefdx2(19,4)=(m*(c9 + c10))/x19^2;
           d2thetarefdx2(19,17)=(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19^2;
           d2thetarefdx2(19,19)=(2*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^3;
        end

        function d2phirefdx2 = d2phirefdx2(model, state, ref)
            d2phirefdx2 = zeros(model.numX, model.numX);

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

            d2phirefdx2(2,19)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2phirefdx2(5,19)=-(m*(c11 + c12))/x19^2;
           d2phirefdx2(18,19)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2phirefdx2(19,2)=-(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19^2;
           d2phirefdx2(19,5)=-(m*(c11 + c12))/x19^2;
           d2phirefdx2(19,18)=-(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19^2;
           d2phirefdx2(19,19)=-(2*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^3;
        end
        
        function d2fdxdp = d2fdxdp(model, state, ref)
            d2fdxdp = zeros(model.numX, model.numX, model.numParams);

            d2fdxdp(10, :, :) = derivativesDroneIBHybrid.d2f10dxdp(model, state, ref);
            d2fdxdp(11, :, :) = derivativesDroneIBHybrid.d2f11dxdp(model, state, ref);
            d2fdxdp(12, :, :) = derivativesDroneIBHybrid.d2f12dxdp(model, state, ref);

            d2fdxdp(13, :, :) = derivativesDroneIBHybrid.d2f13dxdp(model, state, ref);
            d2fdxdp(14, :, :) = derivativesDroneIBHybrid.d2f14dxdp(model, state, ref);
            d2fdxdp(19, :, :) = derivativesDroneIBHybrid.d2f19dxdp(model, state, ref);

        end
        
        function d2f10dxdp = d2f10dxdp(model, state, ref)
            d2f10dxdp = zeros(model.numX, model.numParams);

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

           d2f10dxdp(1,1)=(c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19 - (m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,2)=-(c1*m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,9)=- (c10*m*(lambda1 - c1^2 + 1))/x19 - (c1*c10*m*(c1 + c2))/x19;
           d2f10dxdp(1,10)=- (c9*m*(lambda1 - c1^2 + 1))/x19 - (c1*c9*m*(c1 + c2))/x19;
           d2f10dxdp(1,13)=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
           d2f10dxdp(1,17)=- (m*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2))/x19;
           d2f10dxdp(4,1)=(c1*m*(c9 + c10))/x19 - (m*(c1 + c2)*(c9 + c10))/x19;
           d2f10dxdp(4,2)=-(c1*m*(c9 + c10))/x19;
           d2f10dxdp(4,9)=- (m*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2))/x19;
           d2f10dxdp(4,10)=- (m*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2))/x19;
           d2f10dxdp(4,13)=-(m*(c9 + c10))/x19;
           d2f10dxdp(7,1)=-c2;
           d2f10dxdp(7,2)=-c1;
           d2f10dxdp(7,13)=-1;
           d2f10dxdp(10,1)=-1;
           d2f10dxdp(10,2)=-1;
           d2f10dxdp(13,2)=lambda1;
           d2f10dxdp(13,13)=c2;
           d2f10dxdp(17,1)=(c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19 - (m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           d2f10dxdp(17,2)=-(c1*m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           d2f10dxdp(17,10)=(lambda5*m*(lambda1 - c1^2 + 1))/x19 + (c1*lambda5*m*(c1 + c2))/x19;
           d2f10dxdp(17,13)=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
           d2f10dxdp(17,17)=(c10*m*(lambda1 - c1^2 + 1))/x19 + (c1*c10*m*(c1 + c2))/x19;
           d2f10dxdp(19,1)=(c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 - (m*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,2)=-(c1*m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,9)=- (m*(lambda1 - c1^2 + 1)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2 - (c1*m*(c1 + c2)*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2;
           d2f10dxdp(19,10)=- (m*(lambda1 - c1^2 + 1)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2 - (c1*m*(c1 + c2)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2;
           d2f10dxdp(19,13)=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
           d2f10dxdp(19,17)=- (m*(lambda1 - c1^2 + 1)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2 - (c1*m*(c1 + c2)*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2;
        end

        function d2f11dxdp = d2f11dxdp(model, state, ref)
            d2f11dxdp = zeros(model.numX, model.numParams);

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

           d2f11dxdp(2,3)=(m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 - (c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           d2f11dxdp(2,4)=(c3*m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           d2f11dxdp(2,11)=(c12*m*(lambda2 - c3^2 + 1))/x19 + (c3*c12*m*(c3 + c4))/x19;
           d2f11dxdp(2,12)=(c11*m*(lambda2 - c3^2 + 1))/x19 + (c3*c11*m*(c3 + c4))/x19;
           d2f11dxdp(2,14)=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
           d2f11dxdp(2,18)=(m*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4))/x19;
           d2f11dxdp(5,3)=(m*(c3 + c4)*(c11 + c12))/x19 - (c3*m*(c11 + c12))/x19;
           d2f11dxdp(5,4)=(c3*m*(c11 + c12))/x19;
           d2f11dxdp(5,11)=(m*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4))/x19;
           d2f11dxdp(5,12)=(m*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4))/x19;
           d2f11dxdp(5,14)=(m*(c11 + c12))/x19;
           d2f11dxdp(8,3)=-c4;
           d2f11dxdp(8,4)=-c3;
           d2f11dxdp(8,14)=-1;
           d2f11dxdp(11,3)=-1;
           d2f11dxdp(11,4)=-1;
           d2f11dxdp(14,4)=lambda2;
           d2f11dxdp(14,14)=c4;
           d2f11dxdp(18,3)=(m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19 - (c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           d2f11dxdp(18,4)=(c3*m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           d2f11dxdp(18,12)=- (lambda6*m*(lambda2 - c3^2 + 1))/x19 - (c3*lambda6*m*(c3 + c4))/x19;
           d2f11dxdp(18,14)=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
           d2f11dxdp(18,18)=- (c12*m*(lambda2 - c3^2 + 1))/x19 - (c3*c12*m*(c3 + c4))/x19;
           d2f11dxdp(19,3)=(m*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 - (c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,4)=(c3*m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,11)=(m*(lambda2 - c3^2 + 1)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2 + (c3*m*(c3 + c4)*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2;
           d2f11dxdp(19,12)=(m*(lambda2 - c3^2 + 1)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2 + (c3*m*(c3 + c4)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2;
           d2f11dxdp(19,14)=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
           d2f11dxdp(19,18)=(m*(lambda2 - c3^2 + 1)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2 + (c3*m*(c3 + c4)*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2;
        end

        function d2f12dxdp = d2f12dxdp(model, ~, ~)
            d2f12dxdp = zeros(model.numX, model.numParams);

             c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

           d2f12dxdp(9,5)=-c6;
           d2f12dxdp(9,6)=-c5;
           d2f12dxdp(9,15)=-1;
           d2f12dxdp(12,5)=-1;
           d2f12dxdp(12,6)=-1;
           d2f12dxdp(15,6)=lambda3;
           d2f12dxdp(15,15)=c6;

        end

        function d2f13dxdp = d2f13dxdp(model, state, ref)
            d2f13dxdp = zeros(model.numX, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
           
            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
        
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
            d2f14dxdp = zeros(model.numX, model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;

            m = model.m;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 
        
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
            d2f19dxdp = zeros(model.numX, model.numParams);

            refz = ref.z;
            refzdot = ref.z_dot;
            refzddot = ref.z_ddot;

            c7 = model.c7; 
            c8 = model.c8;
            lambda4 = model.lambda4;
            m = model.m;
            g = model.g;
         
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);
            x16 = state(16);
            x19 = state(19);

           d2f19dxdp(3,7)=(c8*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(3,8)=(m*(c7*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - lambda4*cos(x7)*cos(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(3,16)=(m*(x10*cos(x8)*sin(x7) - c8*cos(x7)*cos(x8) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,7)=(m*(x10*cos(x8)*sin(x7) - c8*cos(x7)*cos(x8) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,8)=(m*(x10*cos(x8)*sin(x7) - c7*cos(x7)*cos(x8) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(6,16)=-m/(4*cos(x7)*cos(x8));
           d2f19dxdp(7,7)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x8)*sin(x7)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) - (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(7,8)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x8)*sin(x7)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - (4*x19*cos(x7)*cos(x8)^2*sin(x7))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(7,16)=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x8)*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
           d2f19dxdp(8,7)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) + cos(x7)*sin(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m) - (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + (c7 + c8)*(refzdot - x6) - c7*(refzdot - x6) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(8,8)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + cos(x7)*sin(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - (4*x19*cos(x7)^2*cos(x8)*sin(x8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) - cos(x7)*cos(x8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m)))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(8,16)=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) + cos(x7)*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(refz - x3 - c7*x16 + x16*(c7 + c8)) - cos(x7)*cos(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
           d2f19dxdp(10,7)=-(m*sin(x7)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(10,8)=-(m*sin(x7)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(10,16)=-(m*sin(x7)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*cos(x7)^2*cos(x8));
           d2f19dxdp(11,7)=-(m*sin(x8)*(refzdot - x6 + (c7 + c8)*(refz - x3) - c7*(refz - x3)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(11,8)=-(m*sin(x8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(11,16)=-(m*sin(x8)*(refz - x3 - c7*x16 + x16*(c7 + c8)))/(4*cos(x7)*cos(x8)^2);
           d2f19dxdp(16,8)=-(lambda4*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(16,16)=-(c8*m*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
           d2f19dxdp(19,7)=-1;
           d2f19dxdp(19,8)=-1;

        end

        function d2fdydp = d2fdydp(model, ~, ~)
            d2fdydp = zeros(model.numX, model.numY, model.numParams);
            d2fdydp(10, 1, 1)=1;
            d2fdydp(10, 1, 2)=1;
       
            d2fdydp(11, 2, 3)=1;
            d2fdydp(11, 2, 4)=1;
        end
        
        function d2thetarefdxdp = d2thetarefdxdp(model, state, ref)
            d2thetarefdxdp = zeros(model.numX, model.numParams);

            refx = ref.x;
            refxdot = ref.x_dot;
        
            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
        
            x1 = state(1);
            x4 = state(4);
            x17 = state(17);
            x19 = state(19);

           d2thetarefdxdp(1,9)=-(c10*m)/x19;
           d2thetarefdxdp(1,10)=-(c9*m)/x19;
           d2thetarefdxdp(1,17)=-m/x19;
           d2thetarefdxdp(4,9)=-m/x19;
           d2thetarefdxdp(4,10)=-m/x19;
           d2thetarefdxdp(17,10)=(lambda5*m)/x19;
           d2thetarefdxdp(17,17)=(c10*m)/x19;
           d2thetarefdxdp(19,9)=-(m*(refxdot - x4 + (c9 + c10)*(refx - x1) - c9*(refx - x1)))/x19^2;
           d2thetarefdxdp(19,10)=-(m*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)))/x19^2;
           d2thetarefdxdp(19,17)=-(m*(refx - x1 - c9*x17 + x17*(c9 + c10)))/x19^2;


        end

        function d2phirefdxdp = d2phirefdxdp(model, state, ref)
            d2phirefdxdp = zeros(model.numX, model.numParams);

            refy = ref.y;
            refydot = ref.y_dot;
        
            m = model.m;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 
        
            x2 = state(2);
            x5 = state(5);
            x18 = state(18);
            x19 = state(19);

           d2phirefdxdp(2,11)=(c12*m)/x19;
           d2phirefdxdp(2,12)=(c11*m)/x19;
           d2phirefdxdp(2,18)=m/x19;
           d2phirefdxdp(5,11)=m/x19;
           d2phirefdxdp(5,12)=m/x19;
           d2phirefdxdp(18,12)=-(lambda6*m)/x19;
           d2phirefdxdp(18,18)=-(c12*m)/x19;
           d2phirefdxdp(19,11)=(m*(refydot - x5 + (c11 + c12)*(refy - x2) - c11*(refy - x2)))/x19^2;
           d2phirefdxdp(19,12)=(m*(refydot - x5 + lambda6*x18 + c11*(refy - x2)))/x19^2;
           d2phirefdxdp(19,18)=(m*(refy - x2 - c11*x18 + x18*(c11 + c12)))/x19^2;


        end
        
        function [model, state] = nudgeParams(model, state, delta, index)
            if index <= 23
                state(index) = state(index) + delta;
            elseif index == 24 
                model.c1 = model.c1 + delta;
            elseif index == 25 
                model.c2 = model.c2 + delta;
            elseif index == 26
                model.c3 = model.c3 + delta;
            elseif index == 27
                model.c4 = model.c4 + delta;
            elseif index == 28
                model.c5 = model.c5 + delta;
            elseif index == 29
                model.c6 = model.c6 + delta;
            elseif index == 30
                model.c7 = model.c7 + delta;
            elseif index == 31
                model.c8 = model.c8 + delta;
            elseif index == 32
                model.c9 = model.c9 + delta;
            elseif index == 33
                model.c10 = model.c10 + delta;
            elseif index == 34
                model.c11 = model.c11 + delta;
            elseif index == 35
                model.c12 = model.c12 + delta;
            elseif index == 36
                model.lambda1 = model.lambda1 + delta;
            elseif index == 37
                model.lambda2 = model.lambda2 + delta;
            elseif index == 38
                model.lambda3 = model.lambda3 + delta;
            elseif index == 39
                model.lambda4 = model.lambda4 + delta;
            elseif index == 40
                model.lambda5 = model.lambda5 + delta;
            elseif index == 41
                model.lambda6 = model.lambda6 + delta;
            end
        end
    end
end


classdef derivatives_drone_IB
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(model.numStates, model.numStates);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = derivatives_drone_IB.df4dx(model, state, ref);
            dfdx(5, :) = derivatives_drone_IB.df5dx(model, state, ref);
            dfdx(6, :) = derivatives_drone_IB.df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = derivatives_drone_IB.df10dx(model, state, ref);
            dfdx(11, :) = derivatives_drone_IB.df11dx(model, state, ref);
            dfdx(12, :) = derivatives_drone_IB.df12dx(model, state, ref);

            dfdx(13, :) = derivatives_drone_IB.df13dx(model, state, ref);
            dfdx(14, :) = derivatives_drone_IB.df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;
        end

        function dUzdx = dUzdx(model, state, ref)
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);

            g = model.g;
            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;

            dUzdx(3) = m/(4*cos(x7)*cos(x8))*((1-c7^2+lambda4)*(-1) + (c7+c8)*(-c7));
            dUzdx(6) = m/(4*cos(x7)*cos(x8))*((c7+c8)*(-1));
            dUzdx(7) = m/(4*cos(x8))*     (g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16) * tan(x7)*sec(x7);
            dUzdx(8) = m/(4*cos(x7))*     (g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16) * tan(x8)*sec(x8);
            dUzdx(16) = m/(4*cos(x7)*cos(x8))*((c7+c8)*lambda4 - c7*lambda4);
        end

        function df4dx = df4dx(model, state, ref)
            df4dx = zeros(1, model.numStates);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);

            m = model.m;

            f4_trig = sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9);
            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);

            df4dx(3) = dU_zdx(3)*f4_trig;
            df4dx(6) = dU_zdx(6)*f4_trig;
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);

            df4dx(7) = 1/m*(dU_zdx(7)*f4_trig + U_z*(cos(x7)*cos(x8)*cos(x9)));
            df4dx(8) = 1/m*(dU_zdx(8)*f4_trig + U_z*(cos(x8)*sin(x9) - sin(x7)*sin(x8)*cos(x9)));

            df4dx(9) = 1/m*(U_z*(sin(x8)*cos(x9)-sin(x7)*cos(x8)*sin(x9)));

            df4dx(16) = dU_zdx(16)*f4_trig;
        end

        function df5dx = df5dx(model, state, ref)
            df5dx = zeros(1, model.numStates);
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x16 = state(16);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            f5_trig = -sin(x8)*cos(x9) + sin(x7)*cos(x8)*sin(x9);
            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);

            df5dx(3) = dU_zdx(3)*f5_trig;
            df5dx(6) = dU_zdx(6)*f5_trig;
            
            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;
            U_z = m/(4*cos(x7)*cos(x8))*(g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            df5dx(7) = 1/m*(dU_zdx(7)*f5_trig + U_z*(cos(x7)*cos(x8)*sin(x9)));
            df5dx(8) = 1/m*(dU_zdx(8)*f5_trig + U_z*(-cos(x8)*cos(x9) - sin(x7)*sin(x8)*sin(x9)));

            df5dx(9) = 1/m*(U_z*(sin(x8)*sin(x9)+sin(x7)*cos(x8)*cos(x9)));

            df5dx(16) = dU_zdx(16)*f5_trig;
        end
        function df6dx = df6dx(model, state, ref)
            df6dx = zeros(1, model.numStates);
            x7 = state(7);
            x8 = state(8);

            m = model.m;

            f6_trig = cos(x7)*cos(x8);
            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);

            df6dx(3) = dU_zdx(3)*f6_trig;
            df6dx(6) = dU_zdx(6)*f6_trig;
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            df6dx(7) = 1/m*(dU_zdx(7)*f6_trig + U_z*(-sin(x7)*cos(x8)));
            df6dx(8) = 1/m*(dU_zdx(8)*f6_trig + U_z*(-cos(x7)*sin(x8)));

            df6dx(16) = dU_zdx(16)*f6_trig;
        end

        function de8dx = de8dx(model)
            de8dx = zeros(1, model.numStates);
            de8dx(4) = -model.c7;
            de8dx(6) = -1;
            de8dx(16) = model.lambda4;
        end

        function de8dotdx = de8dotdx(model, state, ref)
            de8dotdx = derivatives_drone_IB.df6dx(model, state, ref);
            de8dotdx(3) = -model.lambda4;
            de8dotdx(6) = -model.c7;
        end

        function de8ddotdx = de8ddotdx(model, state, ref)
            de8ddotdx = zeros(1, model.numStates);
            c7 = model.c7;
            lambda4 = model.lambda4;

            dz_ddotdx = derivatives_drone_IB.df6dx(model, state, ref);
            dz_dddotdx = derivatives_drone_IB.dzdddotdx(model, state, ref);
            de8ddotdx(3) = c7*(-dz_ddotdx(3))-dz_dddotdx(3);
            de8ddotdx(6) = c7*(-dz_ddotdx(6))- lambda4 - dz_dddotdx(6);
            de8ddotdx(7) = c7*(-dz_ddotdx(7))-dz_dddotdx(7);
            de8ddotdx(8) = c7*(-dz_ddotdx(8))-dz_dddotdx(8);
            de8ddotdx(10) = -dz_dddotdx(10);
            de8ddotdx(11) = -dz_dddotdx(11);
            de8ddotdx(16) = c7*(-dz_ddotdx(16))-dz_dddotdx(16);

        end

        function de10dx = de10dx(model)
            de10dx = zeros(1, model.numStates);
            de10dx(1) = -model.c9;
            de10dx(4) = -1;
            de10dx(17) = model.lambda5;
        end

        function de10dotdx = de10dotdx(model, state, ref)
            de10dotdx = derivatives_drone_IB.df4dx(model, state, ref);
            de10dotdx(1) = -model.lambda5;
            de10dotdx(4) = -model.c9;
        end

        function de10ddotdx = de10ddotdx(model, state, ref)
            de10ddotdx = zeros(1, model.numStates);
            c9 = model.c9;
            lambda5 = model.lambda5;
            dxddot_dx = derivatives_drone_IB.df4dx(model, state, ref);
            dxdddot_dx = derivatives_drone_IB.dxdddotdx(model, state, ref);
            de10ddotdx(3) = -c9*dxddot_dx(3) - dxdddot_dx(3);
            de10ddotdx(4) = -c9*dxddot_dx(4) - lambda5 - dxdddot_dx(4);
            de10ddotdx(6) = -c9*dxddot_dx(6) - dxdddot_dx(6);
            de10ddotdx(7) = -c9*dxddot_dx(7) - dxdddot_dx(7);
            de10ddotdx(8) = -c9*dxddot_dx(8) - dxdddot_dx(8);
            de10ddotdx(9) = -c9*dxddot_dx(9) - dxdddot_dx(9);
            de10ddotdx(16) = -c9*dxddot_dx(16) - dxdddot_dx(16);

            de10ddotdx(10) = - dxdddot_dx(10);
            de10ddotdx(11) = - dxdddot_dx(11);
            de10ddotdx(12) = - dxdddot_dx(12);
        end

        function de12dx = de12dx(model)
            de12dx = zeros(1, model.numStates);
            de12dx(2) = -model.c11;
            de12dx(5) = -1;
            de12dx(18) = model.lambda6;
        end

        function de12dotdx = de12dotdx(model, state, ref)
            de12dotdx =  derivatives_drone_IB.df4dx(model, state, ref);
            de12dotdx(2) = -model.lambda6;
            de12dotdx(5) = -model.c11;
        end

        function de12ddotdx = de12ddotdx(model, state, ref)
            de12ddotdx = zeros(1, model.numStates);

            c11 = model.c11;
            lambda6 = model.lambda5;
            dyddot_dx = derivatives_drone_IB.df5dx(model, state, ref);
            dydddot_dx = derivatives_drone_IB.dydddotdx(model, state, ref);

            de12ddotdx(3) = -c11*dyddot_dx(3) - dydddot_dx(3);
            de12ddotdx(5) = -c11*dyddot_dx(5) - dydddot_dx(5) - lambda6;
            de12ddotdx(6) = -c11*dyddot_dx(6) - dydddot_dx(6);
            de12ddotdx(7) = -c11*dyddot_dx(7) - dydddot_dx(7);
            de12ddotdx(8) = -c11*dyddot_dx(8) - dydddot_dx(8);
            de12ddotdx(9) = -c11*dyddot_dx(9) - dydddot_dx(9);
            de12ddotdx(16) = -c11*dyddot_dx(16) - dydddot_dx(16);

            de12ddotdx(10) = - dydddot_dx(10);
            de12ddotdx(11) = - dydddot_dx(11);
            de12ddotdx(12) = - dydddot_dx(12);
        end

        function dexdx = dexdx(model, ~, ~)
            dexdx = zeros(1, model.numStates);
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 
            de10_dx = derivatives_drone_IB.de10dx(model);
            
            dexdx(1) = (1-c9^2+lambda5)*(-1) + (c9 + c10)*de10_dx(1);
            dexdx(4) = (c9+c10)*de10_dx(4);
            dexdx(17) = (c9+c10)*de10_dx(17) - c9*lambda5;
        end

        function dexdotdx = dexdotdx(model, state, ref)
            dexdotdx = zeros(1, model.numStates);

            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5; 

            dxddot_dx = derivatives_drone_IB.df4dx(model, state, ref);
            de10dot_dx = derivatives_drone_IB.de10dotdx(model, state, ref);
            
            dexdotdx(1) = (c9+c10)*de10dot_dx(1) - c9*lambda5;
            dexdotdx(4) = (1-c9^2+lambda5)*-dxddot_dx(4) + (c9 + c10)*de10dot_dx(4);
        end

        function dexddotdx = dexddotdx(model, state, ref)
            dexddotdx = zeros(1, model.numStates);
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;
            de10ddot_dx = derivatives_drone_IB.de10ddotdx(model, state, ref);
            dxddot_dx = derivatives_drone_IB.df4dx(model, state, ref);

            dexddotdx(3) = (1-c9^2+lambda5)*-dxddot_dx(3) + (c9+c10)*de10ddot_dx(3);
            dexddotdx(4) = -c9*lambda5;
            dexddotdx(6) = (1-c9^2+lambda5)*-dxddot_dx(6) + (c9+c10)*de10ddot_dx(6);
            dexddotdx(7) = (1-c9^2+lambda5)*-dxddot_dx(7) + (c9+c10)*de10ddot_dx(7);
            dexddotdx(8) = (1-c9^2+lambda5)*-dxddot_dx(8) + (c9+c10)*de10ddot_dx(8);
            dexddotdx(9) = (1-c9^2+lambda5)*-dxddot_dx(9) + (c9+c10)*de10ddot_dx(9);
            dexddotdx(10) = (c9+c10)*de10ddot_dx(10);
            dexddotdx(11) = (c9+c10)*de10ddot_dx(11);
            dexddotdx(12) = (c9+c10)*de10ddot_dx(12);
            dexddotdx(16) = (1-c9^2+lambda5)*-dxddot_dx(16) + (c9+c10)*de10ddot_dx(16);
        end

        function deydx = deydx(model, state, ref)
            deydx = zeros(1, model.numStates);
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 
            de12_dx = derivatives_drone_IB.de12dx(model);
            
            deydx(2) = -(1-c11^2+lambda6)*(-1) - (c11 + c12)*de12_dx(2);
            deydx(5) = -(c11+c12)*de12_dx(5);
            deydx(18) = -(c11+c12)*de12_dx(18) + c11*lambda6;
        end

        function deydotdx = deydotdx(model, state, ref)
            deydotdx = zeros(1, model.numStates);

            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6; 

            de12dot_dx = derivatives_drone_IB.de12dotdx(model, state, ref);
            dyddot_dx = derivatives_drone_IB.df5dx(model, state, ref);
            
            deydotdx(2) = -(c11+c12)*de12dot_dx(2) + c11*lambda6;
            deydotdx(5) = -(1-c11^2+lambda6)*-dyddot_dx(5) + (c11 + c12)*de12dot_dx(5);
        end

        function deyddotdx = deyddotdx(model, state, ref)
            deyddotdx = zeros(1, model.numStates);

            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;

            de12ddot_dx = derivatives_drone_IB.de12ddotdx(model, state, ref);
            dyddot_dx = derivatives_drone_IB.df5dx(model, state, ref);

            deyddotdx(3) = -(1-c11^2+lambda6)*-dyddot_dx(3) - (c11+c12)*de12ddot_dx(3);
            deyddotdx(5) = c11*lambda6;
            deyddotdx(6) = -(1-c11^2+lambda6)*-dyddot_dx(6) - (c11+c12)*de12ddot_dx(6);
            deyddotdx(7) = -(1-c11^2+lambda6)*-dyddot_dx(7) - (c11+c12)*de12ddot_dx(7);
            deyddotdx(8) = -(1-c11^2+lambda6)*-dyddot_dx(8) - (c11+c12)*de12ddot_dx(8);
            deyddotdx(9) = -(1-c11^2+lambda6)*-dyddot_dx(9) - (c11+c12)*de12ddot_dx(9);
            deyddotdx(10) = -(c11+c12)*de12ddot_dx(10);
            deyddotdx(11) = -(c11+c12)*de12ddot_dx(11);
            deyddotdx(12) = -(c11+c12)*de12ddot_dx(12);
            deyddotdx(16) = -(1-c11^2+lambda6)*-dyddot_dx(16) - (c11+c12)*de12ddot_dx(16);
        end

        function dezdx = dezdx(model)
            dezdx = (model.c7 + model.c8) * derivatives_drone_IB.de8dx(model);
            dezdx(3) = dezdx(3)-(1-model.c7^2+model.lambda4);
            dezdx(16) = dezdx(16)-model.c7*model.lambda4;
        end

        function dezdotdx = dezdotdx(model, state, ref)
            dezdotdx = (model.c7 + model.c8) * derivatives_drone_IB.de8dotdx(model, state, ref);
            dezdotdx(6) = dezdotdx(6)-(1-model.c7^2+model.lambda4);
            dezdotdx(3) = dezdotdx(3)+model.c7*model.lambda4;
        end

        function dezddotdx = dezddotdx(model, state, ref)
            dezddotdx = zeros(1, model.numStates);

            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            dz_ddotdx = derivatives_drone_IB.df6dx(model, state, ref);
            de8_ddotdx = derivatives_drone_IB.de8ddotdx(model, state, ref);

            dezddotdx(3) = (1-c7^2+lambda4)*(-dz_ddotdx(3)) + (c7+c8)*de8_ddotdx(3);
            dezddotdx(6) = (1-c7^2+lambda4)*(-dz_ddotdx(6)) + (c7+c8)*de8_ddotdx(6) - lambda4;
            dezddotdx(7) = (1-c7^2+lambda4)*(-dz_ddotdx(7)) + (c7+c8)*de8_ddotdx(7);
            dezddotdx(8) = (1-c7^2+lambda4)*(-dz_ddotdx(8)) + (c7+c8)*de8_ddotdx(8);
            dezddotdx(10) = (c7+c8)*de8_ddotdx(10);
            dezddotdx(11) = (c7+c8)*de8_ddotdx(11);
            dezddotdx(16) = (1-c7^2+lambda4)*(-dz_ddotdx(16)) + (c7+c8)*de8_ddotdx(16);
        end

        function dUzdotdx = dUzdotdx(model, state, ref)
            dUzdotdx = zeros(1, model.numStates);
            dez_dotdx = derivatives_drone_IB.dezdotdx(model, state, ref);
            dez_dx = derivatives_drone_IB.dezdx(model);
            
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);
            m = model.m;
            g = model.g;

            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            denom = (cos(x7)*cos(x8))^2;
            dUzdotdx(3) = m/(4*denom) * (dez_dotdx(3)*(cos(x7)*cos(x8)) - dez_dx(3)*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)));
            dUzdotdx(6) = m/(4*denom) * (dez_dotdx(6)*(cos(x7)*cos(x8)));
            dUzdotdx(16) = m/(4*denom) * (-dez_dx(16)*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)));
            
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez_dot = (1-c7^+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            
            Uzdotnum = m/4*(cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)));
            
            dUzdotnumdx7 = m/4*(-sin(x7)*cos(x8)*ez_dot + dez_dotdx(7)*cos(x7)*cos(x8) - ez*(-cos(x7)*cos(x8)+sin(x7)*sin(x8)));
            dUzdotdx(7) = (dUzdotnumdx7*denom-Uzdotnum*(2*cos(x7)*cos(x8)*-sin(x7)))/denom^2;

            dUzdotnumdx8 = m/4*(cos(x7)*-sin(x8)*ez_dot + dez_dotdx(8)*cos(x7)*cos(x8) - ez*(sin(x7)*sin(x8)-cos(x7)*cos(x8)));
            dUzdotdx(8) = (dUzdotnumdx8*denom-Uzdotnum*(2*cos(x7)*cos(x8)*-sin(x8)))/denom^2;
        end

        function dUzdottopdx = dUzdottopdx(model, state, ref)
            dUzdottopdx = zeros(1, model.numStates);
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;
            g = model.g;
            m = model.m;

            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);

            dez_dx = derivatives_drone_IB.dezdx(model);
            dezdot_dx = derivatives_drone_IB.dezdotdx(model, state, ref);

            dUzdottopdx(3) = cos(x7)*cos(x8)*dezdot_dx(3) + dez_dx(3)*(sin(x7)*cos(x8) + cos(x7)*sin(x8));
            dUzdottopdx(6) = cos(x7)*cos(x8)*dezdot_dx(6) ;
            dUzdottopdx(16) = dez_dx(16)*(sin(x7)*cos(x8) + cos(x7)*sin(x8));

            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            
            dUzdottopdx(7) = -sin(x7)*cos(x8)*ez_dot+ez*(cos(x7)*cos(x8) - sin(x7)*sin(x8));
            dUzdottopdx(8) = cos(x7)*-sin(x8)*ez_dot+ez*(cos(x7)*cos(x8) - sin(x7)*sin(x8));
        end

        function dUzdottopdotdx = dUzdottopdotdx(model, state, ref)
            dUzdottopdotdx = zeros(1, model.numStates);

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

            dez_dx = derivatives_drone_IB.dezdx(model);
            dezddot_dx = derivatives_drone_IB.dezddotdx(model, state, ref);

            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            z_dddot = 4/m*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            e8_ddot = c7*(ref.z_ddot-x6)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);

            dUzdottopdotdx(3) = (2*cos(x7)*cos(x8)-2*sin(x7)*sin(x8))*dez_dx(3)+cos(x7)*cos(x8)*dezddot_dx(3);
            dUzdottopdotdx(6) = cos(x7)*cos(x8)*dezddot_dx(6);
            dUzdottopdotdx(10) = cos(x7)*cos(x8)*dezddot_dx(10);
            dUzdottopdotdx(11) = cos(x7)*cos(x8)*dezddot_dx(11);
            dUzdottopdotdx(16) = (2*cos(x7)*cos(x8)-2*sin(x7)*sin(x8))*dez_dx(3)+cos(x7)*cos(x8)*dezddot_dx(16);

            dUzdottopdotdx(7) = dezddot_dx(7)*cos(x7)*cos(x8) - ez_ddot*(-sin(x7)*cos(x8)) ...
                              + dez_dx(7) * 2 * (cos(x7)*cos(x8)-sin(x7)*sin(x8)) ...
                              + ez * 2 * (-sin(x7)*cos(x8) - cos(x7)*sin(x8));
            dUzdottopdotdx(8) = dezddot_dx(8)*cos(x7)*cos(x8) - ez_ddot*(cos(x7)*-sin(x8)) ...
                              + dez_dx(8) * 2 * (cos(x7)*cos(x8)-sin(x7)*sin(x8)) ...
                              + ez * 2 * (cos(x7)*-sin(x8) - sin(x7)*cos(x8));
        end

        function dUzddottopdx = dUzddottopdx(model, state, ref)
            dUzddottopdx  = zeros(1, model.numStates);

            dUzdottopdot_dx = derivatives_drone_IB.dUzdottopdotdx(model, state, ref);
            dUzdottop_dx = derivatives_drone_IB.dUzdottopdx(model, state, ref);
            
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

            dUzddottopdx(3) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(3) - dUzdottop_dx(3)*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            dUzddottopdx(6) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(6) - dUzdottop_dx(6)*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            dUzddottopdx(10) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(10);
            dUzddottopdx(11) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(11);
            dUzddottopdx(16) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(16) - dUzdottop_dx(16)*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;

            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            ez_dot = (1-model.c7^2+model.lambda4)*(ref.z_dot-x6)+(model.c7+model.c8)*e8_dot-model.c7*model.lambda4*e7;
            
           
            z_dddot = 4/m*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            e8_ddot = c7*(ref.z_ddot-x6)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);

            Uzdottopdot = cos(x7)*cos(x8)*ez_ddot + 2*(cos(x7)*cos(x8)-sin(x7)*sin(x8))*ez;
            Uzdottop = cos(x7)*cos(x8)*ez_dot + ez*(sin(x7)*cos(x8)+cos(x7)*sin(x8));

            dUzddottopdx(7) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(7) + 2*cos(x7)*cos(x8)*-sin(x7)*Uzdottopdot ...
                              -dUzdottop_dx(7)*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)) ...
                              -Uzdottop*2*(-sin(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)) ...
                              -Uzdottop*2*(cos(x7)*cos(x8))*(-cos(x7)*cos(x8)+sin(x7)*sin(x8));

            dUzddottopdx(8) = (cos(x7)*cos(x8))^2*dUzdottopdot_dx(8) + 2*cos(x7)*cos(x8)*-sin(x8)*Uzdottopdot ...
                              -dUzdottop_dx(8)*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)) ...
                              -Uzdottop*2*(cos(x7)*-sin(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)) ...
                              -Uzdottop*2*(cos(x7)*cos(x8))*(-sin(x7)*-sin(x8)-cos(x7)*cos(x8));
        end

        function dUzddotdx = dUzddotdx(model, state, ref)
            dUzddotdx = zeros(1, model.numStates);

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

            denom = cos(x7)*cos(x8)^4;

            dUzddottop_dx = derivatives_drone_IB.dUzddottopdx(model, state, ref);
            dUzddotdx(3) = m/4*dUzddottop_dx(3)/denom;
            dUzddotdx(6) = m/4*dUzddottop_dx(6)/denom;
            dUzddotdx(10) = m/4*dUzddottop_dx(10)/denom;
            dUzddotdx(11) = m/4*dUzddottop_dx(11)/denom;
            dUzddotdx(16) = m/4*dUzddottop_dx(16)/denom;
            
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;

            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            ez_dot = (1-model.c7^2+model.lambda4)*(ref.z_dot-x6)+(model.c7+model.c8)*e8_dot-model.c7*model.lambda4*e7;
            
            
            z_dddot = 4/m*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            e8_ddot = c7*(ref.z_ddot-x6)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);

            Uzdottopdot = cos(x7)*cos(x8)*ez_ddot + 2*(cos(x7)*cos(x8)-sin(x7)*sin(x8))*ez;
            Uzdottop = cos(x7)*cos(x8)*ez_dot + ez*(sin(x7)*cos(x8)+cos(x7)*sin(x8));
            Uzddottop = (cos(x7)*cos(x8))^2*Uzdottopdot - Uzdottop*2*(cos(x7)*cos(x8))*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));

            dUzddotdx(7) = m/4*(dUzddottop_dx(7)*denom-Uzddottop*4*(cos(x7)*cos(x8))^3*-sin(x7))/denom^2;
            dUzddotdx(8) = m/4*(dUzddottop_dx(8)*denom-Uzddottop*4*(cos(x7)*cos(x8))^3*-sin(x8))/denom^2;
        end

        function dkomhdx = dkomhdx(model, state, ref)
            dkomhdx = zeros(1, model.numStates);

            m = model.m;
            lambda5 = model.lambda5;
            c9 = model.c9;
            c10 = model.c10;
            x1 = state(1);
            x4 = state(4);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x17 = state(17);


            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            x_ddot = 4*U_z*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(m);
            e9 = ref.x - x1;
          
            e10 = c9*e9 + ref.x_dot + lambda5*x17-x4;
            e10_dot= c9*(ref.x_dot - x4) - ref.x_ddot + lambda5*e9 - x_ddot;
            ex_dot = (1-c9^2+lambda5)*(ref.x_dot-x4)+(c9+c10)*e10_dot+ref.x_dddot-c9*lambda5*e9;
            ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + ref.x_ddot - c9*lambda5*x17;
            
            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);
            dex_dx = derivatives_drone_IB.dexdx(model, state, ref);
            dexdot_dx = derivatives_drone_IB.dexdotdx(model, state, ref);  

            dkomhdx(1) = U_z*dexdot_dx(1) - U_zdot*dex_dx(1);
            dkomhdx(4) = U_z*dexdot_dx(4) - U_zdot*dex_dx(4);
            dkomhdx(17) =  - U_zdot*dex_dx(17);

            dkomhdx(3) = dUz_dx(3)*ex_dot - ex*dUzdot_dx(3);
            dkomhdx(6) = dUz_dx(6)*ex_dot - ex*dUzdot_dx(6);
            dkomhdx(7) = dUz_dx(7)*ex_dot - ex*dUzdot_dx(7);
            dkomhdx(8) = dUz_dx(8)*ex_dot - ex*dUzdot_dx(8);
            dkomhdx(16) = dUz_dx(16)*ex_dot - ex*dUzdot_dx(16);
        end

        function dkomhdotdx = dkomhdotdx(model, state, ref)
            dkomhdotdx = zeros(1, model.numStates);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;
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

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            e9 = ref.x - x1;
            e10 = c9*e9 + ref.x_dot + lambda5*x17-x4;
            ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + ref.x_ddot - c9*lambda5*x17;
            
            x_dddot = 4/model.m*( ...
                U_zdot*(sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9)) + ...
                U_z*(cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 + ...
                     cos(x7)*x10*cos(x8)*cos(x9) - ...
                     sin(x7)*sin(x8)*x11*cos(x9) - ...
                     sin(x7)*cos(x8)*sin(x9)*x12));

            e10_ddot = c9*(ref.x_ddot - x4) - ref.x_dddot + lambda5*(ref.x_dot - x4) - x_dddot;
            ex_ddot = (1-c9^2+lambda5)*(ref.x_ddot - x4) + (c9+c10)*e10_ddot + ref.x_ddddot - c9*lambda5*(ref.x_dot - x4);
           
            z_dddot = 4/(m)*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;   
            e8_ddot = c7*(ref.z_ddot-z_ddot)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);
            U_zdot_top_dot = -sin(x7)*cos(x8)*ez_dot - cos(x7)*sin(x8)*ez_dot + cos(x7)*cos(x8)*ez_ddot + ...
                cos(x7)*cos(x8)*ez - sin(x7)*sin(x8)*ez + sin(x7)*cos(x8)*ez_dot - ...
                sin(x7)*sin(x8)*ez + cos(x7)*cos(x8)*ez + cos(x7)*sin(x8)*ez_dot;
            U_zdot_top = cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            U_zddot = m/4*(cos(x7)*cos(x8)^2*U_zdot_top_dot - U_zdot_top*2*cos(x7)*cos(x8)*(-sin(x7*cos(x8)-cos(x7)*sin(x8))))/(cos(x7)*cos(x8))^4;

            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dex_dx = derivatives_drone_IB.dexdx(model, state, ref);
            dUzddot_dx = derivatives_drone_IB.dUzddotdx(model, state, ref);
            dexddot_dx = derivatives_drone_IB.dexddotdx(model, state, ref);
            
            dkomhdotdx(1) = U_zddot*dex_dx(1);
            dkomhdotdx(4) = U_zddot*dex_dx(4);
            dkomhdotdx(17) = U_zddot*dex_dx(17);

            dkomhdotdx(3) = dUz_dx(3)*ex_ddot + dexddot_dx(3)*U_z - ex*dUzddot_dx(3);
            dkomhdotdx(6) = dUz_dx(6)*ex_ddot + dexddot_dx(6)*U_z - ex*dUzddot_dx(6);
            dkomhdotdx(7) = dUz_dx(7)*ex_ddot + dexddot_dx(7)*U_z - ex*dUzddot_dx(7);
            dkomhdotdx(8) = dUz_dx(8)*ex_ddot + dexddot_dx(8)*U_z - ex*dUzddot_dx(8);

            dkomhdotdx(9) = dexddot_dx(9)*U_z ;
            dkomhdotdx(10) =  dexddot_dx(10)*U_z - ex*dUzddot_dx(10);
            dkomhdotdx(11) =  dexddot_dx(11)*U_z - ex*dUzddot_dx(11);
            dkomhdotdx(12) =  dexddot_dx(12)*U_z ;

            dkomhdotdx(16) = dUz_dx(16)*ex_ddot + dexddot_dx(16)*U_z - ex*dUzddot_dx(16);
        end

        function dciwywdx = dciwywdx(model, state, ref)
            dciwywdx = zeros(1, model.numStates);

            m = model.m;
            lambda6 = model.lambda6;
            x2 = state(2);
            x5 = state(5);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x18 = state(18);

            c11 = model.c11;
            c12 = model.c12;
        

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);

            y_ddot = 4*U_z*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9))/(m);
            e11 = ref.y - x2;
            e12 = c11*e11 + ref.y_dot + model.lambda5*x18 - x5;
            e12_dot= c11*(ref.y_dot - x5) - ref.y_ddot + lambda6*e11 - y_ddot;
            ey     = -(1-c11^2+lambda6)*e11 - (c11+c12)*e12 - ref.y_ddot + c11*lambda6*x18;
            ey_dot = -(1-c11^2+lambda6)*(ref.y_dot-x5)-(c11+c12)*e12_dot-ref.y_dddot+c11*lambda6*e11;
            
            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);
            dey_dx = derivatives_drone_IB.deydx(model, state, ref);
            deydot_dx = derivatives_drone_IB.deydotdx(model, state, ref);  

            dciwywdx(1) = U_z*deydot_dx(1) - U_zdot*dey_dx(1);
            dciwywdx(4) = U_z*deydot_dx(4) - U_zdot*dey_dx(4);
            dciwywdx(17) =  - U_zdot*dey_dx(17);

            dciwywdx(3) = dUz_dx(3)*ey_dot - ey*dUzdot_dx(3);
            dciwywdx(6) = dUz_dx(6)*ey_dot - ey*dUzdot_dx(6);
            dciwywdx(7) = dUz_dx(7)*ey_dot - ey*dUzdot_dx(7);
            dciwywdx(8) = dUz_dx(8)*ey_dot - ey*dUzdot_dx(8);
            dciwywdx(16) = dUz_dx(16)*ey_dot - ey*dUzdot_dx(16);
        end

        function dciwywdotdx = dciwywdotdx(model, state, ref)
            dciwywdotdx = zeros(1, model.numStates);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            lambda4 = model.lambda4;
            lambda6 = model.lambda6;
         
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

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            e11 = ref.y - x2;
            e12 = c11*e11 + ref.y_dot + model.lambda5*x18 - x5;
            ey = (1-c11^2+lambda6)*e11 + (c11+c12)*e12 + ref.y_ddot - c11*lambda6*x18;
            
            y_dddot = 4/model.m*(U_zdot*(-sin(x8)*cos(x9)+sin(x7)*cos(x8)*sin(x9)) ...
                    + U_z *(-cos(x8)*x11*cos(x9)+sin(x8)*sin(x9)*x12 ...
                            +cos(x7)*x10*cos(x9)*sin(x9) ...
                            -sin(x7)*sin(x8)*x11*sin(x9) ...
                            +sin(x7)*cos(x8)*cos(x9)*x12));

            e12_ddot = c11*(ref.y_ddot - x5) - ref.y_dddot + lambda6*(ref.y_dot - x5) - y_dddot;
            ey_ddot = -(1-c11^2+lambda6)*(ref.x_ddot - x5) - (c11+c12)*e12_ddot - ref.x_ddddot + c11*lambda6*(ref.x_dot - x5);
           
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            z_dddot = 4/(m)*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            e8_ddot = c7*(ref.z_ddot-z_ddot)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;

            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;

            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);
            U_zdot_top_dot = -sin(x7)*cos(x8)*ez_dot - cos(x7)*sin(x8)*ez_dot + cos(x7)*cos(x8)*ez_ddot + ...
                cos(x7)*cos(x8)*ez - sin(x7)*sin(x8)*ez + sin(x7)*cos(x8)*ez_dot - ...
                sin(x7)*sin(x8)*ez + cos(x7)*cos(x8)*ez + cos(x7)*sin(x8)*ez_dot;
            U_zdot_top = cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            U_zddot = m/4*(cos(x7)*cos(x8)^2*U_zdot_top_dot - U_zdot_top*2*cos(x7)*cos(x8)*(-sin(x7*cos(x8)-cos(x7)*sin(x8))))/(cos(x7)*cos(x8))^4;

            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dey_dx = derivatives_drone_IB.deydx(model, state, ref);
            dUzddot_dx = derivatives_drone_IB.dUzddotdx(model, state, ref);
            deyddot_dx = derivatives_drone_IB.deyddotdx(model, state, ref);
            
            dciwywdotdx(1) = U_zddot*dey_dx(1);
            dciwywdotdx(4) = U_zddot*dey_dx(4);
            dciwywdotdx(17) = U_zddot*dey_dx(17);

            dciwywdotdx(3) = dUz_dx(3)*ey_ddot + deyddot_dx(3)*U_z - ey*dUzddot_dx(3);
            dciwywdotdx(6) = dUz_dx(6)*ey_ddot + deyddot_dx(6)*U_z - ey*dUzddot_dx(6);
            dciwywdotdx(7) = dUz_dx(7)*ey_ddot + deyddot_dx(7)*U_z - ey*dUzddot_dx(7);
            dciwywdotdx(8) = dUz_dx(8)*ey_ddot + deyddot_dx(8)*U_z - ey*dUzddot_dx(8);

            dciwywdotdx(9) = deyddot_dx(9)*U_z ;
            dciwywdotdx(10) =  deyddot_dx(10)*U_z - ey*dUzddot_dx(10);
            dciwywdotdx(11) =  deyddot_dx(11)*U_z - ey*dUzddot_dx(11);
            dciwywdotdx(12) =  deyddot_dx(12)*U_z ;

            dciwywdotdx(16) = dUz_dx(16)*ey_ddot + deyddot_dx(16)*U_z - ey*dUzddot_dx(16);
        end

        function dthetarefdx = dthetarefdx(model, state, ref)
            dthetarefdx = zeros(1, model.numStates);
            x1 = state(1);
            x4 = state(4);
            x17 = state(17);

            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;

            e9 = ref.x - x1;
            e10 = c9*e9 + ref.x_dot + lambda5*x17-x4;
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);

            dthetarefdx(1) = m/U_z*((1-c9^2+lambda5)*(-1) + (c9+c10)*(-c9));
            dthetarefdx(4) = m/U_z*((c9+c10)*(-1));
            dthetarefdx(17) = m/U_z*((c9+c10)*lambda5-c9*lambda5);

            thetaRef_U_z_coeff = m * ((1-c9^2+lambda5)*e9+(c9+c10)*e10+ref.x_ddot-c9*lambda5*x17);
            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);
            dthetarefdx(3) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(3);
            dthetarefdx(6) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(6);
            dthetarefdx(7) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(7);
            dthetarefdx(8) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(8);
            dthetarefdx(16) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(16);
        end

        function dthetarefdotdx = dthetarefdotdx(model, state, ref)
            dthetarefdotdx = zeros(1, model.numStates);
            
            x1 = state(1);
            x4 = state(4);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x17 = state(17);

            m = model.m;
            c9 = model.c9;
            c10 = model.c10;
            lambda5 = model.lambda5;

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            Uz_dot = derivatives_drone_IB.Uzdot(model, state, ref);

            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);
            de10_dx = derivatives_drone_IB.de10dx(model);
            de10dot_dx = derivatives_drone_IB.de10dotdx(model, state, ref);

            dthetarefdotdx(1)  = m/U_z^2*(U_z*((c9+c10)*de10dot_dx(1) + c9*lambda5) - Uz_dot*(-(1-c9^2+lambda5) + (c9+c10)*de10_dx(1)));
            dthetarefdotdx(4)  = m/U_z^2*(U_z*(-(1-c9^2+lambda5) + (c9+c10)*de10dot_dx(4)) - Uz_dot*((c9+c10)*de10_dx(4)));
            dthetarefdotdx(17) = m/U_z^2*(Uz_dot*((c9+c10)*de10_dx(17) - c9*lambda5));
            
            dthetarefdotdx(9)  = m/U_z^2*(U_z*(c9+c10)*de10dot_dx(9));

            e9 = ref.x - x1;
            e10 = c9*e9 + ref.x_dot + model.lambda5*x17 - x4;
            x_ddot = 4*U_z*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(model.m);
            e10_dot = c9*(ref.x_dot - x4) - ref.x_ddot + lambda5*e9 - x_ddot;
            
            Uzcoeff = (1-c9^2+lambda5)*(ref.x_dot-x4)+(c9+c10)*e10_dot+ref.x_dddot-c9*lambda5*(ref.x-x1);
            Uzdotcoeff = (1-c9^2+lambda5)*(ref.x-x1)+(c9+c10)*e10+ref.x_ddot-c9*lambda5*x17;
            thetaRefNum = m*(U_z*((1-c9^2+lambda5)*(ref.x_dot - x4)+(c9+c10)*e10_dot+ref.x_dddot - c9*lambda5*(ref.x-x1)) ...
                            -Uz_dot*((1-c9^2+lambda5)*(ref.x-x1))+(c9+c10)*e10+ref.x_ddot-c9*lambda5*x17 );
            
            dthetarefdotdx(3)  = m*(dUz_dx(3)*Uzcoeff+U_z*(c9+c10)*de10dot_dx(3)-dUzdot_dx(3)*Uzdotcoeff)-thetaRefNum*2*U_z*dUz_dx(3)/U_z^4;
            dthetarefdotdx(6)  = m*(dUz_dx(6)*Uzcoeff+U_z*(c9+c10)*de10dot_dx(6)-dUzdot_dx(6)*Uzdotcoeff)-thetaRefNum*2*U_z*dUz_dx(6)/U_z^4;
            dthetarefdotdx(7)  = m*(dUz_dx(7)*Uzcoeff+U_z*(c9+c10)*de10dot_dx(7)-dUzdot_dx(7)*Uzdotcoeff)-thetaRefNum*2*U_z*dUz_dx(7)/U_z^4;
            dthetarefdotdx(8)  = m*(dUz_dx(8)*Uzcoeff+U_z*(c9+c10)*de10dot_dx(8)-dUzdot_dx(8)*Uzdotcoeff)-thetaRefNum*2*U_z*dUz_dx(8)/U_z^4;
            dthetarefdotdx(16)  = m*(dUz_dx(16)*Uzcoeff+U_z*(c9+c10)*de10dot_dx(16)-dUzdot_dx(16)*Uzdotcoeff)-thetaRefNum*2*U_z*dUz_dx(16)/U_z^4;
        end

        function dthetarefddotdx = dthetarefddotdx(model, state, ref)
            dthetarefddotdx = zeros(1, model.numStates);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;

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

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);

            e9 = ref.x - x1;
            e10 = c9*e9 + ref.x_dot + lambda5*x17 - x4;
            x_ddot = 4*U_z*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(m);
            e10_dot= c9*(ref.x_dot - x4) - ref.x_ddot + lambda5*e9 - x_ddot;
            ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + ref.x_ddot - c9*lambda5*x17;
            ex_dot = (1-c9^2+lambda5)*(ref.x_dot-x4)+(c9+c10)*e10_dot+ref.x_dddot-c9*lambda5*e9;
            x_dddot = 4/model.m*( ...
                U_zdot*(sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9)) + ...
                U_z*(cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 + ...
                     cos(x7)*x10*cos(x8)*cos(x9) - ...
                     sin(x7)*sin(x8)*x11*cos(x9) - ...
                     sin(x7)*cos(x8)*sin(x9)*x12));

            e10_ddot = c9*(ref.x_ddot - x_ddot) - ref.x_dddot + lambda5*(ref.x_dot - x4) - x_dddot;
            ex_ddot = -(1-c9^2+lambda5)*(ref.x_ddot - x_ddot) - (c9+c10)*e10_ddot - ref.x_ddddot + c9*lambda5*(ref.x_dot - x4);
            
            z_dddot = 4/(m)*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;            
            e8_ddot = c7*(ref.z_ddot-z_ddot)+ref.z_dddot+lambda4*(ref.z_dot-x4)-z_dddot;

            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;

            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);
            U_zdot_top_dot = -sin(x7)*cos(x8)*ez_dot - cos(x7)*sin(x8)*ez_dot + cos(x7)*cos(x8)*ez_ddot + ...
                cos(x7)*cos(x8)*ez - sin(x7)*sin(x8)*ez + sin(x7)*cos(x8)*ez_dot - ...
                sin(x7)*sin(x8)*ez + cos(x7)*cos(x8)*ez + cos(x7)*sin(x8)*ez_dot;
            U_zdot_top = cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            U_zddot = m/4*(cos(x7)*cos(x8)^2*U_zdot_top_dot - U_zdot_top*2*cos(x7)*cos(x8)*(-sin(x7*cos(x8)-cos(x7)*sin(x8))))/(cos(x7)*cos(x8))^4;

            komh =  U_z*ex_dot-ex*U_zdot;
            komh_dot = U_z*ex_ddot-ex*U_zddot;

            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);
            dU_zdotdx = derivatives_drone_IB.dUzdotdx(model, state, ref);
            dkomh_dx = derivatives_drone_IB.dkomhdx(model, state, ref);
            dkomhdot_dx = derivatives_drone_IB.dkomhdotdx(model, state, ref);
                
            thetaddotnum   = U_z^2*komh_dot-komh*2*U_z*U_zdot;
            thetaddotdenom = U_z^4;
            dthetarefddotdx(1) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(1) - dkomh_dx(1)*2*U_z*U_zdot);
            dthetarefddotdx(4) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(4) - dkomh_dx(4)*2*U_z*U_zdot);
            dthetarefddotdx(17) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(17) - dkomh_dx(17)*2*U_z*U_zdot);

            dthetarefddotdx(9) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(10) ) ;
            dthetarefddotdx(10) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(10) ) ;
            dthetarefddotdx(11) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(11) ) ;
            dthetarefddotdx(12) = m/thetaddotdenom*(U_z^2*dkomhdot_dx(12) ) ;

            dthetarefddotdx(3) = m/thetaddotdenom^2*( ...
                (2*U_z*dU_zdx(3)*komh_dot + U_z^2*dkomhdot_dx(3) ...
                - dkomh_dx(3)*2*U_z*U_zdot - komh*2*dU_zdx(3)*U_zdot -komh*2*U_z*dU_zdotdx(3))*thetaddotdenom ...
                - thetaddotnum * 4*U_z^3*dU_zdx(3));

            dthetarefddotdx(6) = m/thetaddotdenom^2*( ...
                (2*U_z*dU_zdx(6)*komh_dot + U_z^2*dkomhdot_dx(6) ...
                - dkomh_dx(6)*2*U_z*U_zdot - komh*2*dU_zdx(6)*U_zdot -komh*2*U_z*dU_zdotdx(6))*thetaddotdenom ...
                - thetaddotnum * 4*U_z^3*dU_zdx(6));

            dthetarefddotdx(7) = m/thetaddotdenom^2*( ...
                (2*U_z*dU_zdx(7)*komh_dot + U_z^2*dkomhdot_dx(7) ...
                - dkomh_dx(7)*2*U_z*U_zdot - komh*2*dU_zdx(7)*U_zdot -komh*2*U_z*dU_zdotdx(7))*thetaddotdenom ...
                - thetaddotnum * 4*U_z^3*dU_zdx(7));

            dthetarefddotdx(8) = m/thetaddotdenom^2*( ...
                (2*U_z*dU_zdx(8)*komh_dot + U_z^2*dkomhdot_dx(8) ...
                - dkomh_dx(8)*2*U_z*U_zdot - komh*2*dU_zdx(8)*U_zdot -komh*2*U_z*dU_zdotdx(8))*thetaddotdenom ...
                - thetaddotnum * 4*U_z^3*dU_zdx(8));

            dthetarefddotdx(16) = m/thetaddotdenom^2*( ...
                (2*U_z*dU_zdx(16)*komh_dot + U_z^2*dkomhdot_dx(16) ...
                - dkomh_dx(16)*2*U_z*U_zdot - komh*2*dU_zdx(16)*U_zdot -komh*2*U_z*dU_zdotdx(16))*thetaddotdenom ...
                - thetaddotnum * 4*U_z^3*dU_zdx(16));

        end

        function dphirefdx = dphirefdx(model, state, ref)
            dphirefdx = zeros(1, model.numStates);
            x2 = state(2);
            x5 = state(5);
            x17 = state(17);
            x18 = state(18);

            m = model.m;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;

            e9 = ref.y - x2;
            e10 = c11*e9 + ref.y_dot + lambda6*x17-x5;
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);

            dphirefdx(2) = m/U_z*(-(1-c11^2+lambda6)*(-1) - (c11+c12)*(-c11));
            dphirefdx(5) = m/U_z*(-(c11+c12)*(-1));
            dphirefdx(18) = m/U_z*(-(c11+c12)*lambda6+c11*lambda6);

            phiRef_U_z_coeff = m * ((1-c11^2+lambda6)*e9+(c11+c12)*e10+ref.y_ddot-c11*lambda6*x18);
            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);
            dphirefdx(3) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(3);
            dphirefdx(6) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(6);
            dphirefdx(7) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(7);
            dphirefdx(8) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(8);
            dphirefdx(16) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(16);
        end

        function dphirefdotdx = dphirefdotdx(model, state, ref)
            dphirefdotdx = zeros(1, model.numStates);
            
            x2 = state(2);
            x5 = state(5);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x18 = state(18);

            m = model.m;
            c11 = model.c11;
            c12 = model.c12;
            lambda6 = model.lambda6;

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            Uz_dot =  derivatives_drone_IB.Uzdot(model, state, ref);

            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx =  derivatives_drone_IB.dUzdotdx(model, state, ref);
            de12_dx =  derivatives_drone_IB.de12dx(model);
            de12dot_dx =  derivatives_drone_IB.de12dotdx(model, state, ref);

            dphirefdotdx(2)  = m/U_z^2*(U_z*((c11+c12)*de12dot_dx(2) + c11*lambda6) - Uz_dot*(-(1-c11^2+lambda6) + (c11+c12)*de12_dx(2)));
            dphirefdotdx(5)  = m/U_z^2*(U_z*(-(1-c11^2+lambda6) + (c11+c12)*de12dot_dx(5)) - Uz_dot*((c11+c12)*de12_dx(5)));
            dphirefdotdx(18) = m/U_z^2*(Uz_dot*((c11+c12)*de12_dx(18) - c11*lambda6));
            
            dphirefdotdx(9)  = m/U_z^2*(U_z*(c11+c12)*de12dot_dx(9));

            e11 = ref.y - x2;
            e12 = c11*e11 + ref.y_dot + model.lambda5*x18 - x5;
            y_ddot = 4*U_z*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9))/(model.m);
            e12_dot = c11*(ref.y_dot - x5) - ref.y_ddot + lambda6*e11 - y_ddot;
            
            Uzcoeff = -(1-c11^2+lambda6)*(ref.y_dot-x5)-(c11+c12)*e12_dot-ref.y_dddot+c11*lambda6*(ref.y-x2);
            Uzdotcoeff = -(1-c11^2+lambda6)*(ref.y-x2)-(c11+c12)*e12-ref.y_ddot+c11*lambda6*x18;
            phiRefNum = m*(U_z*(-(1-c11^2+lambda6)*(ref.y_dot - x5)-(c11+c12)*e12_dot-ref.y_dddot + c11*lambda6*(ref.y-x2)) ...
                            -Uz_dot*(-(1-c11^2+lambda6)*(ref.y-x2))-(c11+c12)*e12-ref.y_ddot+c11*lambda6*x18 );
            
            dphirefdotdx(3)  = m*(dUz_dx(3)*Uzcoeff-U_z*(c11+c12)*de12dot_dx(3)-dUzdot_dx(3)*Uzdotcoeff)-phiRefNum*2*U_z*dUz_dx(3)/U_z^4;
            dphirefdotdx(6)  = m*(dUz_dx(6)*Uzcoeff-U_z*(c11+c12)*de12dot_dx(6)-dUzdot_dx(6)*Uzdotcoeff)-phiRefNum*2*U_z*dUz_dx(6)/U_z^4;
            dphirefdotdx(7)  = m*(dUz_dx(7)*Uzcoeff-U_z*(c11+c12)*de12dot_dx(7)-dUzdot_dx(7)*Uzdotcoeff)-phiRefNum*2*U_z*dUz_dx(7)/U_z^4;
            dphirefdotdx(8)  = m*(dUz_dx(8)*Uzcoeff-U_z*(c11+c12)*de12dot_dx(8)-dUzdot_dx(8)*Uzdotcoeff)-phiRefNum*2*U_z*dUz_dx(8)/U_z^4;
            dphirefdotdx(16)  = m*(dUz_dx(16)*Uzcoeff-U_z*(c11+c12)*de12dot_dx(16)-dUzdot_dx(16)*Uzdotcoeff)-phiRefNum*2*U_z*dUz_dx(16)/U_z^4;
        end

        function dphirefddotdx = dphirefddotdx(model, state, ref)
            dphirefddotdx = zeros(1, model.numStates);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            lambda4 = model.lambda4;
            lambda6 = model.lambda6;

            x3 = state(3);
            x2 = state(2);
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

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);

            e11 = ref.y - x2;
            e12 = c11*e11 + ref.y_dot + lambda6*x18 - x5; 
            
            y_ddot = 4*U_z*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9))/(model.m);
            e12_dot= c11*(ref.y_dot - x5) - ref.y_ddot + lambda6*e11 - y_ddot;
            
            ey     = -(1-c11^2+lambda6)*e11 - (c11+c12)*e12 - ref.y_ddot + c11*lambda6*x18;
            ey_dot = -(1-c11^2+lambda6)*(ref.y_dot-x5)-(c11+c12)*e12_dot-ref.y_dddot+c11*lambda6*e11;

            y_dddot = 4/model.m*(U_zdot*(-sin(x8)*cos(x9)+sin(x7)*cos(x8)*sin(x9)) ...
                    + U_z *(-cos(x8)*x11*cos(x9)+sin(x8)*sin(x9)*x12 ...
                            +cos(x7)*x10*cos(x9)*sin(x9) ...
                            -sin(x7)*sin(x8)*x11*sin(x9) ...
                            +sin(x7)*cos(x8)*cos(x9)*x12));

            e12_ddot = c11*(ref.y_ddot-y_ddot) - ref.y_dddot + lambda6*(ref.y_dot-x5)-y_dddot;
            ey_ddot = (1-c11^2+lambda6)*(ref.y_ddot-y_ddot)+(c11+c12)*e12_ddot+ref.y_ddddot-c11*lambda6*(ref.y_dot-x5);
            
           
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            z_dddot = 4/(m)*(U_zdot*cos(x8)*cos(x7)-U_z*sin(x8)*x11*cos(x7)-U_z*cos(x7)*sin(x7)*x10);
            e8_ddot = c7*(ref.z_ddot-z_ddot)+ref.z_dddot+lambda4*(ref.z_dot-x6)-z_dddot;
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;

            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez_ddot = (1-c7^2+lambda4)*(ref.z_ddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(ref.z_dot-x6);
            U_zdot_top_dot = -sin(x7)*cos(x8)*ez_dot - cos(x7)*sin(x8)*ez_dot + cos(x7)*cos(x8)*ez_ddot + ...
                cos(x7)*cos(x8)*ez - sin(x7)*sin(x8)*ez + sin(x7)*cos(x8)*ez_dot - ...
                sin(x7)*sin(x8)*ez + cos(x7)*cos(x8)*ez + cos(x7)*sin(x8)*ez_dot;
            U_zdot_top = cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            U_zddot = m/4*(cos(x7)*cos(x8)^2*U_zdot_top_dot - U_zdot_top*2*cos(x7)*cos(x8)*(-sin(x7*cos(x8)-cos(x7)*sin(x8))))/(cos(x7)*cos(x8))^4;

            ciwyw =  U_z*ey_dot-ey*U_zdot;
            ciwyw_dot = U_z*ey_ddot-ey*U_zddot;

            dU_zdx = derivatives_drone_IB.dUzdx(model, state, ref);
            dU_zdotdx = derivatives_drone_IB.dUzdotdx(model, state, ref);
            dciwyw_dx = derivatives_drone_IB.dciwywdx(model, state, ref);
            dciwywdot_dx = derivatives_drone_IB.dciwywdotdx(model, state, ref);
                
            phiddotnum   = U_z^2*ciwyw_dot-ciwyw*2*U_z*U_zdot;
            phiddotdenom = U_z^4;
            dphirefddotdx(1) = m/phiddotdenom*(U_z^2*dciwywdot_dx(1) - dciwyw_dx(1)*2*U_z*U_zdot);
            dphirefddotdx(4) = m/phiddotdenom*(U_z^2*dciwywdot_dx(4) - dciwyw_dx(4)*2*U_z*U_zdot);
            dphirefddotdx(17) = m/phiddotdenom*(U_z^2*dciwywdot_dx(17) - dciwyw_dx(17)*2*U_z*U_zdot);

            dphirefddotdx(9) = m/phiddotdenom*(U_z^2*dciwywdot_dx(10) ) ;
            dphirefddotdx(10) = m/phiddotdenom*(U_z^2*dciwywdot_dx(10) ) ;
            dphirefddotdx(11) = m/phiddotdenom*(U_z^2*dciwywdot_dx(11) ) ;
            dphirefddotdx(12) = m/phiddotdenom*(U_z^2*dciwywdot_dx(12) ) ;

            dphirefddotdx(3) = m/phiddotdenom^2*( ...
                (2*U_z*dU_zdx(3)*ciwyw_dot + U_z^2*dciwywdot_dx(3) ...
                - dciwyw_dx(3)*2*U_z*U_zdot - ciwyw*2*dU_zdx(3)*U_zdot -ciwyw*2*U_z*dU_zdotdx(3))*phiddotdenom ...
                - phiddotnum * 4*U_z^3*dU_zdx(3));

            dphirefddotdx(6) = m/phiddotdenom^2*( ...
                (2*U_z*dU_zdx(6)*ciwyw_dot + U_z^2*dciwywdot_dx(6) ...
                - dciwyw_dx(6)*2*U_z*U_zdot - ciwyw*2*dU_zdx(6)*U_zdot -ciwyw*2*U_z*dU_zdotdx(6))*phiddotdenom ...
                - phiddotnum * 4*U_z^3*dU_zdx(6));

            dphirefddotdx(7) = m/phiddotdenom^2*( ...
                (2*U_z*dU_zdx(7)*ciwyw_dot + U_z^2*dciwywdot_dx(7) ...
                - dciwyw_dx(7)*2*U_z*U_zdot - ciwyw*2*dU_zdx(7)*U_zdot -ciwyw*2*U_z*dU_zdotdx(7))*phiddotdenom ...
                - phiddotnum * 4*U_z^3*dU_zdx(7));

            dphirefddotdx(8) = m/phiddotdenom^2*( ...
                (2*U_z*dU_zdx(8)*ciwyw_dot + U_z^2*dciwywdot_dx(8) ...
                - dciwyw_dx(8)*2*U_z*U_zdot - ciwyw*2*dU_zdx(8)*U_zdot -ciwyw*2*U_z*dU_zdotdx(8))*phiddotdenom ...
                - phiddotnum * 4*U_z^3*dU_zdx(8));

            dphirefddotdx(16) = m/phiddotdenom^2*( ...
                (2*U_z*dU_zdx(16)*ciwyw_dot + U_z^2*dciwywdot_dx(16) ...
                - dciwyw_dx(16)*2*U_z*U_zdot - ciwyw*2*dU_zdx(16)*U_zdot -ciwyw*2*U_z*dU_zdotdx(16))*phiddotdenom ...
                - phiddotnum * 4*U_z^3*dU_zdx(16));

        end

        function dxdddotdx = dxdddotdx(model, state, ref)
            dxdddotdx = zeros(1, model.numStates);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
           
            m = model.m;
        
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            Uz_dot = derivatives_drone_IB.Uzdot(model, state, ref);
            
            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);

            Uzdot_coeff = sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9);
            Uz_coeff = cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 ...
                     + cos(x7)*cos(x8)*cos(x9)*x10 ...
                     - sin(x7)*sin(x8)*cos(x9)*x11 ...
                     - sin(x7)*cos(x8)*sin(x9)*x12;

            dxdddotdx(3) = 4/m * (dUzdot_dx(3)*Uzdot_coeff + dUz_dx(3)*Uz_coeff);
            dxdddotdx(6) = 4/m * (dUzdot_dx(6)*Uzdot_coeff + dUz_dx(6)*Uz_coeff);
            dxdddotdx(16) = 4/m * (dUzdot_dx(16)*Uzdot_coeff + dUz_dx(16)*Uz_coeff);

            dxdddotdx(7) = 4/m * (dUzdot_dx(7)*Uzdot_coeff + Uz_dot*cos(x7)*cos(x8)*cos(x9) ...
                               + dUz_dx(7)*Uz_coeff + U_z*( ...
                                    -sin(x7)*cos(x8)*cos(x9)*x10 ...
                                    -cos(x7)*sin(x8)*cos(x9)*x11 ...
                                    -cos(x7)*cos(x8)*sin(x9)*x12...
                               ));

            dxdddotdx(8) = 4/m * (dUzdot_dx(8)*Uzdot_coeff + Uz_dot*(cos(x8)*sin(x9) - sin(x7)*sin(x8)*cos(x9)) ...
                               + dUz_dx(8)*Uz_coeff + U_z*( ...
                                    -sin(x8)*sin(x9)*x11 + cos(x8)*cos(x9)*x12 ...
                                    -cos(x7)*sin(x8)*cos(x9)*x10 ...
                                    -sin(x7)*cos(x8)*cos(x9)*x11 ...
                                    +sin(x7)*sin(x8)*sin(x9)*x12...
                               ));

            dxdddotdx(9) = 4/m * (dUzdot_dx(9)*Uzdot_coeff + Uz_dot*(sin(x8)*cos(x9) - sin(x7)*cos(x8)*sin(x9)) ... 
                               + dUz_dx(8)*Uz_coeff + U_z*( ...
                                     cos(x8)*cos(x9)*x11 - sin(x8)*sin(x9)*x12 ...
                                    -cos(x7)*cos(x8)*sin(x9)*x10...
                                    +sin(x7)*sin(x8)*sin(x9)*x11...
                                    -sin(x7)*cos(x8)*cos(x9)*x12...
                               ));

            dxdddotdx(10) = 4/m*U_z*(cos(x7)*cos(x8)*cos(x9));
            dxdddotdx(11) = 4/m*U_z*(cos(x8)*sin(x9) - sin(x7)*sin(x8)*cos(x9));
            dxdddotdx(12) = 4/m*U_z*(sin(x8)*cos(x9) - sin(x7)*cos(x8)*sin(x9));
        end

        function dydddotdx = dydddotdx(model, state, ref)
            dydddotdx = zeros(1, model.numStates);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);

            x10 = state(10);
            x11 = state(11);
            x12 = state(12);
           
            m = model.m;
        
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            Uz_dot = derivatives_drone_IB.Uzdot(model, state, ref);
            
            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);

            Uzdot_coeff = -sin(x8)*sin(x9) + sin(x7)*cos(x8)*sin(x9);

            Uz_coeff = -cos(x8)*cos(x9)*x11 + sin(x8)*sin(x9)*x12 ...
                     + cos(x7)*cos(x8)*sin(x9)*x10 ...
                     - sin(x7)*sin(x8)*sin(x9)*x11 ...
                     + sin(x7)*cos(x8)*cos(x9)*x12;

            dydddotdx(3) = 4/m * (dUzdot_dx(3)*Uzdot_coeff + dUz_dx(3)*Uz_coeff);
            dydddotdx(6) = 4/m * (dUzdot_dx(6)*Uzdot_coeff + dUz_dx(6)*Uz_coeff);
            dydddotdx(16) = 4/m * (dUzdot_dx(16)*Uzdot_coeff + dUz_dx(16)*Uz_coeff);

            dydddotdx(7) = 4/m * (dUzdot_dx(7)*Uzdot_coeff + Uz_dot*cos(x7)*cos(x8)*sin(x9) ...
                               + dUz_dx(7)*Uz_coeff + U_z*( ...
                                    -sin(x7)*cos(x8)*sin(x9)*x10 ...
                                    -cos(x7)*sin(x8)*sin(x9)*x11 ...
                                    -cos(x7)*cos(x8)*cos(x9)*x12...
                               ));

            dydddotdx(8) = 4/m * (dUzdot_dx(8)*Uzdot_coeff + Uz_dot*(-cos(x8)*cos(x9) - sin(x7)*sin(x8)*sin(x9)) ...
                               + dUz_dx(8)*Uz_coeff + U_z*( ...
                                     sin(x8)*cos(x9)*x11 + cos(x8)*sin(x9)*x12 ...
                                    -cos(x7)*sin(x8)*sin(x9)*x10 ...
                                    -sin(x7)*cos(x8)*sin(x9)*x11 ...
                                    -sin(x7)*sin(x8)*cos(x9)*x12...
                               ));

            dydddotdx(9) = 4/m * (dUzdot_dx(9)*Uzdot_coeff + Uz_dot*(sin(x8)*sin(x9) - sin(x7)*cos(x8)*cos(x9)) ... 
                               + dUz_dx(8)*Uz_coeff + U_z*( ...
                                     cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 ...
                                    +cos(x7)*cos(x8)*cos(x9)*x10...
                                    -sin(x7)*sin(x8)*cos(x9)*x11...
                                    -sin(x7)*cos(x8)*sin(x9)*x12...
                               ));

            dydddotdx(10) = 4/m*U_z*(cos(x7)*cos(x8)*sin(x9));
            dydddotdx(11) = 4/m*U_z*(-cos(x8)*cos(x9) - sin(x7)*sin(x8)*sin(x9));
            dydddotdx(12) = 4/m*U_z*( sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9));
        end

        function dzdddotdx = dzdddotdx(model, state, ref)
            dzdddotdx = zeros(1, model.numStates);
            m = model.m;
            x7 = state(7);
            x8 = state(8);
            x10 = state(10);
            x11 = state(11);

            dUz_dx = derivatives_drone_IB.dUzdx(model, state, ref);
            dUzdot_dx = derivatives_drone_IB.dUzdotdx(model, state, ref);

            dzdddotdx(3) = 4/m*(dUzdot_dx(3)*(cos(x7)*cos(x8) - dUz_dx(3)* (cos(x7)*sin(x8)*x10 + sin(x7)*cos(x8)*x11)));
            dzdddotdx(6) = 4/m*(dUzdot_dx(6)*(cos(x7)*cos(x8) - dUz_dx(6)* (cos(x7)*sin(x8)*x10 + sin(x7)*cos(x8)*x11)));
            dzdddotdx(16) = 4/m*(dUzdot_dx(16)*(cos(x7)*cos(x8) - dUz_dx(16)* (cos(x7)*sin(x8)*x10 + sin(x7)*cos(x8)*x11)));

            U_z = derivatives_drone_IB.Uz(model, state, ref);
            U_zdot = derivatives_drone_IB.Uzdot(model, state, ref);
            dzdddotdx(7) = 4/m*(dUzdot_dx(7)*cos(x7)*cos(x8) + U_zdot*-sin(x7)*cos(x8) ...
                              -dUz_dx(7)*cos(x7)*sin(x8)*x10 - U_z*-sin(x7)*sin(x8)*x10 ...
                              -dUz_dx(7)*cos(x7)*cos(x8)*x11 - U_z*cos(x7)*cos(x8)*x11);

            dzdddotdx(8) = 4/m*(dUzdot_dx(8)*cos(x7)*cos(x8) + U_zdot*cos(x7)*-sin(x8) ...
                              -dUz_dx(8)*cos(x7)*sin(x8)*x10 - U_z*cos(x7)*cos(x8)*x10 ...
                              -dUz_dx(8)*cos(x7)*cos(x8)*x11 - U_z*sin(x7)*-sin(x8)*x11);

            dzdddotdx(10) = 4/m*(-U_z*cos(x7)*sin(x8));
            dzdddotdx(11) = 4/m*(-U_z*sin(x7)*cos(x8));
        end

        function de2dx = de2dx(model, state, ref)
            de2dx = zeros(1, model.numStates);

            c1 = model.c1;
            lambda1 = model.lambda1;
            df13_dx = derivatives_drone_IB.df13dx(model, state, ref);

            dthetarefdot_dx = derivatives_drone_IB.dthetarefdotdx(model, state, ref);
            de2dx(1) = c1*df13_dx(1) + dthetarefdot_dx(1);
            de2dx(3) = c1*df13_dx(3) + dthetarefdot_dx(3);
            de2dx(4) = c1*df13_dx(4) + dthetarefdot_dx(4);
            de2dx(6) = c1*df13_dx(6) + dthetarefdot_dx(6);
            de2dx(7) = c1*df13_dx(7) + dthetarefdot_dx(7);
            de2dx(8) = c1*df13_dx(8) + dthetarefdot_dx(8);
            de2dx(9) =  dthetarefdot_dx(9);
            de2dx(10) = -1;
            de2dx(13) = lambda1;
            de2dx(16) = c1*df13_dx(16) + dthetarefdot_dx(16);
            de2dx(17) = c1*df13_dx(17) ;
        end

        function df10dx = df10dx(model, state, ref)
            df10dx = zeros(1, model.numStates);
            c1 = model.c1;
            c2 = model.c2;
            lambda1 = model.lambda1;

            df13_dx = derivatives_drone_IB.df13dx(model, state, ref);
            de2_dx = derivatives_drone_IB.de2dx(model, state, ref);
            dethetarefddot_dx = derivatives_drone_IB.dthetarefddotdx(model, state, ref);
            df10dx(1) = (1-c1^2+lambda1)*df13_dx(1) + (c1+c2)*de2_dx(1) + dethetarefddot_dx(1);
            df10dx(3) = (1-c1^2+lambda1)*df13_dx(3) + (c1+c2)*de2_dx(3) + dethetarefddot_dx(3);
            df10dx(4) = (1-c1^2+lambda1)*df13_dx(4) + (c1+c2)*de2_dx(4) + dethetarefddot_dx(4);
            df10dx(6) = (1-c1^2+lambda1)*df13_dx(6) + (c1+c2)*de2_dx(6) + dethetarefddot_dx(6);
            df10dx(7) = (1-c1^2+lambda1)*df13_dx(7) + (c1+c2)*de2_dx(7) + dethetarefddot_dx(7);
            df10dx(8) = (1-c1^2+lambda1)*df13_dx(8) + (c1+c2)*de2_dx(8) + dethetarefddot_dx(8);
            df10dx(9) = (c1+c2)*de2_dx(9) + dethetarefddot_dx(9);
            df10dx(10) = (c1+c2)*de2_dx(10) + dethetarefddot_dx(10);
            df10dx(11) = dethetarefddot_dx(11);
            df10dx(12) = dethetarefddot_dx(12);
            df10dx(13) = (c1+c2)*de2_dx(13) + dethetarefddot_dx(13) - c1*lambda1;
            df10dx(16) = (1-c1^2+lambda1)*df13_dx(16) + (c1+c2)*de2_dx(16) + dethetarefddot_dx(16);
            df10dx(17) = (1-c1^2+lambda1)*df13_dx(17) + (c1+c2)*de2_dx(17) + dethetarefddot_dx(17);
        end

        function de4dx = de4dx(model, state, ref)
            de4dx = zeros(1, model.numStates);

            c1 = model.c1;
            lambda1 = model.lambda1;

            df14_dx = derivatives_drone_IB.df14dx(model, state, ref);
            dphirefdot_dx =  derivatives_drone_IB.dphirefdotdx(model, state, ref);
            de4dx(2) = c1*df14_dx(2) + dphirefdot_dx(2);
            de4dx(3) = c1*df14_dx(3) + dphirefdot_dx(3);
            de4dx(5) = c1*df14_dx(5) + dphirefdot_dx(5);
            de4dx(6) = c1*df14_dx(6) + dphirefdot_dx(6);
            de4dx(7) = c1*df14_dx(7) + dphirefdot_dx(7);
            de4dx(8) = c1*df14_dx(8) + dphirefdot_dx(8);
            de4dx(9) = dphirefdot_dx(9);
            de4dx(10) = -1;
            de4dx(13) = lambda1;
            de4dx(16) = c1*df14_dx(16) + dphirefdot_dx(16);
            de4dx(18) = c1*df14_dx(18) ;
        end

        function df11dx = df11dx(model, state, ref)
            df11dx = zeros(1, model.numStates);
            c1 = model.c1;
            c2 = model.c2;
            lambda1 = model.lambda1;

            df14_dx = derivatives_drone_IB.df14dx(model, state, ref);
            de4_dx = derivatives_drone_IB.de4dx(model, state, ref);
            dephirefddot_dx = derivatives_drone_IB.dphirefddotdx(model, state, ref);

            df11dx(2) = (1-c1^2+lambda1)*df14_dx(2) + (c1+c2)*de4_dx(2) + dephirefddot_dx(2);
            df11dx(3) = (1-c1^2+lambda1)*df14_dx(3) + (c1+c2)*de4_dx(3) + dephirefddot_dx(3);
            df11dx(5) = (1-c1^2+lambda1)*df14_dx(5) + (c1+c2)*de4_dx(5) + dephirefddot_dx(5);
            df11dx(6) = (1-c1^2+lambda1)*df14_dx(6) + (c1+c2)*de4_dx(6) + dephirefddot_dx(6);
            df11dx(7) = (1-c1^2+lambda1)*df14_dx(7) + (c1+c2)*de4_dx(7) + dephirefddot_dx(7);
            df11dx(8) = (1-c1^2+lambda1)*df14_dx(8) + (c1+c2)*de4_dx(8) + dephirefddot_dx(8);
            df11dx(9) = (c1+c2)*de4_dx(9) + dephirefddot_dx(9);
            df11dx(10) = (c1+c2)*de4_dx(10) + dephirefddot_dx(10);
            df11dx(11) = dephirefddot_dx(11);
            df11dx(12) = dephirefddot_dx(12);
            df11dx(13) = (c1+c2)*de4_dx(13) + dephirefddot_dx(13) - c1*lambda1;
            df11dx(16) = (1-c1^2+lambda1)*df14_dx(16) + (c1+c2)*de4_dx(16) + dephirefddot_dx(16);
            df11dx(18) = (1-c1^2+lambda1)*df14_dx(18) + (c1+c2)*de4_dx(18) + dephirefddot_dx(18);
        end

        function df12dx = df12dx(model, ~, ~)
            df12dx = zeros(1, model.numStates);

            c5 = model.c5;
            c6 = model.c6;
            lambda3 = model.lambda3;

            df12dx(3) = ((1-c5^2+lambda3)*(-1) + (c5+c6)*(-c5));
            df12dx(6) = ((c5+c6)*(-1));
            df12dx(16) = ((c5+c6)*lambda3 - c5*lambda3);
        end

        function df13dx = df13dx(model, state, ref)
            df13dx = derivatives_drone_IB.dthetarefdx(model, state, ref);
            df13dx(7) = df13dx(7) - 1;
        end

        function df14dx = df14dx(model, state, ref)
            df14dx = derivatives_drone_IB.dphirefdx(model, state, ref);
            df14dx(8) = df14dx(8) - 1;
        end

        function Uz = Uz(model, state, ref)
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);

            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            Uz = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
        end

        function Uzdot = Uzdot(model, state, ref)
            m = model.m;
            g = model.g;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);
            
            U_z = derivatives_drone_IB.Uz(model, state, ref);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            e8_dot = c7*(ref.z_dot - x6)*ref.z_ddot+lambda4*e7-z_ddot;
            ez_dot = (1-c7^2+lambda4)*(ref.z_dot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            Uzdot = m/4*(cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)))/(cos(x7)*cos(x8))^2; 
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


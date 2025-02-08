classdef derivatives_drone_IB
    methods (Static)
        function dfdx = dfdx(model, state, ref)
            dfdx = zeros(model.numStates, model.numStates);
            dfdx(1, 4) = 1;
            dfdx(2, 5) = 1;
            dfdx(3, 6) = 1;

            dfdx(4, :) = df4dx(model, state, ref);
            dfdx(5, :) = df5dx(model, state, ref);
            dfdx(6, :) = df6dx(model, state, ref);

            dfdx(7, 10) = 1;
            dfdx(8, 11) = 1;
            dfdx(9, 12) = 1;

            dfdx(10, :) = df10dx(model, state, ref);
            dfdx(11, :) = df11dx(model, state, ref);
            dfdx(12, :) = df12dx(model, state, ref);

            dfdx(13, :) = df13dx(model, state, ref);
            dfdx(14, :) = df14dx(model, state, ref);
            dfdx(15, 9) = -1;
            dfdx(16, 3) = -1;
            dfdx(17, 1) = -1;
            dfdx(18, 2) = -1;
        end

        function dUzdx = dUzdx(model, state)
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);

            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            dUzdx(3) = m/(4*cos(x7)*cos(x8))*((1-c7^2+lambda4)*(-1) + (c7+c8)*(-c7));
            dUzdx(6) = m/(4*cos(x7)*cos(x8))*((c7+c8)*(-1));
            dUzdx(7) = m/(4*cos(x8))*     (g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16) * tan(x7)*sec(x7);
            dUzdx(8) = m/(4*cos(x7))*     (g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16) * tan(x8)*sec(x8);
            dUzdx(16) = m/(4*cos(x7)*cos(x8))*((c7+c8)*lambda4 - c7*lambda4);
        end

        function df4dx = df4dx(model, state, ref)
            df4dx = zeros(1, model.numStates);
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x9 = state(9);
            x16 = state(16);

            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            f4_trig = sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9);
            dU_zdx = dUzdx(model, state);

            df4dx(3) = dU_zdx(3)*f4_trig;
            df4dx(6) = dU_zdx(6)*f4_trig;
            
            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;
            U_z = m/(4*cos(x7)*cos(x8))*(g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
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
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            f5_trig = -sin(x8)*cos(x9) + sin(x7)*cos(x8)*sin(x9);
            dU_zdx = dUzdx(model, state);

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
            x3 = state(3);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);

            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            lambda4 = model.lambda4;

            f6_trig = cos(x7)*cos(x8);
            dU_zdx = dUzdx(model, state);

            df6dx(3) = dU_zdx(3)*f6_trig;
            df6dx(6) = dU_zdx(6)*f6_trig;
            
            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;
            U_z = m/(4*cos(x7)*cos(x8))*(g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            df6dx(7) = 1/m*(dU_zdx(7)*f6_trig + U_z*(-sin(x7)*cos(x8)));
            df6dx(8) = 1/m*(dU_zdx(8)*f6_trig + U_z*(-cos(x7)*sin(x8)));

            df6dx(16) = dU_zdx(16)*f6_trig;
        end

        function dthetarefdx = dthetarefdx(model, state, ref)
            dthetarefdx = zeros(1, model.numStates);
            x1 = state(1);
            x3 = state(3);
            x4 = state(4);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);
            x17 = state(17);

            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            c9 = model.c9;
            c10 = model.c10;
            lambda4 = model.lambda4;
            lambda5 = model.lambda5;

            e9 = ref.x - x1;
            e10 = c9*e9 + ref.x_dot + lambda5*x17-x4;
            
            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;
            U_z = m/(4*cos(x7)*cos(x8))*(g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            
            dthetarefdx(1) = m/U_z*((1-c9^2+lambda5)*(-1) + (c9+c10)*(-c9));
            dthetarefdx(4) = m/U_z*((c9+c10)*(-1));
            dthetarefdx(17) = m/U_z*((c9+c10)*lambda5-c9*lambda5);

            thetaRef_U_z_coeff = m * ((1-c9^2+lambda5)*e9+(c9+c10)*e10+ref.x_ddot-c9*lambda5*x17);
            dU_zdx = dUzdx(model, state);
            dthetarefdx(3) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(3);
            dthetarefdx(6) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(6);
            dthetarefdx(7) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(7);
            dthetarefdx(8) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(8);
            dthetarefdx(16) = -thetaRef_U_z_coeff/(U_z^2)*dU_zdx(16);
        end

        function dphirefdx = dphirefdx(model, state, ref)
            dphirefdx = zeros(1, model.numStates);
            x2 = state(2);
            x3 = state(3);
            x5 = state(5);
            x6 = state(6);
            x7 = state(7);
            x8 = state(8);
            x16 = state(16);
            x18 = state(18);

            m = model.m;
            c7 = model.c7;
            c8 = model.c8;
            c11 = model.c11;
            c12 = model.c12;
            lambda4 = model.lambda4;
            lambda6 = model.lambda6;

            e9 = ref.y - x2;
            e10 = c11*e9 + ref.y_dot + lambda6*x17-x5;
            
            e7 = ref.z - x3;
            e8 = c7*e7+ref.z_dot + lambda4*x16 - x6;
            U_z = m/(4*cos(x7)*cos(x8))*(g + (1 - c7^2 + lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            
            dphirefdx(2) = m/U_z*(-(1-c11^2+lambda6)*(-1) - (c11+c12)*(-c11));
            dphirefdx(5) = m/U_z*(-(c11+c12)*(-1));
            dphirefdx(18) = m/U_z*(-(c11+c12)*lambda6+c11*lambda6);

            phiRef_U_z_coeff = m * ((1-c11^2+lambda6)*e9+(c11+c12)*e10+ref.y_ddot-c11*lambda6*x18);
            dU_zdx = dUzdx(model, state);
            dphirefdx(3) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(3);
            dphirefdx(6) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(6);
            dphirefdx(7) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(7);
            dphirefdx(8) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(8);
            dphirefdx(16) = -phiRef_U_z_coeff/(U_z^2)*dU_zdx(16);
        end

        function de8dx = de8dx(model)
            de8dx = zeros(1, model.numStates);
            de8dx(4) = -model.c7;
            de8dx(6) = -1;
            de8dx(16) = model.lambda4;
        end

        function de8dotdx = de8dotdx(model, state, ref)
            de8dotdx = df6dx(model, state, ref);
            de8dotdx(3) = -model.lambda4;
            de8dotdx(6) = -model.c7;
        end
        
        function de10dx = de10dx(model)
            de10dx = zeros(1, model.numStates);
            de10dx(1) = -model.c9;
            de10dx(4) = -1;
            de10dx(17) = model.lambda5;
        end

        function de10dotdx = de10dotdx(model, state, ref)
            de10dotdx = df4dx(model, state, ref);
            de10dotdx(1) = -model.lambda5;
            de10dotdx(4) = -model.c9;
        end

        function dezdx = dezdx(model)
            dezdx = (model.c7 + model.c8) * de8dx(model);
            dezdx(3) = dezdx(3)-(1-model.c7^2+model.lambda4);
            dezdx(16) = dezdx(16)-model.c7*model.lambda4;
        end

        function dezdotdx = dezdotdx(model, state, ref)
            dezdotdx = (model.c7 + model.c8) * de8dotdx(model, state, ref);
            dezdotdx(6) = dezdotdx(6)-(1-model.c7^2+model.lambda4);
            dezdotdx(3) = dezdotdx(3)+model.c7*model.lambda4;
        end

        function dUzdotdx = dUzdotdx(model, state, ref)
            dUzdotdx = zeros(1, model.numStates);
            dez_dotdx = dezdotdx(model, state, ref);
            dez_dx = dezdx(model);
            
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
            dUzdotdx(3) = m/(4*denom) * dez_dotdx(3)*(cos(x7)*cos(x8)) - dez_dx(3)*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            dUzdotdx(6) = m/(4*denom) * dez_dotdx(6)*(cos(x7)*cos(x8)) - dez_dx(6)*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            dUzdotdx(16) = m/(4*denom) * dez_dotdx(16)*(cos(x7)*cos(x8)) - dez_dx(16)*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
            
            e7 = ref.z - x3;
            e8 = c7*e7 + ref.z_dot + lambda4*x16 - x6;
            U_z = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);
            z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
            e8_dot = c7*(ref.z_dot - z_dot)*ref.z_ddot+lambda4*e7-z_ddot;
            ez_dot = (1-c7^+lambda4)*(ref.z_dot-z_dot)+(c7+c8)*e8_dot-c7*lambda4*e7;
            ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
            
            Uzdotnum = m/4*(cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)));
            
            dUzdotnumdx7 = m/4*(-sin(x7)*cos(x8)*ez_dot + dez_dotdx(7)*cos(x7)*cos(x8) - ez*(-cos(x7)*cos(x8)+sin(x7)*sin(x8)));
            dUzdotdx(7) = (dUzdotnumdx7*denom-Uzdotnum*(2*cos(x7)*cos(x8)*-sin(x7)))/denom^2;

            dUzdotnumdx8 = m/4*(cos(x7)*-sin(x8)*ez_dot + dez_dotdx(8)*cos(x7)*cos(x8) - ez*(sin(x7)*sin(x8)-cos(x7)*cos(x8)));
            dUzdotdx(8) = (dUzdotnumdx8*denom-Uzdotnum*(2*cos(x7)*cos(x8)*-sin(x8)))/denom^2;
        end

        function df10dx = df10dx(model, state, ref)
            
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
            df13dx = dthetarefdx(model, state, ref);
            df13dx(7) = df13dx(7) - 1;
        end

        function df14dx = df14dx(model, state, ref)
            df14dx = dphirefdx(model, state, ref);
            df14dx(8) = df14dx(8) - 1;
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


classdef derivatives_drone_IB_syms_d2fdx2_11
    methods (Static)
        function d2f11dx2 = d2f11dx2(model, state, ref)
            d2f11dx2 = zeros(model.numStates, model.numStates);
            
            d2f11dx2(1,:)=derivatives_drone_IB_syms_d2fdx2_11_1.d2f11dx2_1(model, state, ref);
            d2f11dx2(2,:)=derivatives_drone_IB_syms_d2fdx2_11_2.d2f11dx2_2(model, state, ref);
            d2f11dx2(3,:)=derivatives_drone_IB_syms_d2fdx2_11_3.d2f11dx2_3(model, state, ref);
            d2f11dx2(4,:)=derivatives_drone_IB_syms_d2fdx2_11_4.d2f11dx2_4(model, state, ref);
            d2f11dx2(5,:)=derivatives_drone_IB_syms_d2fdx2_11_5.d2f11dx2_5(model, state, ref);
            d2f11dx2(6,:)=derivatives_drone_IB_syms_d2fdx2_11_6.d2f11dx2_6(model, state, ref);
            d2f11dx2(7,:)=derivatives_drone_IB_syms_d2fdx2_11_7.d2f11dx2_7(model, state, ref);
            d2f11dx2(8,:)=derivatives_drone_IB_syms_d2fdx2_11_8.d2f11dx2_8(model, state, ref);
            d2f11dx2(9,:)=derivatives_drone_IB_syms_d2fdx2_11_9.d2f11dx2_9(model, state, ref);
            d2f11dx2(10,:)=derivatives_drone_IB_syms_d2fdx2_11_10.d2f11dx2_10(model, state, ref);
            d2f11dx2(11,:)=derivatives_drone_IB_syms_d2fdx2_11_11.d2f11dx2_11(model, state, ref);
            d2f11dx2(12,:)=derivatives_drone_IB_syms_d2fdx2_11_12.d2f11dx2_12(model, state, ref);
            d2f11dx2(13,:)=derivatives_drone_IB_syms_d2fdx2_11_13.d2f11dx2_13(model, state, ref);
            d2f11dx2(14,:)=derivatives_drone_IB_syms_d2fdx2_11_14.d2f11dx2_14(model, state, ref);
            d2f11dx2(16,:)=derivatives_drone_IB_syms_d2fdx2_11_16.d2f11dx2_16(model, state, ref);
            d2f11dx2(17,:)=derivatives_drone_IB_syms_d2fdx2_11_17.d2f11dx2_17(model, state, ref);
            d2f11dx2(18,:)=derivatives_drone_IB_syms_d2fdx2_11_18.d2f11dx2_18(model, state, ref);
            d2f11dx2(19,:)=derivatives_drone_IB_syms_d2fdx2_11_19.d2f11dx2_19(model, state, ref);
            d2f11dx2(20,:)=derivatives_drone_IB_syms_d2fdx2_11_20.d2f11dx2_20(model, state, ref);
            d2f11dx2(21,:)=derivatives_drone_IB_syms_d2fdx2_11_21.d2f11dx2_21(model, state, ref);
            d2f11dx2(22,:)=derivatives_drone_IB_syms_d2fdx2_11_22.d2f11dx2_22(model, state, ref);
            d2f11dx2(23,:)=derivatives_drone_IB_syms_d2fdx2_11_23.d2f11dx2_23(model, state, ref);
            d2f11dx2(24,:)=derivatives_drone_IB_syms_d2fdx2_11_24.d2f11dx2_24(model, state, ref);           
        end
    end
end
classdef derivatives_drone_IB_syms_d2fdx2_10
    methods (Static)
        function d2f10dx2 = d2f10dx2(model, state, ref)
            d2f10dx2 = zeros(model.numStates, model.numStates);
            
            d2f10dx2(1,:)=derivatives_drone_IB_syms_d2fdx2_10_1.d2f10dx2_1(model, state, ref);
            d2f10dx2(2,:)=derivatives_drone_IB_syms_d2fdx2_10_2.d2f10dx2_2(model, state, ref);
            d2f10dx2(3,:)=derivatives_drone_IB_syms_d2fdx2_10_3.d2f10dx2_3(model, state, ref);
            d2f10dx2(4,:)=derivatives_drone_IB_syms_d2fdx2_10_4.d2f10dx2_4(model, state, ref);
            d2f10dx2(5,:)=derivatives_drone_IB_syms_d2fdx2_10_5.d2f10dx2_5(model, state, ref);
            d2f10dx2(6,:)=derivatives_drone_IB_syms_d2fdx2_10_6.d2f10dx2_6(model, state, ref);
            d2f10dx2(7,:)=derivatives_drone_IB_syms_d2fdx2_10_7.d2f10dx2_7(model, state, ref);
            d2f10dx2(8,:)=derivatives_drone_IB_syms_d2fdx2_10_8.d2f10dx2_8(model, state, ref);
            d2f10dx2(9,:)=derivatives_drone_IB_syms_d2fdx2_10_9.d2f10dx2_9(model, state, ref);
            d2f10dx2(10,:)=derivatives_drone_IB_syms_d2fdx2_10_10.d2f10dx2_10(model, state, ref);
            d2f10dx2(11,:)=derivatives_drone_IB_syms_d2fdx2_10_11.d2f10dx2_11(model, state, ref);
            d2f10dx2(12,:)=derivatives_drone_IB_syms_d2fdx2_10_12.d2f10dx2_12(model, state, ref);
            d2f10dx2(13,:)=derivatives_drone_IB_syms_d2fdx2_10_13.d2f10dx2_13(model, state, ref);
            d2f10dx2(14,:)=derivatives_drone_IB_syms_d2fdx2_10_14.d2f10dx2_14(model, state, ref);
            d2f10dx2(16,:)=derivatives_drone_IB_syms_d2fdx2_10_16.d2f10dx2_16(model, state, ref);
            d2f10dx2(17,:)=derivatives_drone_IB_syms_d2fdx2_10_17.d2f10dx2_17(model, state, ref);
            d2f10dx2(18,:)=derivatives_drone_IB_syms_d2fdx2_10_18.d2f10dx2_18(model, state, ref);
            d2f10dx2(19,:)=derivatives_drone_IB_syms_d2fdx2_10_19.d2f10dx2_19(model, state, ref);
            d2f10dx2(20,:)=derivatives_drone_IB_syms_d2fdx2_10_20.d2f10dx2_20(model, state, ref);
            d2f10dx2(21,:)=derivatives_drone_IB_syms_d2fdx2_10_21.d2f10dx2_21(model, state, ref);
            d2f10dx2(22,:)=derivatives_drone_IB_syms_d2fdx2_10_22.d2f10dx2_22(model, state, ref);
            d2f10dx2(23,:)=derivatives_drone_IB_syms_d2fdx2_10_23.d2f10dx2_23(model, state, ref);
            d2f10dx2(24,:)=derivatives_drone_IB_syms_d2fdx2_10_24.d2f10dx2_24(model, state, ref);           
        end
    end
end
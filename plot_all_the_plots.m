function plot_all_the_plots(SC)
    total = int32((SC.stop-SC.start)/SC.step)+1;
    
    % initialize arrays
    x1_val =        zeros([1 total*total]);
    x2_val =        zeros([1 total*total]);

    z =             zeros([1 total*total]);
    ztest =         zeros([1 total*total]);
    z_diff =        zeros([1 total*total]);
    
    zparam =        zeros([1 total*total]);
    zparamtest =    zeros([1 total*total]);
    zparam_diff =   zeros([1 total*total]);
    
    z2 =            zeros([1 total*total]);
    z2test =        zeros([1 total*total]);
    z2_diff =       zeros([1 total*total]);
    
    z2param =       zeros([1 total*total]);
    z2paramtest =   zeros([1 total*total]);
    z2param_diff =  zeros([1 total*total]);

    zdG =           zeros([1 total*total]);
    U =             zeros([1 total*total]);
    L =             zeros([1 total*total]);
    RoA =           zeros([1 total*total]);

    % classes so I can call functions
    tjs = trajectorySensitivities;
    norms = normalStuff;
    splots = splotting;    
    sCtrl = stonkCtrllers;

    % calculate trajectory sensitivities for different initial conditions
    index = 1; 
    for i = SC.start:SC.step:SC.stop
        for j = SC.start:SC.step:SC.stop
            x1_0 = i;
            x2_0 = j;

            initialState.x1 = x1_0;
            initialState.x2 = x2_0;
            initialState.ei = SC.ref-x1_0;
            initialState.e = SC.ref-x1_0;
            initialState.ed = -x2_0;

            % calculate first order derivatives wrt x0
            [stateProgression, ~, ~, ~] = pendulum_trajectory(SC, initialState, SC.controller, @sCtrl.dummyD, false);
            trajSens = tjs.calculateTrajectorySensitivities(stateProgression, SC);
            trajSensTest = tjs.calculateTrajectorySensitivitiesTest(initialState, SC, length(stateProgression));
            G = norms.infimumOfReciprocalOfOneNorm(trajSens);

            % % calculate second derivative wrt x0
            % secondTrajSens = tjs.calculateSecondOrderTrajectorySensitivities(stateProgression, trajSens, SC);
            % secondTrajSensTest = tjs.calculateSecondOrderTrajectorySensitivitiesTest(initialState, SC, length(stateProgression));
            % 
            % % calculate first order trajectory sensitivities wrt p
            % trajSensParam = tjs.calculateTrajectorySensitivitiesParameters(stateProgression, SC);
            % trajSensTestParam = tjs.calculateTrajectorySensitivitiesTestParams(initialState, SC, length(stateProgression));
            % 
            % % calculate second derivative wrt p
            % secondTrajSensParam = tjs.calculateSecondOrderTrajectorySensitivitiesParams(stateProgression, trajSens, trajSensParam, SC);
            % secondTrajSensTestParam = tjs.calculateSecondOrderTrajectorySensitivitiesTestParams(initialState, SC, length(stateProgression));
            % 
            % tp = tjs.findtp(trajSens);
            % 
            % dGp = tjs.calculate_dG(G, trajSens, secondTrajSensParam, tp, SC);

            z(index) = G;
            ztest(index) = norms.infimumOfReciprocalOfOneNorm(trajSensTest);
            % z_diff(index) = norms.supremum(trajSens-trajSensTest);

            % z2(index) = norms.secondOrderMax(secondTrajSens);
            % z2test(index) = norms.secondOrderMax(secondTrajSensTest);
            % z2_diff(index) = norms.secondOrderMax(secondTrajSens-secondTrajSensTest);
            % 
            % zparam(index) = norms.supremum(trajSensParam);
            % zparamtest(index) = norms.supremum(trajSensTestParam);
            % zparam_diff(index) = norms.supremum(trajSensParam-trajSensTestParam);
            % 
            % z2param(index) = norms.secondOrderMax(secondTrajSensParam);
            % z2paramtest(index) = norms.secondOrderMax(secondTrajSensTestParam);
            % z2param_diff(index) = norms.secondOrderMax(secondTrajSensParam-secondTrajSensTestParam);

            % zdG(index) = norm(dGp);
            % U(index) = dGp(1);
            % L(index) = dGp(2);
            
            x1_val(index) = i;
            x2_val(index) = j;

            if (abs(stateProgression(size(stateProgression, 2)).x1) < pi) 
                RoA(index) = 1;
            else 
                RoA(index) = -1;
            end 
            index = index+1;
        end 
    end 
    
    index = index-1;
    splots.splot(1, x1_val(:, 1:index), x2_val(:, 1:index), z(:, 1:index), "infimum inverse trajectory sensitivity");
    splots.splot(2, x1_val(:, 1:index), x2_val(:, 1:index), ztest(:, 1:index), "infimum inverse trajectory sensitivity test");
    % splots.splot(3, x1_val(:, 1:index), x2_val(:, 1:index), z_diff(:, 1:index), "max diff");
    % disp("max z diff: " + max(z_diff));

    % splots.splot(4, x1_val(:, 1:index), x2_val(:, 1:index), z2(:, 1:index), "second order trajectory sensitivity");
    % splots.splot(5, x1_val(:, 1:index), x2_val(:, 1:index), z2test(:, 1:index), "second order trajectory sensitivity test");
    % splots.splot(6, x1_val(:, 1:index), x2_val(:, 1:index), z2_diff(:, 1:index), "max diff second order trajectory sensitivity");
    % % disp("max z2 diff: " + max(z2_diff));
    % 
    % splots.splot(7, x1_val(:, 1:index), x2_val(:, 1:index), zparam(:, 1:index), "supremum trajectory sensitivity param");
    % splots.splot(8, x1_val(:, 1:index), x2_val(:, 1:index), zparamtest(:, 1:index), "supremum trajectory sensitivity test param");
    % splots.splot(9, x1_val(:, 1:index), x2_val(:, 1:index), zparam_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % 
    % splots.splot(10, x1_val(:, 1:index), x2_val(:, 1:index), z2param(:, 1:index), "second order trajectory sensitivity param");
    % splots.splot(11, x1_val(:, 1:index), x2_val(:, 1:index), z2paramtest(:, 1:index), "second order trajectory sensitivity test param");
    % splots.splot(12, x1_val(:, 1:index), x2_val(:, 1:index), z2param_diff(:, 1:index), "supremum diff trajectory sensitivity param");
    % zlim([0, 30]);

    % splots.splot(13, x1_val(:, 1:index), x2_val(:, 1:index), zdG(:, 1:index), "dG");
    % splots.squiver(14, x1_val(:, 1:index), x2_val(:, 1:index), U(:, 1:index), L(:, 1:index), "u l dG");

    splots.splot(15, x1_val(:, 1:index), x2_val(:, 1:index), RoA(:, 1:index), "RoA");
end
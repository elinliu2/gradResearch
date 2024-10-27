function plot_all_the_plots(SC)
    total = int32((SC.stop-SC.start)/SC.step)+1;
    
    % initialize arrays
    x1_val =        zeros([1 total*total]);
    x2_val =        zeros([1 total*total]);
    z =             zeros([1 total*total]);
    ztest =         zeros([1 total*total]);
    z_param =       zeros([1 total*total]);
    z_paramtest =   zeros([1 total*total]);
    z_paramdiff =   zeros([1 total*total]);
    z_diff =        zeros([1 total*total]);
    z_sup =         zeros([1 total*total]);
    z_suptest =     zeros([1 total*total]);
    z2 =            zeros([1 total*total]);
    z2test =        zeros([1 total*total]);
    z2param =       zeros([1 total*total]);
    z2testparam =   zeros([1 total*total]);
    z2paramdiff =   zeros([1 total*total]);
    z2_diff =       zeros([1 total*total]);
    zdG =           zeros([1 total*total]);
    zdGtest =       zeros([1 total*total]);
    zdGdiff =       zeros([1 total*total]);
    U =             zeros([1 total*total]);
    V =             zeros([1 total*total]);
    RoA =           zeros([1 total*total]);

    % classes so I can call functions
    deriv = derivatives;
    norms = normalStuff;
    splots = splotting;    
    sCtrl = stonkCtrllers;

    % calculate trajectory sensitivities for different initial conditions
    index = 1; 
    for i = SC.start:SC.step:SC.stop
        for j = SC.start:SC.step:SC.stop
            x1_0 = i;
            x2_0 = j;
            % SC.u = i; 
            % SC.l = j;
            % calculate first order derivatives wrt x0
            [x1, x2] = pendulum_trajectory(SC, x1_0, x2_0, @sCtrl.dummyU, @sCtrl.dummyD, false);
            % trajSens = deriv.calculateTrajectorySensitivities(x1, x2, SC);
    
            % validate results for first derivative
            % trajSensTest = deriv.calculateTrajectorySensitivitiesTest(i, j, SC);
    
            % calculate second derivative
            % secondTrajSens = deriv.calculateSecondOrderTrajectorySensitivities(x1, x2, trajSens, SC);
    
            % validate results for second derivative
            % secondTrajSensTest = deriv.calculateSecondOrderTrajectorySensitivitiesTest(i, j, SC);
    
            % calculate dG
            % tp = deriv.findtp(trajSens);
            % dG = deriv.calculate_dG(trajSens, secondTrajSens, tp);
            % dGtest = deriv.calculate_dGtest(i, j, initdxdxo, SC);
            
            % G = norms.infimumOfReciprocalOfOneNorm(trajSens);
            % z(index) = G
            % ztest(index) = norms.infimumOfReciprocalOfOneNorm(trajSensTest);
            % z_diff(index) = norms.infimum(1./trajSens);
            % 
            % z_sup(index) = norms.supremum(trajSens);
            % z_suptest(index) = norms.supremum(trajSensTest);
    
            % z2(index) = norms.secondOrderMax(secondTrajSens);
            % z2test(index) = norms.secondOrderMax(secondTrajSensTest);
            % z2_diff(index) = norms.secondOrderMax(secondTrajSens-secondTrajSensTest);
            
            % zdG(index) = norm(dG);
            % zdGtest(index) = norm(dGtest);
            % zdGdiff(index) = norm(dG - dGtest);

            % calculate first order trajectory sensitivities wrt p
           
            % trajSensParam = deriv.calculateTrajectorySensitivitiesParameters(x1, x2, SC);

            % validate results wrt p
            % trajSensTestParam = deriv.calculateTrajectorySensitivitiesTestParams(i, j, 0, 0, SC);
            
            % z_param(index) = norms.infimumOfReciprocalOfOneNorm(trajSensParam);
            % z_paramtest(index) = norms.infimumOfReciprocalOfOneNorm(trajSensTestParam);
            % z_paramdiff(index) = norms.supremum(trajSensParam-trajSensTestParam);

            % secondTrajSensParam = deriv.calculateSecondOrderTrajectorySensitivitiesParams(x1, x2, trajSens, trajSensParam, SC);
            % secondTrajSensTestParam = deriv.calculateSecondOrderTrajectorySensitivitiesTestParams(i, j, SC);
            
            % z2param(index) = norms.secondOrderMax(secondTrajSensParam);
            % z2testparam(index) = norms.secondOrderMax(secondTrajSensTestParam);
            % z2paramdiff(index) = norms.secondOrderMax(secondTrajSensParam-secondTrajSensTestParam);

            % tp = deriv.findtp(trajSens);
            % dGp = deriv.calculate_dG(G, trajSensParam, secondTrajSensParam, tp);
            % 
            % zdG(index) = norm(dGp);
            % 
            % U(index) = dGp(1);
            % V(index) = dGp(2);

            if abs(x1(size(x1, 2))) < pi
                RoA(index) = 1;
            else
                RoA(index) = -1;
            end
       
            x1_val(index) = i;
            x2_val(index) = j;
            index = index+1;
        end 
    end 
    
    index = index-1;
    % splots.splot(1, x1_val(:, 1:index), x2_val(:, 1:index), z(:, 1:index), "infimum inverse trajectory sensitivity");
    % splots.splot(2, x1_val(:, 1:index), x2_val(:, 1:index), ztest(:, 1:index), "infimum inverse trajectory sensitivity test");
    % splots.splot(3, x1_val(:, 1:index), x2_val(:, 1:index), z_diff(:, 1:index), "wrong infimum");
    % % disp("max z diff: " + max(z_diff));
    % 
    % splots.splot(4, x1_val(:, 1:index), x2_val(:, 1:index), 1./z_sup(:, 1:index), "inverse supremum trajectory sensitivity");
    % splots.splot(4, x1_val, x2_val, 1./z_suptest, "inverse supremum trajectory sensitivity test");
    
    % splots.splot(5, x1_val, x2_val, z2, "second order trajectory sensitivity");
    % splots.splot(6, x1_val, x2_val, z2test, "second order trajectory sensitivity test");
    % splots.splot(7, x1_val, x2_val, z2_diff, "max diff second order trajectory sensitivity");
    % disp("max z2 diff: " + max(z2_diff));
    
    % splots.splot(8, x1_val, x2_val, zdG, "dG");
    % splots.splot(9, x1_val, x2_val, zdGtest, "dG test");
    % 
    % splots.splot(10, x1_val, x2_val, zdG, "dG");
    % zlim([0, 100]);
    % splots.splot(11, x1_val, x2_val, zdGtest, "dG test");
    % zlim([0, 100]);
    % 
    % splots.splot(12, x1_val, x2_val, zdGdiff, "dG diff");
    % 
    % figure(13)
    % hold off
    % histogram(zdGdiff/zdG);

    % splots.splot(14, x1_val(:, 1:index), x2_val(:, 1:index), z_param(:, 1:index), "infimum inverse trajectory sensitivity param");
    % splots.splot(15, x1_val(:, 1:index), x2_val(:, 1:index), z_paramtest(:, 1:index), "infimum inverse trajectory sensitivity test param");
    % splots.splot(16, x1_val(:, 1:index), x2_val(:, 1:index), z_paramdiff(:, 1:index), "supremum diff trajectory sensitivity param");

    % splots.splot(17, x1_val(:, 1:index), x2_val(:, 1:index), z2param(:, 1:index), "second order trajectory sensitivity param");
    % zlim([0, 600]);
    % splots.splot(18, x1_val(:, 1:index), x2_val(:, 1:index), z2testparam(:, 1:index), "second order trajectory sensitivity test param");
    % zlim([0, 600]);
    % splots.splot(19, x1_val(:, 1:index), x2_val(:, 1:index), z2paramdiff(:, 1:index), "supremum diff trajectory sensitivity param");
    % zlim([0, 30]);

    splots.splot(20, x1_val(:, 1:index), x2_val(:, 1:index), RoA(:, 1:index), "RoA");

    % splots.squiver(21, x1_val, x2_val, U, V, "u l trajectory");
end
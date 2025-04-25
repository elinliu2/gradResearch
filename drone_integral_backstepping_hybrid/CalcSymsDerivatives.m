clc;
syms c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 
syms m g Jtheta Jphi Jpsi K_psi l
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23
syms refx refxdot refxddot refxdddot refxddddot
syms refy refydot refyddot refydddot refyddddot
syms refz refzdot refzddot refzdddot refzddddot
syms refpsi refpsidot refpsiddot refpsidddot refpsiddddot

e7 = refz - x3;
e8 = c7*e7 + refzdot + lambda4*x16 - x6;
U_z = x19;

e9 = refx - x1;
e10 = c9*e9 + refxdot + lambda5*x17 - x4;
theta_ref = m/U_z*((1-c9^2  + lambda5)*e9 + (c9 +c10)*e10 + refxddot - c9*lambda5*x17);
e11 = refy - x2;
e12 = c11*e11 + refydot + lambda6*x18 - x5;
phi_ref =   m/U_z*( -(1-c11^2 + lambda6)*e11 - (c11+c12)*e12 - refyddot + c11*lambda6*x18);

theta_ref_dot = x20;
phi_ref_dot = x21;
theta_ref_ddot = x22;
phi_ref_ddot = x23;

e1 = theta_ref - x7;
e2 = c1*e1 + theta_ref_dot + lambda1*x13 - x10;
e3 = phi_ref - x8;
e4 = c3*e3 + phi_ref_dot + lambda2*x14 - x11;
e5 = refpsi - x9;
e6 = c5*e5 + refpsidot + lambda3*x15 - x12;

z_ddot = 4*U_z*cos(x7)*cos(x8)/m-g;
e8_dot = c7*(refzdot - x6)+refzddot+lambda4*e7-z_ddot;
ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
ez_dot = ((1-c7^2+lambda4)*(refzdot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7);
U_zdot = m/4*( cos(x7)*cos(x8)*ez_dot + ez*(-sin(x7)*cos(x8)*x10-cos(x7)*sin(x8)*x11) )/(cos(x7)*cos(x8))^2;

U_theta = Jtheta/(2*l)*((1-c1^2+lambda1)*e1 + (c1+c2)*e2 + theta_ref_ddot - c1*lambda1*x13);
U_phi = Jphi/(2*l)*((1-c3^2+lambda2)*e3 + (c3+c4)*e4 + phi_ref_ddot - c3*lambda2*x14);
U_psi = Jpsi/(4*K_psi)*((1-c5^2+lambda3)*e5 + (c5+c6)*e6 + refpsiddot - c5*lambda3*x15);

f1 = x4; 
f2 = x5;
f3 = x6;

f4 = 4*U_z/m*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9));
f5 = 4*U_z/m*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9));
f6 = 4*U_z/m*(cos(x7)*cos(x8)) - g;

f7 = x10; 
f8 = x11;
f9 = x12;

f10 = 1/(Jtheta) * l * (2*U_theta); 
f11 = 1/(Jphi) * l * (2*U_phi);
f12 = 1/(Jpsi) * K_psi * (4*U_psi);

f13 = theta_ref - x7;
f14 = phi_ref - x8;
f15 = refpsi - x9;
f16 = refz - x3;
f17 = refx - x1;
f18 = refy - x2;

f19 = U_zdot;

  func = [f1; f2; f3; f4; f5; f6; f7; f8; f9; f10; f11; f12; f13; f14; f15; f16; f17; f18; f19; theta_ref; phi_ref];
  vars = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18 x19 x20 x21 x22 x23];
  params = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6];
  
filename = "C:\Users\elinl\OneDrive - University of Waterloo\research\derivatives.txt";

% for i = 1:23
%     derivative = diff(f19, vars(i));
%     if(derivative ~= 0)
%         disp(strcat("           dfdx(", string(i) ,")=", string(derivative),";"))
%         % lines = strcat("           df24dx(", string(i) ,")=", string(derivative),";");
%         % writelines(lines,filename,WriteMode="append")
%     end
% end
 
for i = 1:18
    derivative = diff(phi_ref, params(i));
    if(derivative ~= 0)
        disp(strcat("           dgdq(2,", string(23+i) ,")=", string(derivative),";"))
        % lines = strcat("           df1dp(", string(i) ,")=", string(derivative),";");
        % writelines(lines,filename,WriteMode="append")
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% calc d2fdx2

% df4dx_7=(4*x19*cos(x7)*cos(x8)*cos(x9))/m;
% df4dx_8=(4*x19*(cos(x8)*sin(x9) - cos(x9)*sin(x7)*sin(x8)))/m;
% df4dx_9=(4*x19*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
% df4dx_19=(4*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
% dfunc = [ df4dx_7 df4dx_8 df4dx_9 df4dx_19];
% dfunc_indicies = [7;8;9;19];

% df5dx_7=(4*x19*cos(x7)*cos(x8)*sin(x9))/m;
% df5dx_8=-(4*x19*(cos(x8)*cos(x9) + sin(x7)*sin(x8)*sin(x9)))/m;
% df5dx_9=(4*x19*(sin(x8)*sin(x9) + cos(x8)*cos(x9)*sin(x7)))/m;
% df5dx_19=-(4*(cos(x9)*sin(x8) - cos(x8)*sin(x7)*sin(x9)))/m;
% dfunc = [df5dx_7 df5dx_8 df5dx_9 df5dx_19  ];
% dfunc_indicies = [7;8;9;19];

% df6dx_7=-(4*x19*cos(x8)*sin(x7))/m;
% df6dx_8=-(4*x19*cos(x7)*sin(x8))/m;
% df6dx_19=(4*cos(x7)*cos(x8))/m;
% dfunc = [df6dx_7 df6dx_8 df6dx_19 ];
% dfunc_indicies = [7;8;19];

% df10dx_1=- (m*(lambda1 - c1^2 + 1)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19 - (c1*m*(c1 + c2)*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
% df10dx_4=- (m*(c9 + c10)*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2)*(c9 + c10))/x19;
% df10dx_7=c1^2 - lambda1 - c1*(c1 + c2) - 1;
% df10dx_10=- c1 - c2;
% df10dx_13=lambda1*(c1 + c2) - c1*lambda1;
% df10dx_17=- (m*(c9*lambda5 - lambda5*(c9 + c10))*(lambda1 - c1^2 + 1))/x19 - (c1*m*(c1 + c2)*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
% df10dx_19=- (m*(lambda1 - c1^2 + 1)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2 - (c1*m*(c1 + c2)*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
% 
% dfunc = [df10dx_1 df10dx_4 df10dx_7 df10dx_10 df10dx_13 df10dx_17 df10dx_19 ];
% dfunc_indicies = [1;4;7;10;13;17;19];

% df11dx_2=(m*(lambda2 - c3^2 + 1)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19 + (c3*m*(c3 + c4)*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
% df11dx_5=(m*(c11 + c12)*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4)*(c11 + c12))/x19;
% df11dx_8=c3^2 - lambda2 - c3*(c3 + c4) - 1;
% df11dx_11=- c3 - c4;
% df11dx_14=lambda2*(c3 + c4) - c3*lambda2;
% df11dx_18=(m*(c11*lambda6 - lambda6*(c11 + c12))*(lambda2 - c3^2 + 1))/x19 + (c3*m*(c3 + c4)*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
% df11dx_19=(m*(lambda2 - c3^2 + 1)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2 + (c3*m*(c3 + c4)*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
% 
% dfunc = [df11dx_2 df11dx_5 df11dx_8 df11dx_11 df11dx_14 df11dx_18 df11dx_19 ];
% dfunc_indicies = [2;5;8;11;14;18;19];

% df12dx_9=c5^2 - lambda3 - c5*(c5 + c6) - 1;
% df12dx_12=- c5 - c6;
% df12dx_15=lambda3*(c5 + c6) - c5*lambda3;
% 
% dfunc = [df12dx_9 df12dx_12 df12dx_15 ];
% dfunc_indicies = [9;12;15];

% df13dx_1=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
% df13dx_4=-(m*(c9 + c10))/x19;
% df13dx_7=-1;
% df13dx_17=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
% df13dx_19=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
% 
% dfunc = [df13dx_1 df13dx_4 df13dx_7 df13dx_17 df13dx_19];
% dfunc_indicies = [1;4;7;17;19];

% df14dx_2=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
% df14dx_5=(m*(c11 + c12))/x19;
% df14dx_8=-1;
% df14dx_18=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
% df14dx_19=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
% 
% dfunc = [df14dx_2 df14dx_5 df14dx_8 df14dx_18 df14dx_19];
% dfunc_indicies = [2;5;8;18;19];

% df19dx_3=(m*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(lambda4 - c7^2 + c7*(c7 + c8) + 1) + cos(x7)*cos(x8)*(c7*lambda4 - lambda4*(c7 + c8))))/(4*cos(x7)^2*cos(x8)^2);
% df19dx_6=(m*((c7 + c8)*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)) - cos(x7)*cos(x8)*(lambda4 - c7^2 + c7*(c7 + c8) + 1)))/(4*cos(x7)^2*cos(x8)^2);
% df19dx_7=- (m*((x10*cos(x7)*cos(x8) - x11*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x8)*sin(x7)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)*cos(x8)^2*sin(x7)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x7)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^3*cos(x8)^2);
% df19dx_8=- (m*((x11*cos(x7)*cos(x8) - x10*sin(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) + cos(x7)*sin(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3)) - (4*x19*cos(x7)^2*cos(x8)*sin(x8)*(c7 + c8))/m))/(4*cos(x7)^2*cos(x8)^2) - (m*sin(x8)*((x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8))*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16) - cos(x7)*cos(x8)*((refzdot - x6)*(lambda4 - c7^2 + 1) + (c7 + c8)*(g + refzddot + c7*(refzdot - x6) + lambda4*(refz - x3) - (4*x19*cos(x7)*cos(x8))/m) - c7*lambda4*(refz - x3))))/(2*cos(x7)^2*cos(x8)^3);
% df19dx_10=-(m*sin(x7)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)^2*cos(x8));
% df19dx_11=-(m*sin(x8)*(g + (c7 + c8)*(refzdot - x6 + lambda4*x16 + c7*(refz - x3)) + (refz - x3)*(lambda4 - c7^2 + 1) - c7*lambda4*x16))/(4*cos(x7)*cos(x8)^2);
% df19dx_16=(m*(c7*lambda4 - lambda4*(c7 + c8))*(x10*cos(x8)*sin(x7) + x11*cos(x7)*sin(x8)))/(4*cos(x7)^2*cos(x8)^2);
% df19dx_19=- c7 - c8;
% 
% dfunc = [df19dx_3 df19dx_6 df19dx_7 df19dx_8 df19dx_10 df19dx_11 df19dx_16 df19dx_19];
% dfunc_indicies = [3;6;7;8;10;11;16;19];

% dthetarefdx_1=-(m*(lambda5 - c9^2 + c9*(c9 + c10) + 1))/x19;
% dthetarefdx_4=-(m*(c9 + c10))/x19;
% dthetarefdx_17=-(m*(c9*lambda5 - lambda5*(c9 + c10)))/x19;
% dthetarefdx_19=-(m*(refxddot + (c9 + c10)*(refxdot - x4 + lambda5*x17 + c9*(refx - x1)) + (refx - x1)*(lambda5 - c9^2 + 1) - c9*lambda5*x17))/x19^2;
% 
% dfunc = [dthetarefdx_1 dthetarefdx_4 dthetarefdx_17 dthetarefdx_19 ];
% dfunc_indicies = [1;4;17;19];

% dphirefdx_2=(m*(lambda6 - c11^2 + c11*(c11 + c12) + 1))/x19;
% dphirefdx_5=(m*(c11 + c12))/x19;
% dphirefdx_18=(m*(c11*lambda6 - lambda6*(c11 + c12)))/x19;
% dphirefdx_19=(m*(refyddot + (c11 + c12)*(refydot - x5 + lambda6*x18 + c11*(refy - x2)) + (refy - x2)*(lambda6 - c11^2 + 1) - c11*lambda6*x18))/x19^2;
% 
% dfunc = [dphirefdx_2 dphirefdx_5 dphirefdx_18 dphirefdx_19 ];
% dfunc_indicies = [2; 5; 18; 19];


% for i = 1:length(dfunc_indicies)
%     index = dfunc_indicies(i);
%     for j = 1:19
%         derivative = diff(dfunc(i), vars(j));
%         if(derivative ~= 0)
%             disp(strcat("           d2f11dx2(", num2str(index), ",", num2str(j), ")=", string(derivative),";"));
%         end
%     end
% end


% for i = 1:length(dfunc_indicies)
%     index = dfunc_indicies(i);
%     for j = 1:18
%         derivative = diff(dfunc(i), params(j));
%         if(derivative ~= 0)
%             disp(strcat("           d2thetarefdxdp(", num2str(index), ",", num2str(j), ")=", string(derivative),";"))
%         end
%     end
% end





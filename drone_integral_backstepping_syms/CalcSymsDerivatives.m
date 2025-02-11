syms m c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 g Jtheta Jphi Jpsi K_psi l
syms x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18
syms refx refxdot refxddot refxdddot refxddddot
syms refy refydot refyddot refydddot refyddddot
syms refz refzdot refzddot refzdddot refzddddot
syms refpsi refpsidot refpsiddot refpsidddot refpsiddddot

e7 = refz - x3;
e8 = c7*e7 + refzdot + lambda4*x16 - x6;
Uz = m/(cos(x7)*cos(x8)*4)*(g + (1-c7^2+lambda4)*e7 + (c7+c8)*e8 - c7*lambda4*x16);

e9 = refx - x1;
e10 = c9*e9 + refxdot + lambda5*x17 - x4;
theta_ref = m/Uz*((1-c9^2  + lambda5)*e9 + (c9 +c10)*e10 + refxddot - c9*lambda5*x17);

e11 = refy - x2;
e12 = c11*e11 + refydot + lambda6*x18 - x5;
phi_ref =   m/Uz*( -(1-c11^2 + lambda6)*e11 - (c11+c12)*e12 - refyddot + c11*lambda6*x18);

z_ddot = 4*Uz*cos(x7)*cos(x8)/m-g;
e8_dot = c7*(refzdot - x6)*refzddot+lambda4*e7-z_ddot;
ez_dot = ((1-c7^2-lambda4)*(refzdot-x6)+(c7+c8)*e8_dot-c7*lambda4*e7);
ez     = (1-c7^2+lambda4)*e7+(c7+c8)*e8-c7*lambda4*x16+g;
U_zdot = m/4*(cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8)))/(cos(x7)*cos(x8))^2;

x_ddot = 4*Uz*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9))/(m);
e10_dot= c9*(refxdot - x4) - refxddot + lambda5*e9 - x_ddot;
ex_dot = (1-c9^2+lambda5)*(refxdot-x4)+(c9+c10)*e10_dot+refxdddot-c9*lambda5*e9;
ex     = (1-c9^2+lambda5)*e9 + (c9+c10)*e10 + refxddot - c9*lambda5*x17;

y_ddot = 4*Uz*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9));
e12_dot= c11*(refydot - x5) - refyddot + lambda6*e11 - y_ddot;
ey_dot = -(1-c11^2+lambda6)*(refydot-x5)-(c11+c12)*e12_dot-refydddot+c11*lambda6*e11;
ey     = -(1-c11^2+lambda6)*e11 - (c11+c12)*e12 - refyddot + c11*lambda6*x18;

theta_ref_dot = m * (Uz*ex_dot - ex*U_zdot)/Uz^2;
phi_ref_dot = m * (Uz*ey_dot - ey*U_zdot)/Uz^2;

x_dddot = 4/m*( ...
    U_zdot*(sin(x8)*sin(x9) + sin(x7)*cos(x8)*cos(x9)) + ...
    Uz*(cos(x8)*sin(x9)*x11 + sin(x8)*cos(x9)*x12 + ...
         cos(x7)*x10*cos(x8)*cos(x9) - ...
         sin(x7)*sin(x8)*x11*cos(x9) - ...
         sin(x7)*cos(x8)*sin(x9)*x12));

e10_ddot = c9*(refxddot - x_ddot) - refxdddot + lambda5*(refxdot - x4) - x_dddot;
ex_ddot = -(1-c9^2+lambda5)*(refxddot - x_ddot) - (c9+c10)*e10_ddot - refxddddot + c9*lambda5*(refxdot - x4);

z_dddot = 4/(m)*(U_zdot*cos(x8)*cos(x7)-Uz*sin(x8)*x11*cos(x7)-Uz*cos(x7)*sin(x7)*x10);
e8_ddot = c7*(refzddot-z_ddot)+refzdddot+lambda4*(refzdot-x6)-z_dddot;
ez_ddot = (1-c7^2-lambda4)*(refzddot-z_ddot)+(c7+c8)*e8_ddot-c7*lambda4*(refzdot-x6);
U_zdot_top_dot = -sin(x7)*cos(x8)*ez_dot - cos(x7)*sin(x8)*ez_dot + cos(x7)*cos(x8)*ez_ddot + ...
    cos(x7)*cos(x8)*ez - sin(x7)*sin(x8)*ez + sin(x7)*cos(x8)*ez_dot - ...
    sin(x7)*sin(x8)*ez + cos(x7)*cos(x8)*ez + cos(x7)*sin(x8)*ez_dot;
U_zdot_top = cos(x7)*cos(x8)*ez_dot - ez*(-sin(x7)*cos(x8)-cos(x7)*sin(x8));
U_zddot = m/4*(cos(x7)*cos(x8)^2*U_zdot_top_dot - U_zdot_top*2*cos(x7)*cos(x8)*(-sin(x7*cos(x8)-cos(x7)*sin(x8))))/(cos(x7)*cos(x8))^4;

komh_dot = (Uz*ex_ddot)-(ex*U_zddot);
komh = Uz*ex_dot - ex*U_zdot;
theta_ref_ddot = m*(Uz^2*komh_dot-komh*2*Uz*U_zdot)/Uz^4;

y_dddot = 4/m*(U_zdot*(-sin(x8)*cos(x9)+sin(x7)*cos(x8)*sin(x9)) ...
        + Uz *(-cos(x8)*x11*cos(x9)+sin(x8)*sin(x9)*x12 ...
                +cos(x7)*x10*cos(x9)*sin(x9) ...
                -sin(x7)*sin(x8)*x11*sin(x9) ...
                +sin(x7)*cos(x8)*cos(x9)*x12));
e12_ddot = c11*(refyddot-y_ddot) - refydddot + lambda6*(refydot-x5)-y_dddot;
ey_ddot = (1-c11^2+lambda6)*(refyddot-y_ddot)+(c11+c12)*e12_ddot+refyddddot-c11*lambda6*(refydot-x5);
ciwyw_dot = (U_zdot*ey_dot+Uz*ey_ddot)-(U_zdot*ey_dot+ey*U_zddot);
ciwyw = Uz*ey_dot - ey*U_zdot;
phi_ref_ddot = m*(Uz^2*ciwyw_dot-ciwyw*2*Uz*U_zdot)/Uz^4;

e1 = theta_ref - x7;
e2 = c1*e1 + theta_ref_dot + lambda1*x13 - x10;
U_theta = Jtheta/(2*l)*((1-c1^2+lambda1)*e1 + (c1+c2)*e2 + theta_ref_ddot - c1*lambda1*x13);

e3 = phi_ref - x8;
e4 = c3*e3 + phi_ref_dot + lambda2*x14 - x11;
U_phi = Jphi/(2*l)*((1-c3^2+lambda2)*e3 + (c3+c4)*e4 + phi_ref_ddot - c3*lambda2*x14);

e5 = refpsi - x9;
e6 = c5*e5 + refpsidot + lambda3*x15 - x12;
U_psi = Jpsi/(4*K_psi)*((1-c5^2+lambda3)*e5 + (c5+c5)*e6 + refpsiddot - c5*lambda3*x15);


f4 = 4*Uz/m*(sin(x7)*cos(x8)*cos(x9) + sin(x8)*sin(x9));
f5 = 4*Uz/m*(sin(x7)*cos(x8)*sin(x9) - sin(x8)*cos(x9));
f6 = 4*Uz/m*((cos(x7)*cos(x8)) - g);

f10 = 1/(Jtheta) * l * (2*U_theta); 
f11 = 1/(Jphi) * l * (2*U_phi);
f12 = 1/(Jpsi) * K_psi * (4*U_psi);

f13 = theta_ref - x7;
f14 = phi_ref - x8;
f15 = refpsi - x9;
f16 = refz - x3;
f17 = refx - x1;
f18 = refy - x2;

% vars = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
% 
% for i = 1:18
%     disp(i);
%     disp(diff(f15, vars(i)));
% end

vars = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14 x15 x16 x17 x18];
coeff = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 lambda1 lambda2 lambda3 lambda4 lambda5 lambda6 m g l Jtheta Jphi Jpsi K_psi];
for i = 1:25
    disp(i);
    disp(diff(f14, coeff(i)));
end




n = 23;
p = 18;
z0    = sym('z0', [n 1]);
zk    = sym('zk', [n 1]);
pk    = sym('pk', [p 1]);
pkp1  = sym('pkp1', [p 1]);
DG    = sym('DG', [1 n+p]);
alpha = sym('alpha', [1 1]);

vz = DG*[eye(n); zeros(p, n)];
vp = DG*[zeros(n, p); eye(p)];

% zkp1_1 = [eye(n) zeros(n, 1)]*inv([eye(n) transpose(vz); vz 0])*[z0; DG*[zk;pk]-vp*pkp1];
% zkp1_2 = z0 + transpose(vz)/(vz*transpose(vz))*(-vz*z0+DG*[zk;pk]-vp*pkp1);
% simplify(zkp1_1-zkp1_2)

% dfdp1 = transpose(transpose(z0 + transpose(vz)/(vz*transpose(vz))*(-vz*z0+DG*[zk;pk]-vp*pk) - z0)*(-transpose(vz)*vp/(vz*transpose(vz))));
% dfdp2 = transpose(vp)*vz/(vz*transpose(vz))*(z0-zk);
% simplify(dfdp1 - dfdp2)

% diff1 = z0 - (z0 + transpose(vz)/(vz*transpose(vz))*(-vz*z0+DG*[zk;pk]-vp*(pk + alpha*transpose(vp)*vz/(vz*transpose(vz))*(z0-zk))));
% diff2 = transpose(vz)/(vz*transpose(vz))*(vz*z0 - vz*zk - vp*pk + vp*(pk + alpha*transpose(vp)*vz/(vz*transpose(vz))*(z0-zk)));
% diff3 = transpose(vz)/(vz*transpose(vz))*(vz*z0 - vz*zk  + alpha*vp*transpose(vp)*vz/(vz*transpose(vz))*(z0-zk));
% diff4 = 1/(vz*transpose(vz))^2*(vz*transpose(vz)*transpose(vz)*vz + alpha*transpose(vz)*vp*transpose(vp)*vz)*(z0-zk);
% simplify(diff3 - diff4)

simplify(norm( (vz*transpose(vz)*transpose(vz)*vz)/(vz*transpose(vz))^2 ))
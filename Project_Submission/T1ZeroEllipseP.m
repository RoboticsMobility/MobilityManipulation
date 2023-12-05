function qd = T1ZeroEllipseP(t,q)

% q = [x , y , phi, theR, theL, the1, the2]
L1 = 1.5; L2 = 1;
r = 0.5;
bw = 0.6;

a = 2.5;
b = 0.5;

dXed = [-a*cosd(60)*sin((2*pi/20)*t)-b*sind(60)*cos((2*pi/20)*t);
    -a*sind(60)*sin((2*pi/20)*t)+b*cosd(60)*cos((2*pi/20)*t)]*(2*pi/20);

Xe = [L1*cos(q(6))+L2*cos(q(6)+q(7));
      L1*sin(q(6))+L2*sin(q(6)+q(7))];

J1 = (r/2)*[1 1;0 0];
J2 = (r/2*bw)*[Xe(2) -Xe(2);-Xe(1) Xe(1)];
J3 = [-L1*sin(q(6))-L2*sin(q(6)+q(7)) -L2*sin(q(6)+q(7));L1*cos(q(6))+L2*cos(q(6)+q(7)) L2*cos(q(6)+q(7))];

R = [cos(q(3)) -sin(q(3)); sin(q(3)) cos(q(3))];

J = R*[J1+J2 J3];

qtemp = pinv(J)*dXed;

qnull = pinv([0 0 1 0])*0;

[qd(4:7,1)] = qtemp + (eye(4) - J\J)*qnull;

Jb = [(r/2)*cos(q(3)) (r/2)*cos(q(3));
      (r/2)*sin(q(3)) (r/2)*sin(q(3));
      r/(2*bw) -r/(2*bw)];

[qd(1:3,1)] = Jb*qd(4:5,1);
end
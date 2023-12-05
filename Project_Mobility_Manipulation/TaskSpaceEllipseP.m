function qd = TaskSpaceEllipseP(t,q)

% q = [x , y , phi, theR, theL, the1, the2]
L1 = 1.5; L2 = 1;
r = 0.5;
bw = 0.6;

ec = [0.5 0.5];
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

R = [cosd(60) -sind(60);
    sind(60) cosd(60)];
    
xe = a*cos((2*pi/20)*t);
ye = b*sin((2*pi/20)*t);

Xedes = [ec]' + R*[xe;ye];

Xeact = [q(1)+L1*cos(q(6))+L2*cos(q(6)+q(7));
      q(2)+L1*sin(q(6))+L2*sin(q(6)+q(7))];
  
[qd(4:5,1)] = qtemp(1:2);
[qd(8:9,1)] = qtemp(3:4);

K = [1/(5*t+1) 0;
     0 1/(10*t+1)];

qdtemp = pinv(J)*(dXed + K*(Xedes-Xeact));

[qd(6:7,1)] = qdtemp(3:4);

Jb = [(r/2)*cos(q(3)) (r/2)*cos(q(3));
      (r/2)*sin(q(3)) (r/2)*sin(q(3));
      r/(2*bw) -r/(2*bw)];

[qd(1:3,1)] = Jb*qd(4:5,1);
end
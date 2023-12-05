function qd = JointSpaceSquare(t,q)

% q = [x , y , phi, theR, theL, the1, the2]
L1 = 1.5; L2 = 1;
r = 0.5;
bw = 0.6;

s = 3.5;
sc = [1 1];

ang = rad2deg((2*pi/20)*t);

ang = deg2rad(wrapTo360(ang));

if ang <= (pi/4)
    Xed = [0;
        s*sec(ang)^2]*(2*pi/20);
elseif ((pi/4) < ang) &&  (ang <= (3*pi/4))
    Xed = [s*csc(ang)^2;
        0]*(2*pi/20);
elseif ((3*pi/4) < ang) && (ang <= (5*pi/4))
    Xed = [0;
        -s*sec(ang)^2]*(2*pi/20);
elseif ((5*pi/4) < ang) && (ang <= (7*pi/4))
    Xed = [-s*csc(ang)^2;
        0]*(2*pi/20);
elseif ((7*pi/4) < ang) && (ang <= (2*pi))
    Xed = [0;
        s*sec(ang)^2]*(2*pi/20);
end

J = [r/2*cos(q(3)) r/2*cos(q(3)) -L1*sin(q(6)) -L2*sin(q(7));
     r/2*sin(q(3)) r/2*sin(q(3)) L1*cos(q(6)) L2*cos(q(7))];

qtemp = pinv(J)*Xed;

dv = [0;0.66*2*qtemp(2);2*qtemp(3);0.25*2*qtemp(4)];

qdes = qtemp + (eye(4) - J\J)*dv;

[qd(4:5,1)] = qdes(1:2);
[qd(8:9,1)] = qdes(3:4);

K = (1/3)*eye(2);

[qd(6:7,1)] = qdes(3:4) + K*(q(8:9)-q(6:7));

Jb = [(r/2)*cos(q(3)) (r/2)*cos(q(3));
      (r/2)*sin(q(3)) (r/2)*sin(q(3));
      r/(2*bw) -r/(2*bw)];

[qd(1:3,1)] = Jb*qd(4:5,1);
end
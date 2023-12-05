function qint = IniInvKinP(x)
L1 = 1.5; L2 = 1;

t1 = deg2rad(15);
t2 = deg2rad(30);

phi = deg2rad(30);

R = [cos(phi), -sin(phi)
     sin(phi), cos(phi)];

xc = x - R*[L1*cos(t1)+L2*cos(t1+t2);
            L1*sin(t1)+L2*sin(t1+t2)];
      
qint = [xc(1);xc(2);deg2rad(30);deg2rad(5);deg2rad(5);t1;t2];
end
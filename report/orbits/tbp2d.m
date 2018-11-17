function xdot=tbp2d(x,t)
global m1 m2 m3 alpha;
xdot=zeros(12,1);
q1  = x(1:2);
q2  = x(3:4);
q3  = x(5:6);
p   = x(7:12);

q12 = q1-q2;
q13 = q1-q3;
q23 = q2-q3;
aux12 = alpha/(norm(q12))^(alpha+2);
aux13 = alpha/(norm(q13))^(alpha+2);
aux23 = alpha/(norm(q23))^(alpha+2);

xdot(1:6)   =  p;
xdot(7:8) = -m2*aux12*q12-m3*aux13*q13;
xdot(9:10) =  m1*aux12*q12-m3*aux23*q23;
xdot(11:12) =  m2*aux23*q23+m1*aux13*q13;

endfunction;

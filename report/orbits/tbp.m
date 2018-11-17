function xdot=tbp(x,t)
global m1 m2 m3 alpha;
xdot=zeros(18,1);
q1  = x(1:3);
q2  = x(4:6);
q3  = x(7:9);
p   = x(10:18);

q12 = q1-q2;
q13 = q1-q3;
q23 = q2-q3;

aux12 = alpha/(norm(q12))^(alpha+2);
aux13 = alpha/(norm(q13))^(alpha+2);
aux23 = alpha/(norm(q23))^(alpha+2);

xdot(1:9)   =  p; 
xdot(10:12) = -m2*aux12*q12-m3*aux13*q13;
xdot(13:15) =  m1*aux12*q12-m3*aux23*q23;
xdot(16:18) =  m2*aux23*q23+m1*aux13*q13;

endfunction;

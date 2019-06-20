l = 0.13;
m = l*(cos(q1) - 1 - cos(q2));
n = l*(sin(q1) - sin(q2));
% A = l^2 - (m^2 + n^2 + l^2);
% C9 = 2*l*sqrt(m^2 + n^2);
% Beta1 = atan2(m/n);
% q3 = atan2(A/C9)-Beta1;
% q4 = atan2((n+l*sin(q3))/(m+l*cos(q3)));


xFlag = 0;
yFlag = 0;
if q1 <= 0
q1 = 0.1*pi;
end
if q2 <= 0
q2 = 0.1*pi;
end
if q2 < q1
q2 = q1;
end
xb = l*cos(q2);
yb = l*sin(q2);
xd = l + l*cos(q1);
yd = l*sin(q1);
dsquare = (xb - xd)^2 + (yb - yd)^2;
%avoiding singularity,
if dsquare >= (2*l)^2
dsquare = (2*l)^2 - l/20;
x = sqrt(2*(l^2) + 2*(l^2)*cos(q1));
alpha = asin((l*sin(q1))/x);
beta = acos(((x^2) - (3*(l^2)))/(2*l*x));
%we are assuming that q2 is being driven by q1, which can also be the
%other way around.
q2 = alpha + beta;
xb = l*cos(q2);
yb = l*sin(q2);
slope = (sin(q1) - sin(q2))/(1 + cos(q1) - cos(q2));
q4 = atan2(slope, 1);
q3 = pi + q4;
%for x and y,
x = xb + l*cos(q4);
y = yb + l*sin(q4); % its actually yb - l*sin(q4) but since q4 is already denoted as negative
dq4_1 = -cot(q1)/((sec(q4))^2);
dq4_2 = -cot(q2)/((sec(q4))^2);
dx_1 = -l*sin(q4)*dq4_1;
dx_2 = -l*sin(q2) -l*sin(q4)*dq4_2;
dy_1 = l*cos(q4)*dq4_1;
dy_2 = l*cos(q2) + l*cos(q4)*dq4_2;
q4d = (cos(q1)*q1d - cos(q2)*q2d)/(2*cos(q4));
q4dd = (cos(q1)*q1dd - cos(q2)*q2dd + sin(q2)*(q2d^2) - sin(q1)*(q1d^2) + 2*sin(q4)*(q4d^2)); 
q3d = q4d;
q3dd = q4dd;
else
%by Heron's formula
k = 0.25 * sqrt(((2*l)^2 - dsquare) * (dsquare));
m = sqrt((2*l)^2 - dsquare);
%considering only the y-coordinate that is above the line joining the
%centers of the corresponding two circles.
x = (0.5 * (xb + xd)) + (2 * (yb - yd) * (k/dsquare));
y = (0.5 * (yb + yd)) + (2 * (xd - xb) * (k/dsquare));
%establishing constraints:
%for x,
if x >= 2*l + l*cos(q1)
x = 2*l + l*cos(q1);
xFlag = 1;
elseif x <= l*cos(q1)
x = l*cos(q1);
xFlag = 1;
end
%for y,
if y >= l + l*sin(q2)
y = l + l*sin(q2);
yFlag = 1;
elseif y <= l*sin(q2) - l
y = l*sin(q2) -l;
yFlag = 1;
end
%now for computing q3 and q4,
q4 = asin((y - (l*sin(q2)))/l);
q3 = acos((x - l - l*cos(q1))/l);
d = sqrt(dsquare);
ddq_1 = ((l^2)/d)*(sin(q1 - q2) - sin(q1));
ddq_2 = ((l^2)/d)*(sin(q2) - sin(q1 - q2));
dL_1 = 2*(l^2)*(sin(q1) - sin(q1 - q2));
dL_2 = 2*(l^2)*(sin(q1 - q2) - sin(q2));
dK_1 = ((d/(2*m))*dL_1 - m*ddq_1) / (4*dsquare);
dK_2 = ((d/(2*m))*dL_2 - m*ddq_2) / (4*dsquare);
if xFlag == 0
dx_1 = -0.5*l*sin(q1) - 2* (k/dsquare) *l*cos(q1) + 2*l*(sin(q2) - sin(q1))*dK_1;
dx_2 = -0.5*l*sin(q2) + 2* (k/dsquare) *l*cos(q2) + 2*l*(sin(q2) - sin(q1))*dK_2;
else
dx_1 = -l*sin(q1);
dx_2 = 0;
end
if yFlag == 0
dy_1 = 0.5*l*cos(q1) - 2* (k/dsquare) *l*sin(q1) + 2*l*(1 + cos(q1) - cos(q2))*dK_1;
dy_2 = 0.5*l*cos(q2) + 2* (k/dsquare) *l*sin(q2) + 2*l*(1 + cos(q1) - cos(q2))*dK_2;
else
dy_1 = 0;
dy_2 = l*cos(q2);
end
%for q3d and q4d,
a = sin(q2)*q2d - sin(q1)*q1d;
b = cos(q1)*q1d - cos(q2)*q2d;
c = sin(q3)*cos(q4) - sin(q4)*cos(q3);
q3d = (cos(q4)*a + sin(q4)*b)/c;
q4d = (sin(q3)*q3d + sin(q1)*q1d - sin(q2)*q2d)/sin(q4);
e = cos(q1)*(q1d^2) - cos(q2)*(q2d^2) + cos(q3)*(q3d^2) - cos(q4)*(q4d^2);
f = sin(q1)*(q1d^2) - sin(q2)*(q2d^2) + sin(q3)*(q3d^2) - sin(q4)*(q4d^2);
g = (cos(q1)*sin(q4) - cos(q4)*sin(q1)) * q1dd;
h = (cos(q2)*sin(q4) - cos(q4)*sin(q2)) * q2dd;
q3dd = (g - h - sin(q4)*f - cos(q4)*e) / (cos(q4)*sin(q3) - sin(q4)*cos(q3));
q4dd = (- f - cos(q2)*q2dd + cos(q1)*q1dd + cos(q3)*q3dd) / (cos(q4));
end

%q5 = 4*pi - (q1 + q2 + q3 + q4);
q = [q1; q2; q3; q4];
qd = [q1d;q2d;q3d;q4d];
qdd = [q1dd;q2dd;q3dd;q4dd];

s1 = sin(q1);
s2 = sin(q2);
s3 = sin(q3);
s4 = sin(q4);
c1 = cos(q1);
c2 = cos(q2);
c3 = cos(q3);
c4 = cos(q4);
M = [0.491777E-3+0.555539E-3.*c1.^2+0.555539E-3.*s1.^2,0.E-307,  0.E-307+0.277583E-3.*c1.*c3+0.277583E-3.*s1.*s3,0.E-307;0.E-307,   0.491777E-3+0.517514E-3.*c2.^2+0.517514E-3.*s2.^2,0.E-307,0.E-307+   0.25857E-3.*c2.*c4+0.25857E-3.*s2.*s4;0.E-307+0.277583E-3.*c1.*c3+  0.277583E-3.*s1.*s3,0.E-307,0.228464E-3+0.138791E-3.*c3.^2+   0.138791E-3.*s3.^2,0.E-307;0.E-307,0.E-307+0.25857E-3.*c2.*c4+   0.25857E-3.*s2.*s4,0.E-307,0.235023E-3+0.129285E-3.*c4.^2+  0.129285E-3.*s4.^2];
C = [0.E-307,0,0.E-307+0.277583E-3.*c3.*q3d.*s1+(-0.277583E-3).*c1.* q3d.*s3,0;0,0.E-307,0,0.E-307+0.25857E-3.*c4.*q4d.*s2+(   -0.25857E-3).*c2.*q4d.*s4;0.E-307+(-0.277583E-3).*c3.*q1d.*s1+   0.277583E-3.*c1.*q1d.*s3,0,0.E-307,0;0,0.E-307+(-0.25857E-3).*c4.*   q2d.*s2+0.25857E-3.*c2.*q2d.*s4,0,0.E-307];
JEta = [0.13E0.*s1,(-0.13E0).*s2,0.13E0.*s3,(-0.13E0).*s4;(-0.13E0).*c1,   0.13E0.*c2,(-0.13E0).*c3,0.13E0.*c4];
JEtad = [0.13E0.*c1.*q1d,(-0.13E0).*c2.*q2d,0.13E0.*c3.*q3d,(-0.13E0).*   c4.*q4d;0.13E0.*q1d.*s1,(-0.13E0).*q2d.*s2,0.13E0.*q3d.*s3,(   -0.13E0).*q4d.*s4];
Minv = inv(M);
A = JEta*Minv*JEta';
Ainv = inv(A);
B = diag(ones(4,1)) - JEta'*Ainv*JEta*Minv;
Binv = pinv(B);
Q = Binv*(M*qdd + JEta'*Ainv*JEtad*qd) + C*qd;
JEtaphi = JEta(:,3:4);
JEtatheta = JEta(:,1:2);
Jphitheta = -inv(JEtaphi)*JEtatheta;
Jqtheta = [diag(ones(2,1));Jphitheta];
Tau = Jqtheta'*Q;

function [outputset] = frontright4(t,state,n,tdindex,actibreak)
global xffrontright yffrontright omega;
if tdindex==0
    outputset(1:4,1)=0;
    return;
end
datafile;
%t = t - floor(n/2)*2*pi/omega - pi/omega;
%t = t + (t<0)*2*pi/omega;
n=1;

x=state(1);     xdot=state(2);
y=state(3);     ydot=state(4);
theta=state(5); thetadot=state(6);



% front leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta
qxfront=x-xffrontright-((-1)^n)*d1front*cos(theta)-d2front*sin(theta);
qyfront=y-yffrontright-((-1)^n)*d1front*sin(theta)+d2front*cos(theta);
qbodyfront=[cos(theta) sin(theta); -sin(theta) cos(theta)]*[qxfront; qyfront];
qxbodyfront=qbodyfront(1,1);
qybodyfront=qbodyfront(2,1);
qmagfront=sqrt(qxfront^2+qyfront^2);
gamma2front=acos((l1front^2+l2front^2-qmagfront^2)/(2*l1front*l2front));
phi2front=acos((l1front^2+qmagfront^2-l2front^2)/(2*l1front*qmagfront));
if (((-1)^n)*qxbodyfront)>=0
    gamma1front=pi/2+((-1)^n)*atan(qybodyfront/qxbodyfront)-phi2front; 
elseif  ((((-1)^n)*qxbodyfront)<0 & (qybodyfront>0))
    gamma1front=3*pi/2+((-1)^n)*atan(qybodyfront/qxbodyfront)-phi2front;
elseif  ((((-1)^n)*qxbodyfront)<0 & (qybodyfront<=0))
    gamma1front=-pi/2+((-1)^n)*atan(qybodyfront/qxbodyfront)-phi2front;
end
phi1front=acos((l2front^2+qmagfront^2-l1front^2)/(2*l2front*qmagfront));
psudeltafront=gamma2front-gamma1front-((-1)^n)*theta;

% Angular rates
gamma1dotfront=((-1)^n)*sin(psudeltafront)*xdot/(l1front*sin(gamma2front)) + ...
    cos(psudeltafront)*ydot/(l1front*sin(gamma2front)) + (-((-1)^n) + ...
    (cos(psudeltafront)/(l1front*sin(gamma2front)))*(-((-1)^n)*d1front*cos(theta) - d2front*sin(theta)) ...
    + (((-1)^n)*sin(psudeltafront)/(l1front*sin(gamma2front)))*(((-1)^n)*d1front*sin(theta) - ...
    d2front*cos(theta)))*thetadot;
gamma2dotfront=(((-1)^n)*xdot/sin(gamma2front))*(sin(gamma2front-psudeltafront)/l2front + ...
    sin(psudeltafront)/l1front) + (cos(psudeltafront)/l1front - ...
    cos(gamma2front - psudeltafront)/l2front)*ydot/sin(gamma2front) + ...
    (-((-1)^n)*(sin(gamma2front-psudeltafront)/l2front + ...
    sin(psudeltafront)/l1front)*(d2front*cos(theta)-((-1)^n)*d1front*sin(theta)) ...
    - ((-1)^n)*(cos(psudeltafront)/l1front - cos(gamma2front - psudeltafront)/l2front)*...
    (d1front*cos(theta)+((-1)^n)*d2front*sin(theta)))*thetadot/sin(gamma2front);

% Muscle forces and calculation of joint torque
tau1front=hiptorquevector(t,gamma1front,gamma1dotfront,1,actibreak(1:2));
tau2front=kneetorquevector(t,gamma2front,gamma2dotfront,1,actibreak(3:4));

% Get Fx,Fy
Fxfront=(((-1)^n)/sin(gamma2front))*(tau2front*sin(gamma2front-psudeltafront)/l2front+(tau1front+tau2front)*sin(psudeltafront)/l1front);
Fyfront=(1/sin(gamma2front))*(-tau2front*cos(gamma2front-psudeltafront)/l2front+(tau1front+tau2front)*cos(psudeltafront)/l1front);
outputset(1,1)=Fxfront;
outputset(2,1)=Fyfront;
outputset(3,1)=tau1front;
outputset(4,1)=tau2front;



end
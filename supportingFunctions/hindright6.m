function [outputset] = hindright6(t,state,n,tdindex,actibreak)
global xfrearright yfrearright omega;
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


% hind leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta
qxrear=x-xfrearright-((-1)^n)*d1rear*cos(theta)-d2rear*sin(theta);
qyrear=y-yfrearright-((-1)^n)*d1rear*sin(theta)+d2rear*cos(theta);
qbodyrear=[cos(theta) sin(theta); -sin(theta) cos(theta)]*[qxrear; qyrear];
qxbodyrear=qbodyrear(1,1);
qybodyrear=qbodyrear(2,1);
qmagrear=sqrt(qxrear^2+qyrear^2);
gamma2rear=acos((l1rear^2+l2rear^2-qmagrear^2)/(2*l1rear*l2rear));
phi2rear=acos((l1rear^2+qmagrear^2-l2rear^2)/(2*l1rear*qmagrear));
if (((-1)^n)*qxbodyrear)>=0
    gamma1rear=pi/2+((-1)^n)*atan(qybodyrear/qxbodyrear)-phi2rear; 
elseif  ((((-1)^n)*qxbodyrear)<0 & (qybodyrear>0))
    gamma1rear=3*pi/2+((-1)^n)*atan(qybodyrear/qxbodyrear)-phi2rear;
elseif  ((((-1)^n)*qxbodyrear)<0 & (qybodyrear<=0))
    gamma1rear=-pi/2+((-1)^n)*atan(qybodyrear/qxbodyrear)-phi2rear;
end
phi1rear=acos((l2rear^2+qmagrear^2-l1rear^2)/(2*l2rear*qmagrear));
psudeltarear=gamma2rear-gamma1rear-((-1)^n)*theta;

% Angular rates
gamma1dotrear=((-1)^n)*sin(psudeltarear)*xdot/(l1rear*sin(gamma2rear)) + ...
    cos(psudeltarear)*ydot/(l1rear*sin(gamma2rear)) + (-((-1)^n) + ...
    (cos(psudeltarear)/(l1rear*sin(gamma2rear)))*(-((-1)^n)*d1rear*cos(theta) - d2rear*sin(theta)) ...
    + (((-1)^n)*sin(psudeltarear)/(l1rear*sin(gamma2rear)))*(((-1)^n)*d1rear*sin(theta) - ...
    d2rear*cos(theta)))*thetadot;
gamma2dotrear=(((-1)^n)*xdot/sin(gamma2rear))*(sin(gamma2rear-psudeltarear)/l2rear + ...
    sin(psudeltarear)/l1rear) + (cos(psudeltarear)/l1rear - ...
    cos(gamma2rear - psudeltarear)/l2rear)*ydot/sin(gamma2rear) + ...
    (-((-1)^n)*(sin(gamma2rear-psudeltarear)/l2rear + ...
    sin(psudeltarear)/l1rear)*(d2rear*cos(theta)-((-1)^n)*d1rear*sin(theta)) ...
    - ((-1)^n)*(cos(psudeltarear)/l1rear - cos(gamma2rear - psudeltarear)/l2rear)*...
    (d1rear*cos(theta)+((-1)^n)*d2rear*sin(theta)))*thetadot/sin(gamma2rear);

% Muscle forces and calculation of joint torques
tau1rear=hiptorquevector(t,gamma1rear,gamma1dotrear,3,actibreak(1:2));
tau2rear=kneetorquevector(t,gamma2rear,gamma2dotrear,3,actibreak(3:4));

% Get Fx,Fy
Fxrear=(((-1)^n)/sin(gamma2rear))*(tau2rear*sin(gamma2rear-psudeltarear)/l2rear+(tau1rear+tau2rear)*sin(psudeltarear)/l1rear);
Fyrear=(1/sin(gamma2rear))*(-tau2rear*cos(gamma2rear-psudeltarear)/l2rear+(tau1rear+tau2rear)*cos(psudeltarear)/l1rear);
outputset(1,1)=Fxrear;
outputset(2,1)=Fyrear;
outputset(3,1)=tau1rear;
outputset(4,1)=tau2rear;


end

function [outputset] = getvarmusc(t,state,n,acti)

%global xffront yffront xfmid yfmid xfrear yfrear;

global xffrontleft yffrontleft xfmidright yfmidright xfrearleft yfrearleft;
global xffrontright yffrontright xfmidleft yfmidleft xfrearright yfrearright;

xffront=xffrontleft*(mod(n,2)==0)+xffrontright*(mod(n,2)==1);
yffront=yffrontleft*(mod(n,2)==0)+yffrontright*(mod(n,2)==1);
xfmid=xfmidleft*(mod(n,2)==1)+xfmidright*(mod(n,2)==0);
yfmid=yfmidleft*(mod(n,2)==1)+yfmidright*(mod(n,2)==0);
xfrear=xfrearleft*(mod(n,2)==0)+xfrearright*(mod(n,2)==1);
yfrear=yfrearleft*(mod(n,2)==0)+yfrearright*(mod(n,2)==1);


datamusc;

x=state(1);     xdot=state(2);
y=state(3);     ydot=state(4);
theta=state(5); thetadot=state(6);


% front leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta

qxfront=x-xffront-((-1)^n)*d1front*cos(theta)-d2front*sin(theta);
qyfront=y-yffront-((-1)^n)*d1front*sin(theta)+d2front*cos(theta);

%qxfront=comXvector-xfootfrontvector-((-1).^stancecount).*d1front.*cos(thetavector)-d2front.*sin(thetavector);
%qyfront=comYvector-yfootfrontvector-((-1).^stancecount).*d1front.*sin(thetavector)+d2front.*cos(thetavector);

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
tau1front=hiptorquevector(t,gamma1front,gamma1dotfront,1,acti(1:2,1));
tau2front=kneetorquevector(t,gamma2front,gamma2dotfront,1,acti(3:4,1));

% Get Fx,Fy
Fxfront=(((-1)^n)/sin(gamma2front))*(tau2front*sin(gamma2front-psudeltafront)/l2front+(tau1front+tau2front)*sin(psudeltafront)/l1front);
Fyfront=(1/sin(gamma2front))*(-tau2front*cos(gamma2front-psudeltafront)/l2front+(tau1front+tau2front)*cos(psudeltafront)/l1front);
outputset(1,1)=Fxfront;
outputset(2,1)=Fyfront;
outputset(7,1)=tau1front;
outputset(8,1)=tau2front;

% middle leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta
qxmid=x-xfmid-((-1)^(n+1))*d1mid*cos(theta)-d2mid*sin(theta);
qymid=y-yfmid-((-1)^(n+1))*d1mid*sin(theta)+d2mid*cos(theta);
qbodymid=[cos(theta) sin(theta); -sin(theta) cos(theta)]*[qxmid; qymid];
qxbodymid=qbodymid(1,1);
qybodymid=qbodymid(2,1);
qmagmid=sqrt(qxmid^2+qymid^2);
gamma2mid=acos((l1mid^2+l2mid^2-qmagmid^2)/(2*l1mid*l2mid));
phi2mid=acos((l1mid^2+qmagmid^2-l2mid^2)/(2*l1mid*qmagmid));
if  (((-1)^(n+1))*qxbodymid)>=0
    gamma1mid=pi/2+((-1)^(n+1))*atan(qybodymid/qxbodymid)-phi2mid;
elseif  ((((-1)^(n+1))*qxbodymid)<0 & (qybodymid>0))
    gamma1mid=3*pi/2+((-1)^(n+1))*atan(qybodymid/qxbodymid)-phi2mid;
elseif  ((((-1)^(n+1))*qxbodymid)<0 & (qybodymid<=0))
    gamma1mid=-pi/2+((-1)^(n+1))*atan(qybodymid/qxbodymid)-phi2mid;   
end
phi1mid=acos((l2mid^2+qmagmid^2-l1mid^2)/(2*l2mid*qmagmid));
psudeltamid=gamma2mid-gamma1mid-((-1)^(n+1))*theta;

% Angular rates
gamma1dotmid=-((-1)^n)*sin(psudeltamid)*xdot/(l1mid*sin(gamma2mid)) + ...
    cos(psudeltamid)*ydot/(l1mid*sin(gamma2mid)) + (((-1)^n) + ...
    (cos(psudeltamid)/(l1mid*sin(gamma2mid)))*(((-1)^n)*d1mid*cos(theta) - d2mid*sin(theta)) ...
    + (-((-1)^n)*sin(psudeltamid)/(l1mid*sin(gamma2mid)))*(-((-1)^n)*d1mid*sin(theta) - ...
    d2mid*cos(theta)))*thetadot;
gamma2dotmid=(-((-1)^n)*xdot/sin(gamma2mid))*(sin(gamma2mid-psudeltamid)/l2mid + ...
    sin(psudeltamid)/l1mid) + (cos(psudeltamid)/l1mid - ...
    cos(gamma2mid - psudeltamid)/l2mid)*ydot/sin(gamma2mid) + ...
    (((-1)^n)*(sin(gamma2mid-psudeltamid)/l2mid + ...
    sin(psudeltamid)/l1mid)*(d2mid*cos(theta)+((-1)^n)*d1mid*sin(theta)) ...
    + ((-1)^n)*(cos(psudeltamid)/l1mid - cos(gamma2mid - psudeltamid)/l2mid)*...
    (d1mid*cos(theta)-((-1)^n)*d2mid*sin(theta)))*thetadot/sin(gamma2mid);

% Muscle forces and calculation of joint torques
tau1mid=hiptorquevector(t,gamma1mid,gamma1dotmid,2,acti(5:6,1));
tau2mid=kneetorquevector(t,gamma2mid,gamma2dotmid,2,acti(7:8,1));

% Get Fx,Fy
Fxmid=(((-1)^(n+1))/sin(gamma2mid))*(tau2mid*sin(gamma2mid-psudeltamid)/l2mid+(tau1mid+tau2mid)*sin(psudeltamid)/l1mid);
Fymid=(1/sin(gamma2mid))*(-tau2mid*cos(gamma2mid-psudeltamid)/l2mid+(tau1mid+tau2mid)*cos(psudeltamid)/l1mid);
outputset(3,1)=Fxmid;
outputset(4,1)=Fymid;
outputset(9,1)=tau1mid;
outputset(10,1)=tau2mid;

% hind leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta
qxrear=x-xfrear-((-1)^n)*d1rear*cos(theta)-d2rear*sin(theta);
qyrear=y-yfrear-((-1)^n)*d1rear*sin(theta)+d2rear*cos(theta);
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
tau1rear=hiptorquevector(t,gamma1rear,gamma1dotrear,3,acti(9:10,1));
tau2rear=kneetorquevector(t,gamma2rear,gamma2dotrear,3,acti(11:12,1));

% Get Fx,Fy
Fxrear=(((-1)^n)/sin(gamma2rear))*(tau2rear*sin(gamma2rear-psudeltarear)/l2rear+(tau1rear+tau2rear)*sin(psudeltarear)/l1rear);
Fyrear=(1/sin(gamma2rear))*(-tau2rear*cos(gamma2rear-psudeltarear)/l2rear+(tau1rear+tau2rear)*cos(psudeltarear)/l1rear);
outputset(5,1)=Fxrear;
outputset(6,1)=Fyrear;
outputset(11,1)=tau1rear;
outputset(12,1)=tau2rear;


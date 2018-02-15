function [outputset] = midleft5(t,state,n,tdindex,actibreak)
global xfmidleft yfmidleft omega;
if tdindex==0
    outputset(1:4,1)=0;
    return;
end
datafile;
%t = t - floor(n/2)*2*pi/omega - pi/omega; % This is right since completely feedforward
%t = t + (t<0)*2*pi/omega;
n=1;

x=state(1);     xdot=state(2);
y=state(3);     ydot=state(4);
theta=state(5); thetadot=state(6);

% middle leg
% Get gamma1,gamma2,qx,qy,qmag,phi1,phi2,psudelta
qxmid=x-xfmidleft-((-1)^(n+1))*d1mid*cos(theta)-d2mid*sin(theta);
qymid=y-yfmidleft-((-1)^(n+1))*d1mid*sin(theta)+d2mid*cos(theta);
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
tau1mid=hiptorquevector(t,gamma1mid,gamma1dotmid,2,actibreak(1:2));
tau2mid=kneetorquevector(t,gamma2mid,gamma2dotmid,2,actibreak(3:4));

% Get Fx,Fy
Fxmid=(((-1)^(n+1))/sin(gamma2mid))*(tau2mid*sin(gamma2mid-psudeltamid)/l2mid+(tau1mid+tau2mid)*sin(psudeltamid)/l1mid);
Fymid=(1/sin(gamma2mid))*(-tau2mid*cos(gamma2mid-psudeltamid)/l2mid+(tau1mid+tau2mid)*cos(psudeltamid)/l1mid);
outputset(1,1)=Fxmid;
outputset(2,1)=Fymid;
outputset(3,1)=tau1mid;
outputset(4,1)=tau2mid;

end


global gam1footfront gam1footmid gam1footrear;
global gam2footfront gam2footmid gam2footrear;
global Vchangefront Vchangemid Vchangerear;


bodylength = 0.044;

% Foot TD
footplacement = [footModel(0,1,4); footModel(pi,2,4); footModel(0,3,4); 
    footModel(pi,4,4); footModel(0,5,4); footModel(pi,6,4)];
dataplace = [real(footplacement).'; imag(footplacement).'];
modelplace = [0.02 0.007 -0.01; -0.011 0.013 -0.013]./bodylength;

xbfootfront=mean([-dataplace(2,1) dataplace(2,6)])*bodylength - 0.007; % - 0.006, - 0.008
ybfootfront=mean([dataplace(1,1) dataplace(1,6)])*bodylength; % - 0.004
xbfootmid=mean([dataplace(2,2) -dataplace(2,5)])*bodylength;
ybfootmid=mean([dataplace(1,2) dataplace(1,5)])*bodylength;
xbfootrear=mean([-dataplace(2,3) dataplace(2,4)])*bodylength;
ybfootrear=mean([dataplace(1,3) dataplace(1,4)])*bodylength;

footplace = [xbfootfront xbfootmid xbfootrear; ybfootfront ybfootmid ybfootrear]./bodylength;

% xbfootfront=-0.011;
% ybfootfront=0.02; % 0.03
% xbfootmid=0.013;
% ybfootmid=0.007; % 0.017
% xbfootrear=-0.013;
% ybfootrear=-0.01;

datamusc;

qxfootfront = xbfootfront + d1front;
qyfootfront = ybfootfront - d2front;
qxfootmid = xbfootmid - d1mid;
qyfootmid = ybfootmid - d2mid;
qxfootrear = xbfootrear + d1rear;
qyfootrear = ybfootrear - d2rear;

qmagfootfront = sqrt(qxfootfront^2 + qyfootfront^2);
qmagfootmid = sqrt(qxfootmid^2 + qyfootmid^2);
qmagfootrear = sqrt(qxfootrear^2 + qyfootrear^2);

gam2footfront = acos( ( l1front^2 + l2front^2 - qmagfootfront^2 )/( 2*l1front*l2front ) );
gam2footmid = acos( ( l1mid^2 + l2mid^2 - qmagfootmid^2 )/( 2*l1mid*l2mid ) );
gam2footrear = acos( ( l1rear^2 + l2rear^2 - qmagfootrear^2 )/( 2*l1rear*l2rear ) );

phi2footfront = acos( ( l1front^2 + qmagfootfront^2 - l2front^2 )/( 2*l1front*qmagfootfront ) );
phi2footmid = acos( ( l1mid^2 + qmagfootmid^2 - l2mid^2 )/( 2*l1mid*qmagfootmid ) );
phi2footrear = acos( ( l1rear^2 + qmagfootrear^2 - l2rear^2 )/( 2*l1rear*qmagfootrear ) );

gam1footfront = pi/2 + atan(qyfootfront/qxfootfront) - phi2footfront;
gam1footmid = pi/2 + atan(qyfootmid/qxfootmid) - phi2footmid;
gam1footrear = pi/2 + atan(qyfootrear/qxfootrear) - phi2footrear;

Vchangefront = -(0.1 - gam1footfront)/2 + 0.25; % 1.75
Vchangemid = -(0.15 - gam1footmid)/1.5 + 0.25; % 2.25
Vchangerear = -(0.15 - gam1footrear)/1.28 + 0.25; % 2.5


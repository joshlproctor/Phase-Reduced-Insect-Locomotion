function fun = evathetamusc(input,plotfig,output)
global Vdes omega rightthetainit rightthetadotinit;
datamusc;
% Vdes=0.25;
% omega=2*pi*12;

[frontxf frontyf midxf midyf rearxf rearyf]=footpP;
% if Vdes<0.25
%     frontxf=xbfootfront; frontyf=ybfootfront;
%     midxf=xbfootmid; midyf=ybfootmid;
%     rearxf=xbfootrear; rearyf=ybfootrear;
% else
%     frontxf=xbfootfront+0.01*(Vdes-0.25);
%     frontyf=ybfootfront+0.03*(Vdes-0.25);
%     midxf=xbfootmid+0.02*(Vdes-0.25);
%     midyf=ybfootmid+0.04*(Vdes-0.25);
%     rearxf=xbfootrear-0.01*(Vdes-0.25);
%     rearyf=ybfootrear+0.045*(Vdes-0.25);
% end

n=0;
xinit=0;    yinit=0;
thetainit=input(1,1);
thetadotinit=input(1,2);
xfootfr=((-1)^n)*frontxf*cos(thetainit) - frontyf*sin(thetainit);     % w.r.t. the starting point of left stance   
yfootfr=((-1)^n)*frontxf*sin(thetainit) + frontyf*cos(thetainit);     % w.r.t. the starting point of left stance
xfootm=((-1)^n)*midxf*cos(thetainit) - midyf*sin(thetainit);     % w.r.t. the starting point of left stance   
yfootm=((-1)^n)*midxf*sin(thetainit) + midyf*cos(thetainit);     % w.r.t. the starting point of left stance
xfoothind=((-1)^n)*rearxf*cos(thetainit) - rearyf*sin(thetainit);     % w.r.t. the starting point of left stance   
yfoothind=((-1)^n)*rearxf*sin(thetainit) + rearyf*cos(thetainit);     % w.r.t. the starting point of left stance
time=linspace(0,pi/omega,10);
x=-(Axrear/(m*omega^2)).*sin(omega.*time);                  % Not for the right stance
y=(Aymid/(4*m*omega^2)).*sin(2.*omega.*time) + Vdes.*time;  % Not for the right stance
Cone=-xfootfr*Ayfront+xfoothind*Ayrear-yfootfr*Axfront+yfootm*Axmid-yfoothind*Axrear;
lefttheta=(3*Axrear*Aymid/(8*m*I*omega^2)).*(cos(omega.*time)./(omega^2) - cos(3.*omega.*time)./(9*omega^2)) ...
    + xfootm.*Aymid.*sin(2.*omega.*time)./(4*I*omega^2) - Cone.*sin(omega.*time)./(I*omega^2) ...
    - (Vdes*Axrear/I).*(2.*cos(omega.*time)./omega^3 + time.*sin(omega.*time)./omega^2) ...
    + (thetadotinit - xfootm*Aymid/(2*I*omega) + Cone/(I*omega)).*time ...
    - Axrear*Aymid/(3*m*I*omega^4) + 2*Vdes*Axrear/(I*omega^3) + thetainit;
leftthetadot=(3*Axrear*Aymid/(8*m*I*omega^2)).*(sin(3.*omega.*time)./(3*omega) - sin(omega.*time)./omega) ...
    + xfootm.*Aymid.*cos(2.*omega.*time)./(2*I*omega) - Cone.*cos(omega.*time)./(I*omega) ...
    + (Vdes*Axrear/I).*(-time.*cos(omega.*time)./omega + sin(omega.*time)./omega^2) ...
    + thetadotinit - xfootm*Aymid/(2*I*omega) + Cone/(I*omega);        
leftmoment=-xfootm.*Aymid.*sin(2.*omega.*time)+Axrear.*Vdes.*time.*sin(omega.*time) ...
    -(3*Axrear*Aymid/(4*m*omega^2)).*sin(omega.*time).*sin(2.*omega.*time)+Cone.*sin(omega.*time);

if output==1
    rightthetainit=lefttheta(length(lefttheta));
    rightthetadotinit=leftthetadot(length(leftthetadot));
end
n=1;
xinit=x(length(x));
yinit=y(length(y));
thetainit=lefttheta(length(lefttheta));
thetadotinit=leftthetadot(length(leftthetadot));
xfootfr=((-1)^n)*frontxf*cos(thetainit) - frontyf*sin(thetainit);     % w.r.t. the starting point of right stance   
yfootfr=((-1)^n)*frontxf*sin(thetainit) + frontyf*cos(thetainit);     % w.r.t. the starting point of right stance
xfootm=((-1)^n)*midxf*cos(thetainit) - midyf*sin(thetainit);     % w.r.t. the starting point of right stance   
yfootm=((-1)^n)*midxf*sin(thetainit) + midyf*cos(thetainit);     % w.r.t. the starting point of right stance
xfoothind=((-1)^n)*rearxf*cos(thetainit) - rearyf*sin(thetainit);     % w.r.t. the starting point of right stance   
yfoothind=((-1)^n)*rearxf*sin(thetainit) + rearyf*cos(thetainit);     % w.r.t. the starting point of right stance
Ctwo=-xfootfr*Ayfront+xfoothind*Ayrear+yfootfr*Axfront-yfootm*Axmid+yfoothind*Axrear;
righttheta=-(3*Axrear*Aymid/(8*m*I*omega^2)).*(cos(omega.*time)./(omega^2) - cos(3.*omega.*time)./(9*omega^2)) ...
    + xfootm.*Aymid.*sin(2.*omega.*time)./(4*I*omega^2) - Ctwo.*sin(omega.*time)./(I*omega^2) ...
    + (Vdes*Axrear/I).*(2.*cos(omega.*time)./omega^3 + time.*sin(omega.*time)./omega^2) ...
    + (thetadotinit - xfootm*Aymid/(2*I*omega) + Ctwo/(I*omega)).*time ...
    + Axrear*Aymid/(3*m*I*omega^4) - 2*Vdes*Axrear/(I*omega^3) + thetainit;
rightthetadot=-(3*Axrear*Aymid/(8*m*I*omega^2)).*(sin(3.*omega.*time)./(3*omega) - sin(omega.*time)./omega) ...
    + xfootm.*Aymid.*cos(2.*omega.*time)./(2*I*omega) - Ctwo.*cos(omega.*time)./(I*omega) ...
    - (Vdes*Axrear/I).*(-time.*cos(omega.*time)./omega + sin(omega.*time)./omega^2) ...
    + thetadotinit - xfootm*Aymid/(2*I*omega) + Ctwo/(I*omega);        
rightmoment=-xfootm.*Aymid.*sin(2.*omega.*time)-Axrear.*Vdes.*time.*sin(omega.*time) ...
    +(3*Axrear*Aymid/(4*m*omega^2)).*sin(omega.*time).*sin(2.*omega.*time)+Ctwo.*sin(omega.*time);

theta=[lefttheta';righttheta'];
thetadot=[leftthetadot';rightthetadot'];
moment=[leftmoment';rightmoment'];
time=[time';(pi/omega + time)'];

thetafunction=fit(time,theta,'smoothingspline');
thetadotfunction=fit(time,thetadot,'smoothingspline');

thetainteg=integrate(thetafunction,2*pi/omega,0);
thetadotinteg=integrate(thetadotfunction,2*pi/omega,0);

if plotfig==1
    figure
    plot(time,theta,'b-');
    figure    
    plot(time,thetadot,'r-');
    figure 
    plot(time,moment,'k-');
end

fun=[thetainteg thetadotinteg];
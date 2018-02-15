function torque = kneetorquevector(time,jointangle,jointrate,num,actfun)
global omega vmax moment_arm tors_k cdamp highvmax Vdes muscSn FreqCheck SpringCheck;
global tors_k cdamp
% global pertarray randspike param expno;
neuraldata;

% c='kbcgmry'; 
vmax=75; moment_arm=0.4; 
%tors_k=2e-5; cdamp=4e-6; 
highvmax=75;      
%Vdes=0.26;
if Vdes<0.375
    omega=2*pi*(13+(Vdes-0.375)/0.0375);
else
    omega=2*pi*13;
end

if SpringCheck==1
    
elseif SpringCheck==0
     tors_k=2e-5;
    cdamp=4e-6;
end

vmax = vmax.*(Vdes<0.4).*(num==1) + highvmax.*(Vdes>=0.4).*(num==1) + ...
    vmax.*(num~=1);

%% Hill-type muscle
% par=param(num,:) + pertarray(num,:).*(randspike==1);

steep=1; swing=1;
time=time.*1000;

if swing==1
    extsig=1+exp( steep.*(time-1000*pi/omega) ).*(num~=1);
%     flexsig=1+exp( -steep.*(time-1000*pi/omega) );
    flexsig=1;
end

%% Activation

% Extensor calculations

exteact = actfun(1,1);

% Flexor calculations

flexact = actfun(2,1);


%% Length-dependent force 
rest_len=8.95; % In millimetres
% rest_angle=(0.37/0.0035)*(pi/180); % In radians
rest_angle = 1.4*(num==1) + 1.2*(num==2) + 1.1*(num==3);
% rest_len_179=4; % In millimetres
opt_len=1.11*rest_len; % Assumed from info about muscle 179
p1=-10.49;  p2=21.08;   p3=-10.48;  p4=0.8551; % Metathoracic 179 (Ahn et al., 2006)

% Parallel passive element
par_k=4;

% For the extensor
exteref = rest_len + moment_arm*rest_angle;
extelength = exteref - moment_arm.*jointangle; 
% extelength=rest_len + (0.37-0.0035.*jointangle.*(180/pi)).*rest_len_179.*moment_arm; % wierdly changed
x=extelength./opt_len;
exteFlen=p1.*x.^3 + p2.*x.^2 + p3.*x + p4;

exteparallforce=(x>1).*(par_k.*((x-1).^2));

% For the flexor
flexref = rest_len - moment_arm*rest_angle;
flexlength = flexref + moment_arm.*jointangle;
% attach_len = moment_arm*cos(rest_angle) + ...
%     sqrt( rest_len^2 - (moment_arm*sin(rest_angle))^2 );
% flexlength = sqrt( moment_arm^2 + attach_len^2 - ...
%     2.*moment_arm.*attach_len.*cos(jointangle) );
x=flexlength./opt_len;
flexFlen=p1.*x.^3 + p2.*x.^2 + p3.*x + p4;

flexparallforce=(x>1).*(par_k.*((x-1).^2));

%% Velocity-dependent force (177c)
% vmax=49; % In mm/s
c=2.4298; % Metathoracic (Ahn and Full, 2002)
f=1.6;  btilda=1/c;

% For the extensor - vel is defined as -ldot
extevelocity = moment_arm.*jointrate;
% extevelocity=0.0035.*rest_len_179.*jointrate.*(180/pi).*moment_arm; % wierdly changed 
vel=extevelocity;
exteFvel=(vel>=0).*((vmax-vel)./(vmax+c.*vel)) + (vel<0).*(((f-1)*vmax-f.*(1+btilda).*vel)./((f-1)*vmax-(1+btilda).*vel));
% For the flexor
flexvelocity = -moment_arm.*jointrate;
% flexvelocity=-moment_arm.*attach_len.*sin(jointangle).*jointrate./flexlength;
vel=flexvelocity;
flexFvel=(vel>=0).*((vmax-vel)./(vmax+c.*vel)) + (vel<0).*(((f-1)*vmax-f.*(1+btilda).*vel)./((f-1)*vmax-(1+btilda).*vel));


%% Final force calculations

% Isometric force
musclemass=21.2; % In mg for 177c (Ahn and Full, 2002)
area=musclemass*0.01/rest_len; % length in mm, mass in mg, area in cm^2, density assumed to be 1 gm/cm^3 
max_stress=25.6; % In N/cm^2
F0=max_stress*area; % In N




%%Define Vel dependent part of Force Scaling.
% if Vdes==.26
if FreqCheck==0
 par = muscS(3:4,1)*(num==1) + muscS(7:8,1)*(num==2) + muscS(11:12,1)*(num==3);
end
if FreqCheck==1
 par = muscSn(3:4,1)*(num==1) + muscSn(7:8,1)*(num==2) + muscSn(11:12,1)*(num==3);
end
%
 
 % end
% if Vdes~=.26
% 
% FrontmuscSe=muscS(3,1);
% FrontmuscSf=muscS(4,1);
% 
% MidmuscSe=muscS(7,1);
% MidmuscSf=muscS(8,1);
% 
% RearmuscSe=muscS(11,1);
% RearmuscSf=muscS(12,1);
% 
% par = [FrontmuscSe;FrontmuscSf]*(num==1) + [MidmuscSe;MidmuscSf]*(num==2) + [RearmuscSe;RearmuscSf]*(num==3);
% end

% Instantaneous muscle stress/force
% For the extensor - was 1
exteforceperunitarea=par(1,1).*max_stress.*(exteact.*exteFlen.*exteFvel + exteparallforce)./extsig;
% exteforceperunitarea=par(5).*max_stress.*(exteact.*exteFlen + exteparallforce);
exteforce=exteforceperunitarea.*area;
% For the flexor - was 0.2
flexforceperunitarea=par(2,1).*max_stress.*(flexact.*flexFlen.*flexFvel + flexparallforce)./flexsig;
% flexforceperunitarea=par(6).*max_stress.*(flexact.*flexFlen + flexparallforce);
flexforce=flexforceperunitarea.*area;
% theangle=acos( (flexlength.^2 + moment_arm^2 - attach_len^2)./ ...
%     (2.*flexlength.*moment_arm) );
% act_moment_arm=(moment_arm/1000).*sin(theangle); % Moment arm is in mm

% tors_k=2e-5;
passtorque = tors_k.*(rest_angle-jointangle) - cdamp.*jointrate;

% Joint torque
% torque=exteforce.*moment_arm.*0.001 - flexforce.*act_moment_arm;
torque=exteforce.*moment_arm.*0.001 - flexforce.*moment_arm.*0.001 + passtorque;

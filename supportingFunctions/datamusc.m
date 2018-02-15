% Roach
m=0.0025;  %0.0025  % Mass
I=2.04e-7;          % Moment of Inertia about COM

% front leg
l1front=0.015; % 0.0075, 0.012, 0.016, 0.013
l2front=0.016; % 0.0096, 0.009, 0.016, 0.017
d1front=0.001; % 0.0035
d2front=0.018; % 0.014 
% middle leg
l1mid=0.013; % 0.01, 0.008, 0.015, 0.015
l2mid=0.019; % 0.014, 0.014, 0.019, 0.023
d1mid=0.001; % 0.0035
d2mid=0.009; % 0.007 
% hind leg
l1rear=0.015; % 0.01, 0.016, 0.016
l2rear=0.026; % 0.02, 0.026, 0.028
d1rear=0.001; % 0.0035
d2rear=0; % -0.002 0.002 0.003

k1front=2e-4;
k2front=2e-4;
k1mid=2e-4;
k2mid=2e-4;
k1rear=2e-4;
k2rear=2e-4;

Axfront=0.0051;     %Foot reaction force in the x direction
Axmid=0.0051;
Axrear=0.0032;
Ayfront=0.0049;     %Foot reaction force in the y direction
Aymid=0.004;
Ayrear=0.0049;


% Axfront = Axfront*0.5;
% Axmid = Axmid*0.5;
% Axrear = Axrear*0.5;
% Ayfront = Ayfront*0.5;
% Aymid = Aymid*0.5;
% Ayrear = Ayrear*0.5;

% The protocol for keeping the feet forward for higher speeds affects
% stability. Scaling of the foot position for larger body lengths affects
% stability. Damping affects stability. Passive torsional stiffness affects
% stability. Arbitrary change in the front foot position affects stability.
% Change in thigh and shank lengths affect stability. Have to study the
% affect of changing hip positions on stability.
% Last time, the foot positions were changed to Shai's data in order to get
% the hind knee and hip torques positive.
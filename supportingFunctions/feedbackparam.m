

%Feedback parameters
EpostPhase_FBe=zeros(24,1);
EpostPhase_FBi=-70*ones(24,1);

Hip1=0;
Hip2=0;
Hip3=0;
Knee1=0;
Knee2=0;
Knee3=0;

% Hip1=0;
% Hip2=.4e-4;
% Hip3=-.6e-4;
% Knee1=1e-4;
% Knee2=-3e-4;
% Knee3=2.5e-4;
Hip1=-2.5e-5;
Hip2=5e-5;
Hip3=-3.5e-5;
Knee1=8e-5;
Knee2=3e-5;
Knee3=2e-5;

Hip=[Hip1;Hip2;Hip3;Hip1;Hip2;Hip3];
Knee=[Knee1;Knee2;Knee3;Knee1;Knee2;Knee3];

Hip1y1e=.1;
Hip1x1e=(6e-5-Hip1);
Hip1y2e=.4;
Hip1x2e=(8e-5-Hip1);

Hip2y1e=.1;
Hip2x1e=(15e-5-Hip2);
Hip2y2e=.4;
Hip2x2e=(20e-5-Hip2);

Hip3y1e=.1;
Hip3x1e=(3e-5-Hip3);
Hip3y2e=.4;
Hip3x2e=(5e-5-Hip3);

Hip1y1i=.2;
Hip1x1i=(6e-5-Hip1);
Hip1y2i=.8;
Hip1x2i=(8e-5-Hip1);

Hip2y1i=.2;
Hip2x1i=(15e-5-Hip2);
Hip2y2i=.8;
Hip2x2i=(20e-5-Hip2);

Hip3y1i=.2;
Hip3x1i=(3e-5-Hip3);
Hip3y2i=.8;
Hip3x2i=(5e-5-Hip3);

Knee1y1e=.1;
Knee1x1e=(16e-5-Knee1);
Knee1y2e=.4;
Knee1x2e=(3e-4-Knee1);

Knee2y1e=.1;
Knee2x1e=(11e-5-Knee2);
Knee2y2e=.4;
Knee2x2e=(1.8e-4-Knee2);

Knee3y1e=.1;
Knee3x1e=(10e-5-Knee3);
Knee3y2e=.4;
Knee3x2e=(2e-4-Knee3);

Knee1y1i=.2;
Knee1x1i=(16e-5-Knee1);
Knee1y2i=.8;
Knee1x2i=(3e-4-Knee1);

Knee2y1i=.2;
Knee2x1i=(11e-5-Knee2);
Knee2y2i=.8;
Knee2x2i=(1.8e-4-Knee2);

Knee3y1i=.2;
Knee3x1i=(10e-5-Knee3);
Knee3y2i=.8;
Knee3x2i=(2e-4-Knee3);

SlopeKnee1i=(Knee1y1i-Knee1y2i)/(Knee1x1i-Knee1x2i);
SlopeKnee2i=(Knee2y1i-Knee2y2i)/(Knee2x1i-Knee2x2i);
SlopeKnee3i=(Knee3y1i-Knee3y2i)/(Knee3x1i-Knee3x2i);
SlopeKnee1e=(Knee1y1e-Knee1y2e)/(Knee1x1e-Knee1x2e);
SlopeKnee2e=(Knee2y1e-Knee2y2e)/(Knee2x1e-Knee2x2e);
SlopeKnee3e=(Knee3y1e-Knee3y2e)/(Knee3x1e-Knee3x2e);

SlopeKneei=[SlopeKnee1i;SlopeKnee2i;SlopeKnee3i;SlopeKnee1i;SlopeKnee2i;SlopeKnee3i];
SlopeKneee=[SlopeKnee1e;SlopeKnee2e;SlopeKnee3e;SlopeKnee1e;SlopeKnee2e;SlopeKnee3e];


SlopeHip1i=(Hip1y1i-Hip1y2i)/(Hip1x1i-Hip1x2i);
SlopeHip2i=(Hip2y1i-Hip2y2i)/(Hip2x1i-Hip2x2i);
SlopeHip3i=(Hip3y1i-Hip3y2i)/(Hip3x1i-Hip3x2i);
SlopeHip1e=(Hip1y1e-Hip1y2e)/(Hip1x1e-Hip1x2e);
SlopeHip2e=(Hip2y1e-Hip2y2e)/(Hip2x1e-Hip2x2e);
SlopeHip3e=(Hip3y1e-Hip3y2e)/(Hip3x1e-Hip3x2e);

SlopeHipi=[SlopeHip1i;SlopeHip2i;SlopeHip3i;SlopeHip1i;SlopeHip2i;SlopeHip3i];
SlopeHipe=[SlopeHip1e;SlopeHip2e;SlopeHip3e;SlopeHip1e;SlopeHip2e;SlopeHip3e];


yintKnee1i=Knee1y1i-SlopeKnee1i*Knee1x1i;
yintKnee1e=Knee1y1e-SlopeKnee1e*Knee1x1e;
yintKnee2i=Knee2y1i-SlopeKnee2i*Knee1x1i;
yintKnee2e=Knee2y1e-SlopeKnee2e*Knee1x1e;
yintKnee3i=Knee3y1i-SlopeKnee3i*Knee1x1i;
yintKnee3e=Knee3y1e-SlopeKnee3e*Knee1x1e;

yintHip1i=Hip1y1i-SlopeHip1i*Hip1x1i;
yintHip1e=Hip1y1e-SlopeHip1e*Hip1x1e;
yintHip2i=Hip2y1i-SlopeHip2i*Hip1x1i;
yintHip2e=Hip2y1e-SlopeHip2e*Hip1x1e;
yintHip3i=Hip3y1i-SlopeHip3i*Hip1x1i;
yintHip3e=Hip3y1e-SlopeHip3e*Hip1x1e;

yintKneei=[yintKnee1i;yintKnee2i;yintKnee3i;yintKnee1i;yintKnee2i;yintKnee3i];
yintKneee=[yintKnee1e;yintKnee2e;yintKnee3e;yintKnee1e;yintKnee2e;yintKnee3e];

yintHipi=[yintHip1i;yintHip2i;yintHip3i;yintHip1i;yintHip2i;yintHip3i];
yintHipe=[yintHip1e;yintHip2e;yintHip3e;yintHip1e;yintHip2e;yintHip3e];



% SlopeHip1=(.2)/(6e-5);
% SlopeHip2=(.2)/(2.5e-4);
% SlopeHip3=(.2)/(10e-5);
% intKnee1=0;
% intKnee2=0;
% intKnee3=0;
% intHip1=0;
% intHip2=0;
% intHip3=0;

% SlopeKnee1=(.2)/(1e-4);
% SlopeKnee2=(.2)/(2.5e-5);
% SlopeKnee3=(.2)/(7.5e-5);
% SlopeHip1=(.2)/(6e-5);
% SlopeHip2=(.2)/(1.4e-4);
% SlopeHip3=(.2)/(4e-5);
% intKnee1=0;
% intKnee2=0;
% intKnee3=0;
% intHip1=0;
% intHip2=0;
% intHip3=0;

ChangeSlope1=1;
ChangeSlope2=1;
ChangeSlope3=1;
ChangeSlope4=1;
ChangeSlope5=1;
ChangeSlope6=1;


ChangeSlope=[ChangeSlope1; ChangeSlope2; ChangeSlope3; ChangeSlope4; ChangeSlope5; ChangeSlope6];
function main()


%% The beginning
clear all;

addpath('supportingFunctions')

global Iext_CPG Iext_MN;
%global xffront yffront xfmid yfmid xfrear yfrear comXimp comYimp; 
global  comXimp comYimp;
global ximp yimp
global Vdes omega;
% global input timevector workrearhip univcount turn; 
global vmax moment_arm tors_k cdamp highvmax;
% global randspike arrowx arrowy plotfoot leftthetainit leftthetadotinit;
% global expno gaussfoot param pert;
global Fxfrontarray Fyfrontarray Fxmidarray Fymidarray Fxreararray Fyreararray tarray activarray3 activarray4 activarray11 activarray12 activarray19 activarray20 tarray2;
global no mm ttau ttau2 teta teta2 Theta_MNR ttau3 ttau4 teta3 teta4
global kneeTorquearray hipTorquearray
global Deltaarray1 Deltaarray2 nvectorarray activarray

global xffrontleft yffrontleft xfmidright yfmidright xfrearleft yfrearleft;
global xffrontright yffrontright xfmidleft yfmidleft xfrearright yfrearright;
global Ffrontleftvector Fmidrightvector Frearleftvector Ffrontrightvector Fmidleftvector Frearrightvector
global Fxfrontleftvector Fyfrontleftvector Fxmidrightvector Fymidrightvector Fxrearleftvector Fyrearleftvector
global thetasave Raghu_MNPhases gsynfbet gsynfbit savevector

%The top for-loop analyzes the impact of different synaptic strengths given
%the lateral perturbations. You can change how this $strengths$ vector
%varies 

%strengths=[0 0.05 .1 .15 .2 .25 .3 .35 .4 .45 .5];
strengths = 0.01;
gsynfbevector=[ones(24,1)*strengths];
gsynfbivector=[zeros(24,length(strengths))]; 

for Pertstep=1:length(strengths)
    
    savevector=[];
    
footdatamusc; 
Phaseneuraldata3;                 %% Parameters for the entire system loaded in neuraldata
datamusc;
%Series of vectors created to keep track of specific variables F=Force;
%x,y=direction ; front,mid,rear=which leg

tarray=[];

Fxfrontvector = []; Fyfrontvector = [];
Fxmidvector = []; Fymidvector = [];
Fxrearvector = []; Fyrearvector = [];
activarray3=[];  activarray4=[]; tarray2=[];
kneeTorquearray=[]; hipTorquearray=[];
activarray=[];
thetasave=[];
inputvector=[];
inputMNvector=[];

% c='kbcgmry'; 
vmax=75; moment_arm=0.4; tors_k=2e-5; cdamp=4e-6; highvmax=75;      
Vdes=0.26;
if Vdes<0.375
    omega=2*pi*(13+(Vdes-0.375)/0.0375);
else
    omega=2*pi*13;
end
% timeperiod = 2*pi/omega*1000;% In milliseconds

%% Some CPG and MN parameters

Iext_CPG=36.39; % 35.6, 35.8  % External input current to CPG
Iext_MN = [36.48 36.24 37.9 33.87   36.48 33.37 39.3 36.33   38.3 36.24 37.32 36.42 ...
        36.48 36.24 37.9 33.87   36.48 33.37 39.3 36.33    38.3 36.24 37.32 36.42].';  %%External input current to Motoneurons
mns=24;     %%Number of Motoneurons


%% Calculating initial theta
options=optimset('TolX',1e-10);
%[fixpoint,funcval]=fsolve(@evathetamusc,[0 0],options,0,0);
fixpoint =[0.0744 0.8505];
leftthetainit=fixpoint(1,1);

%evathetamusc(fixpoint,0,1);  % The fit function has been changed by Matlab


%% Initialize vectors

no=ones(24,1);
mm=ones(24,1);
ttau=1000*ones(24,1);
ttau2=1000*ones(24,1); 
ttau3=1000*ones(24,1);
ttau4=1000*ones(24,1);
teta=1000*ones(24,1); 
teta2=1000*ones(24,1);
teta3=1000*ones(24,1); 
teta4=1000*ones(24,1);
Theta_MNR=ones(24,1);

% univcount=1; turn=0; gaussfoot=0; 
n=0; 
% seemovie=1; trailer=1; normal=0; report=1; 
% expno=0; randspike=0; workloop=0;


dV=0; ddelta=0; dthetadot=0;
input(1,1)=0.2526 + dV; input(1,2)=-0.0217 + ddelta;
input(1,3)=leftthetainit; input(1,4)=-1.8643 + dthetadot;

comXimp=0;   comYimp=0;
[xffront,yffront,xfmid,yfmid,xfrear,yfrear]=getfootmusc(input,comXimp,comYimp,n);


comXvector=[]; comXdotvector=[]; comYvector=[]; comYdotvector=[];
thetavector=[]; thetadotvector=[]; timevector=[]; 
xfootfrontvector=[]; yfootfrontvector=[]; 
xfootmidvector=[]; yfootmidvector=[]; 
xfootrearvector=[]; yfootrearvector=[];
nvector=[]; stancecount=[];
Deltaarray1=[];  Deltaarray2=[]; nvectorarray=[];



ttausave=[];
ttau2save=[];
ttau3save=[];
ttau4save=[];
tetasave=[];
teta2save=[];
teta3save=[];
teta4save=[];


CPGvector=[];  MNvector=[];  Sensvector=[];      %%%%%NEW%%%%%%%

Ffrontleftvector=[]; 
Fmidrightvector=[]; 
Frearleftvector=[]; 
Ffrontrightvector=[]; 
Fmidleftvector=[]; 
Frearrightvector=[];

Fxfrontleftvector=[]; 
Fyfrontleftvector=[]; 
Fxmidrightvector=[]; 
Fymidrightvector=[]; 
Fxrearleftvector=[]; 
Fyrearleftvector=[];

xinit=[3.4627 6.1645 2.5837 0.0410 2.2502 5.8529 2.2502 4.4580 3.8858 4.6312 ...
    5.8557 5.9589 2.0828 3.5633 5.8557 3.7711 5.9505 3.3748 5.9505 1.3016 0.7889 1.4794 2.5836 ...
    3.4184 5.0950];

input=[0.3162    0.2366    0.7280   -3.5411];
SensInit=[];
hh=1;
%% The stance to stance loop

tendlast=0;
firsttime=1;

 LO=.4489;
% 
% 
ExtPhase=mod((1-[.9988;.9531;.0413;.2896;.2923;0]),1);
FlexPhase=mod((1-[.5720;.4810;.4692;.3088;.9887;.4600]),1);
%xinit(1)=LLO;
xinit(2)=mod(ExtPhase(4),1);
xinit(3)=mod(FlexPhase(4),1);
xinit(4)=mod(ExtPhase(1),1);
xinit(5)=mod(FlexPhase(1),1);
xinit(6)=mod(ExtPhase(5),1);
xinit(7)=mod(FlexPhase(5),1);
xinit(8)=mod(ExtPhase(2),1);
xinit(9)=mod(FlexPhase(2),1);
xinit(10)=mod(ExtPhase(6),1);
xinit(11)=mod(FlexPhase(6),1);
xinit(12)=mod(ExtPhase(3),1);
xinit(13)=mod(FlexPhase(3),1);


xinit(14)=mod(xinit(2)+.5,1);
xinit(15)=mod(xinit(3)+.5,1);
xinit(16)=mod(xinit(4)+.5,1);
xinit(17)=mod(xinit(5)+.5,1);
xinit(18)=mod(xinit(6)+.5,1);
xinit(19)=mod(xinit(7)+.5,1);
xinit(20)=mod(xinit(8)+.5,1);
xinit(21)=mod(xinit(9)+.5,1);
xinit(22)=mod(xinit(10)+.5,1);
xinit(23)=mod(xinit(11)+.5,1);
xinit(24)=mod(xinit(12)+.5,1);
xinit(25)=mod(xinit(13)+.5,1);


xinit(1)=.56;
%xinit(26)=mod(xinit(1)+.5,1);
xinit(26)=.56;

ttausave=[];
ttau2save=[];
ttau3save=[];
ttau4save=[];
tetasave=[];
teta2save=[];
teta3save=[];
teta4save=[];

% L stance
xffrontleft=0; yffrontleft=0; 
xfmidright=0; yfmidright=0; 
xfrearleft=0; yfrearleft=0;
% R stance
xffrontright=0; yffrontright=0; 
xfmidleft=0; yfmidleft=0; 
xfrearright=0; yfrearright=0;

xfootvector=[]; yfootvector=[];
TDindexvector=[];
contactvector=[];

TDindex=[1 1 1 0 0 0];

ximp=0;   yimp=0;
state6 = fs2spfile(input,n);
  
    gsynfbet=gsynfbevector(:,Pertstep);
    gsynfbit=gsynfbivector(:,Pertstep);
 
 for j=1:50
    clear t; clear x; clear te xe ie;
    t=[];  x=[];
    options=odeset('Events',@evtNconst2,'RelTol',1e-6,'AbsTol',1e-7);
    ic=fs2spmusc(input,n);
    tspan=[0 52];
    
    x0 = [xinit ic.*[1 1e-3 1 1e-3 1 1e-3] SensInit'].';   %%%%%NEW%%%%%%%
    
    
     for Legindex=1:6
           if TDindex(Legindex)==1
            getfoot(state6,ximp,yimp,Legindex); 
           end
     end
    
    for jj=1:40
    
        
    [t1,x1,te,xe,ie] = ode45(@eqns2,tspan,x0,options,n,tendlast,TDindex);
     
    state6=x1(end,27:32);
    ximp=x1(end,27);
    yimp=x1(end,29);
    
    lengthie=length(ie);
    %ie
        t=[t; t1];
        x=[x; x1];
        x0 = [x1(end,:) SensInit'].';   %%%%%NEW%%%%%%%
        tspan=[t1(end) 52];
        
        
ttausave=[ttausave ttau*ones(1,length(t1))];
ttau2save=[ttau2save ttau2*ones(1,length(t1))];
ttau3save=[ttau3save ttau3*ones(1,length(t1))];
ttau4save=[ttau4save ttau4*ones(1,length(t1))];
tetasave=[tetasave teta*ones(1,length(t1))];
teta2save=[teta2save teta2*ones(1,length(t1))];
teta3save=[teta3save teta3*ones(1,length(t1))];
teta4save=[teta4save teta4*ones(1,length(t1))]; 

     xfootvector=[xfootvector; ones(size(t1))*[xffrontleft xfmidright xfrearleft xffrontright xfmidleft xfrearright]]; %#ok<AGROW>
     yfootvector=[yfootvector; ones(size(t1))*[yffrontleft yfmidright yfrearleft yffrontright yfmidleft yfrearright]]; %#ok<AGROW>
     TDindexvector=[TDindexvector; ones(size(t1))*TDindex];
        
    for kk=1:lengthie
     if ie(kk)==3
       x0(1)=0;
     end
    end
    
        for kk=1:lengthie
     if ie(kk)==28
        x0(26)=0;
     end
    end
    
    for kk=1:lengthie
        for jjj=1:24
        if ie(kk)==jjj+3
        if n==0
            tf=t1(end);
        else
        tf=timevector(end)+t1(end);
        end
        ttau(jjj)=tf;
        x1(end,jjj+1);
        
        teta(jjj)=tf+.5;
    	no(jjj)=no(jjj)+1;
        ttau2(jjj)=ttau(jjj)+ONSET(jjj);
        teta2(jjj)=ttau(jjj)+ONSET(jjj)+.5;
        x0(jjj+1)=0;
        
        if (jjj==1)||(jjj==2)||(jjj==9)||(jjj==10)||(jjj==13)||(jjj==14)||(jjj==21)||(jjj==22)  %%||(jjj==3)||(jjj==4)||(jjj==15)||(jjj==16)
        ttau3(jjj)=ttau2(jjj)+ONSET(jjj);
        teta3(jjj)=ttau2(jjj)+ONSET(jjj)+.5;
        end
        
        if (jjj==1)||(jjj==2)||(jjj==13)||(jjj==14)
        ttau4(jjj)=ttau3(jjj)+ONSET(jjj);
        teta4(jjj)=ttau3(jjj)+ONSET(jjj)+.5;
        end
        
ttausave(:,end)=ttau;
ttau2save(:,end)=ttau2;
ttau3save(:,end)=ttau3;
ttau4save(:,end)=ttau4;
tetasave(:,end)= teta;
teta2save(:,end)=teta2;
teta3save(:,end)=teta3;
teta4save(:,end)=teta4;

         end

        end
    end
    
    %% This clause changes the gait for the set duty cycle.

    for kk=1:lengthie
        if (ie(kk)==2)||(ie(kk)==1) 
    if (ie(kk)==2)
        x(end,26)=-.95;
    elseif (ie(kk)==1)
        x(end,1)=-.45;
    end
    
    hh=hh+1;
    
    break
        end
    end
    
    if ((ie(kk)==2)||(ie(kk)==1))&&(hh==2)
        TDindex=1-TDindex;     
        hh=1;
        
        break
    end
   
    
    end
    
    comXimp=x(end,27);
    comYimp=x(end,29);
    
    ximp=x(end,27);
    yimp=x(end,29);
    trans = x(end,27:32).*[1 1000 1 1000 1 1000];
    output=ep2fsmusc(trans,n);
    state6=x1(end,27:32);
    
    %xinit = x(end,1:240);
    
    xinit=x(end,1:26);

    %disp('is it the first time?')
    if firsttime==1
        firsttime=2;
        timevector=[];
    else
        t=t+timevector(length(timevector));
       
%     Here the TD/LO time instant is repeated once, at each TD/LO transition
    end
    
    comXvector=[comXvector; x(:,27)]; %#ok<AGROW>
    comXdotvector=[comXdotvector; x(:,28).*1000]; %#ok<AGROW>
    comYvector=[comYvector; x(:,29)]; %#ok<AGROW>
    comYdotvector=[comYdotvector; x(:,30).*1000]; %#ok<AGROW>
    thetavector=[thetavector; x(:,31)]; %#ok<AGROW>
    thetadotvector=[thetadotvector; x(:,32).*1000]; %#ok<AGROW>
    
    % foot position and others
%     xfootfrontvector=[xfootfrontvector; xffront*ones(size(t))]; %#ok<AGROW>
%     yfootfrontvector=[yfootfrontvector; yffront*ones(size(t))]; %#ok<AGROW>
%     xfootmidvector=[xfootmidvector; xfmid*ones(size(t))]; %#ok<AGROW>
%     yfootmidvector=[yfootmidvector; yfmid*ones(size(t))]; %#ok<AGROW>
%     xfootrearvector=[xfootrearvector; xfrear*ones(size(t))]; %#ok<AGROW>
%     yfootrearvector=[yfootrearvector; yfrear*ones(size(t))]; %#ok<AGROW>
%     
    
    nvector=[nvector; ((-1)^n)*ones(size(t))]; %#ok<AGROW>
    stancecount=[stancecount; n*ones(size(t))]; %#ok<AGROW>
    
    CPGvector=[CPGvector; x(:,1)];       %%%%%NEW%%%%%%%
    MNvector=[MNvector; x(:,2:25)];       %%%%%NEW%%%%%%%
    %Sensvector=[Sensvector; x(:,247:258)]; %%%%%NEW%%%%%%%
    contactvector=[contactvector; ones(size(t))*TDindex];
    
    
    tendlast = t(end,1);
    
    timevector=[timevector; t]; %#ok<AGROW>
      
    disp(n);
    
    n=n+1;
    input=output;
    inputvector=[inputvector; input];
    if mod(n,2)==0
    inputMNvector=[inputMNvector; MNvector(end,:)];
    end
    
    [xffront,yffront,xfmid,yfmid,xfrear,yfrear]=getfootmusc(input,comXimp,comYimp,n);
       
 end
 
 
end

% This is where you would save the output from a simulation.  The output
% can then be used to parse the 

figure
plot(comXvector, comYvector,'k','LineWidth',2)
axis square
set(gca,'FontSize',22,'LineWidth',2)
xlabel('X')
ylabel('Y')
title('COM trajectory with perturbation')

save TESTCase;

end



%end

%% The dynamic equations
function [sys] = eqns2(t,x,n,tendlast,TDindex)

global Iext_CPG Iext_MN;
global Fxfrontarray Fyfrontarray Fxmidarray Fymidarray Fxreararray Fyreararray tarray activarray3 activarray4 activarray11 activarray12 activarray19 activarray20 tarray2;
global no mm ttau ttau2 teta teta2 Theta_MNR ttau3 ttau4 teta3 teta4
global kneeTorquearray hipTorquearray
global Deltaarray1 Deltaarray2 nvectorarray activarray
global Ffrontleft Fmidright Frearleft Ffrontright Fmidleft Frearright
global Ffrontleftvector Fmidrightvector Frearleftvector Ffrontrightvector Fmidleftvector Frearrightvector
global Fxfrontleftvector Fyfrontleftvector Fxmidrightvector Fymidrightvector Fxrearleftvector Fyrearleftvector
global thetasave gsynfbet gsynfbit savevector

Phaseneuraldata3; datamusc;

% Time is in ms, voltage is in mV, conductances are in mS/cm^2, slope
% coefficients are in mV/s, capacitance is in microF/cm^2
  
tarray = [tarray; t + tendlast];    %%Keep track of inner time

Theta_CPG=mod(x(1),1);
Theta_MN=mod(x(2:25,1),1);
Theta_CPG2=mod(x(26),1);
x(1,1)=mod(x(1,1),1);

comX = x(27,1);            %% Center of Mass trajectory in X
comXdot = x(28,1);         %% Center of Mass velocity in X
comY = x(29,1);            %% Center of Mass trajectory in Y
comYdot = x(30,1);         %% Center of Mass velocity in Y
theta = x(31,1);           %% Body orientation theta
thetadot = x(32,1);        %% Angular velocity of the Body 


six = x(27:32,1).*[1; 1000; 1; 1000; 1; 1000];

thetasave=[thetasave; theta];


Thetadot_CPG=w0_CPG/2/pi;
Thetadot_CPG2=w0_CPG/2/pi;

smalleps=.1;
Theta_CPGf=mod([ones(12,1)*Theta_CPG;ones(12,1)*(Theta_CPG+.5)],1);
phasediffcalc=mod(Theta_MN-Theta_CPGf,1);

%Index_MN
for jj=1:24

H_fu(jj,1)=interp1q(newvec,DiffH_Func(:,jj),phasediffcalc(jj));
EPRCe(jj,1)=interp1q(newvec,DiffEPRCe_MN(:,jj),Theta_MN(jj));
EPRCi(jj,1)=interp1q(newvec,DiffEPRCi_MN(:,jj),Theta_MN(jj));

end
gg=.1;

%Symmetric or Unsymmetric feedback values
%gsyn=[gg*3; gg*2; gg; gg; gg*3; .1; gg; gg; gg; gg*2; gg; gg; 
%      gg*3; gg*2; gg; gg; gg*3; .1; gg; gg; gg; gg*2; gg; gg];

gsyn=[gg; gg; gg; gg; gg; gg; gg; gg; gg; gg; gg; gg; 
      gg; gg; gg; gg; gg; gg; gg; gg; gg; gg; gg; gg];


%% 24 muscles - excitation dynamics

a1=6*2;                           %% Parameters for the muscle model
a2=0.1804*4;                      %% a1-a3 are for Beta
a3=0.8923*3.2; 

a4=19*2;                          %%Parameters for the muscle model
a5=0.6023*4;                      %% a4-a6 are for Gammeta
a6=5.2469*3.2; 

rho=1.272;                        %% rho and a0 are parameters for the activation equation
a0=0;

r(1)=0;
r(2)=( -a1 + sqrt(a1^2-4*a2) )/2;
r(3)=( -a1 - sqrt(a1^2-4*a2) )/2;
r(4)=( -a4 + sqrt(a4^2-4*a5) )/2;
r(5)=( -a4 - sqrt(a4^2-4*a5) )/2;
A=zeros(1,5);

for k=1:5
    A(k)=inv( ( r(k)-r(mod(k+1,5)+(k+1==5)*5) )*( r(k)-r(mod(k+2,5)+(k+2==5)*5) )*( r(k)-r(mod(k+3,5)+(k+3==5)*5) )*( r(k)-r(mod(k+4,5)+(k+4==5)*5) ) );
end

gamma=0;

for j=1:24
   
for k=1:5
                first=0;
                second=0;
                if (tarray(end)-ttau(j))>=0 
                first = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-ttau(j)) );
                end
                if (tarray(end)-teta(j))>=0
                second = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-teta(j)) );
                end
                gamma = gamma + first - second;
 end

 
for k=1:5
                first2=0;
                second2=0;
                if (tarray(end)-ttau2(j))>=0
                first2 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-ttau2(j)) );
                end
                if (tarray(end)-teta2(j))>=0
                second2 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-teta2(j)) );
                end
                gamma= gamma + first2 - second2;          
end

if (j==1)||(j==2)||(j==9)||(j==10)||(j==13)||(j==14)||(j==21)||(j==22)||(j==3)||(j==4)||(j==15)||(j==16)
    
for k=1:5
                first3=0;
                second3=0;
                if (tarray(end)-ttau3(j))>=0
                first3 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-ttau3(j)) );
                end
                if (tarray(end)-teta3(j))>=0
                second3 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-teta3(j)) );
                end
                gamma= gamma + first3 - second3;          
end 

end

if (j==1)||(j==2)||(j==13)||(j==14)
    
for k=1:5
                first4=0;
                second4=0;
                if (tarray(end)-ttau4(j))>=0
                first4 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-ttau4(j)) );
                end
                if (tarray(end)-teta4(j))>=0
                second4 = (a6*a3*A(k)).*exp( r(k).*(tarray(end)-teta4(j)) );
                end
                gamma= gamma + first4 - second4;          
end 

end


exteact(j)=(a0+(rho.*gamma).^2)./(1+(rho.*gamma).^2);
gamma=0;

end

Theta_MNR=Theta_MN;


%% Activation function

activ= [exteact' ];


%% CoM dynamics

tins = t(1,1)*0.001; % t in ms taken care of 

actibreak2 = activ(1:12,1)*(mod(n,2)==0) + activ(13:24,1)*(mod(n,2)==1);
actibreak=activ;
% 
% activarray=[activarray actibreak];

outarray=getvarmusc(tins,six,n,actibreak2); % Have to define n

Fxfront=outarray(1,1);  %Fxfrontarray = [Fxfrontarray; Fxfront];
Fyfront=outarray(2,1); %Fyfrontarray = [Fyfrontarray; Fyfront];
Fxmid=outarray(3,1);   %Fxmidarray = [Fxmidarray; Fxmid];
Fymid=outarray(4,1);    %Fymidarray = [Fymidarray; Fymid];
Fxrear=outarray(5,1);   %Fxreararray = [Fxreararray; Fxrear];
Fyrear=outarray(6,1);   %Fyreararray = [Fyreararray; Fyrear];
tau1front=outarray(7,1);
tau1mid=outarray(9,1);
tau1rear=outarray(11,1);

% Left stance
outarray=frontleft1(tins,six,n,TDindex(1,1),actibreak(1:4));
Fxfrontleft=outarray(1,1);
Fyfrontleft=outarray(2,1);
tau1frontleft=outarray(3,1);
tau2frontleft=outarray(4,1);
Ffrontleft=sqrt(Fxfrontleft^2+Fyfrontleft^2);

outarray=midright2(tins,six,n,TDindex(1,2),actibreak(5:8));
Fxmidright=outarray(1,1);
Fymidright=outarray(2,1);
tau1midright=outarray(3,1);
tau2midright=outarray(4,1);
Fmidright=sqrt(Fxmidright^2+Fymidright^2);

outarray=hindleft3(tins,six,n,TDindex(1,3),actibreak(9:12));
Fxrearleft=outarray(1,1);
Fyrearleft=outarray(2,1);
tau1rearleft=outarray(3,1);
tau2rearleft=outarray(4,1);
Frearleft=sqrt(Fxrearleft^2+Fyrearleft^2);


% Right stance

outarray=frontright4(tins,six,n,TDindex(1,4),actibreak(13:16)); 
Fxfrontright=outarray(1,1);
Fyfrontright=outarray(2,1);
tau1frontright=outarray(3,1);
tau2frontright=outarray(4,1);
Ffrontright=sqrt(Fxfrontright^2+Fyfrontright^2);

outarray=midleft5(tins,six,n,TDindex(1,5),actibreak(17:20));
Fxmidleft=outarray(1,1);
Fymidleft=outarray(2,1);
tau1midleft=outarray(3,1);
tau2midleft=outarray(4,1);
Fmidleft=sqrt(Fxmidleft^2+Fymidleft^2);

outarray=hindright6(tins,six,n,TDindex(1,6),actibreak(21:24));
Fxrearright=outarray(1,1);
Fyrearright=outarray(2,1);
tau1rearright=outarray(3,1);
tau2rearright=outarray(4,1);
Frearright=sqrt(Fxrearright^2+Fyrearright^2);

 
%size(outarray)

%  Need to Establish Sensory Neurons Governing Equations
%  And Add the appropriate external current depending on the Torque at each
%  of the Six knees. 

kneeTorque=[tau2frontleft; tau2midright; tau2rearleft; tau2frontright; tau2midleft; tau2rearright];
hipTorque=[tau1frontleft; tau1midright; tau1rearleft; tau1frontright; tau1midleft; tau1rearright];

[Delta_FBi,Delta_FBe]= GetDelta(kneeTorque,hipTorque,TDindex);


kk=mod(n,2);   %%Evaluates whether we are on the left or right Tripod

Thetadot_MN=w0_CPG/2/pi+smalleps.*(Delta+H_fu-gsynfbet.*Delta_FBe.*EPRCe-gsynfbit.*Delta_FBi.*EPRCi);

% Defining the lateral Perturbation experiment given by extforce started 
% at stance 10. 

extforce=0;

if (n==20)
    if (tins>=0.006)&&(tins<=0.009)
        extforce=0.42*(tins-0.006)/0.003;
    elseif (tins>=0.009)&&(tins<=0.01)
        extforce=0.42-0.42*(tins-0.009)/0.001;
    elseif (tins<0.006)
        arrowx=six(1);
        arrowy=six(3);
    end
end

% Dynamic Equations for the cockroach body moving in a plane

n=0;
leftstancemoment = (((-1)^n)*(-tau1frontleft+tau1midright-tau1rearleft)+ ...
    ((-1)^(n+1))*Fyfrontleft*(d1front*cos(theta)+((-1)^n)*d2front*sin(theta)) - ...
    Fxfrontleft*(d2front*cos(theta)-((-1)^n)*d1front*sin(theta)) + ...
    ((-1)^n)*Fymidright*(d1mid*cos(theta)+((-1)^(n+1))*d2mid*sin(theta)) - ...
    Fxmidright*(d2mid*cos(theta)-((-1)^(n+1))*d1mid*sin(theta)) + ...
    ((-1)^(n+1))*Fyrearleft*(d1rear*cos(theta)+((-1)^n)*d2rear*sin(theta)) - ...
    Fxrearleft*(d2rear*cos(theta)-((-1)^n)*d1rear*sin(theta)))/I;


n=1;
rightstancemoment = (((-1)^n)*(-tau1frontright+tau1midleft-tau1rearright)+ ...
    ((-1)^(n+1))*Fyfrontright*(d1front*cos(theta)+((-1)^n)*d2front*sin(theta)) - ...
    Fxfrontright*(d2front*cos(theta)-((-1)^n)*d1front*sin(theta)) + ...
    ((-1)^n)*Fymidleft*(d1mid*cos(theta)+((-1)^(n+1))*d2mid*sin(theta)) - ...
    Fxmidleft*(d2mid*cos(theta)-((-1)^(n+1))*d1mid*sin(theta)) + ...
    ((-1)^(n+1))*Fyrearright*(d1rear*cos(theta)+((-1)^n)*d2rear*sin(theta)) - ...
    Fxrearright*(d2rear*cos(theta)-((-1)^n)*d1rear*sin(theta)))/I;

%  Here are the governing equations for the mechanical system.  This is
%  where hill climbing or descending can directly implemented.

comXdoubdot=1e-6*(Fxfrontleft+Fxmidright+Fxrearleft+Fxfrontright+Fxmidleft+ ...
     Fxrearright-extforce*cos(six(5)))/m;
comYdoubdot=1e-6*(Fyfrontleft+Fymidright+Fyrearleft+Fyfrontright+Fymidleft+ ...
     Fyrearright-extforce*sin(six(5)))/m;
thetadoubdot=1e-6*(leftstancemoment + rightstancemoment);
 
 
%% Derivatives

sys = [Thetadot_CPG; Thetadot_MN; Thetadot_CPG2;
    comXdot; comXdoubdot; comYdot; comYdoubdot; thetadot; thetadoubdot];

if nargout>1
    FxfrontT = Fxfrontarray;
    FyfrontT = Fyfrontarray;
    FxmidT = Fxmidarray;
    FymidT = Fymidarray;
    FxrearT = Fxreararray;
    FyrearT = Fyreararray;
    tarrayT = tarray;
end

end

%% The infinity function
function inftemp = inffunc(v,k0,vTH)
inftemp = 1./( 1 + exp(-2.*k0.*(v-vTH)) );
end


%% The tau function
function tau = taufunc(v,k0,vTH)
tau = sech( k0.*(v-vTH)./2 );
end

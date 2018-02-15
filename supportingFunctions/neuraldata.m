
% All these parameters are assumed to be of the same value for all the CPG
% neurons. CPGs %%%%

cpgs = 6;                       %% Number of CPGs
vec_CPG  = ones(cpgs,1);        %% Vector of CPGs:  Refer to NamingDiagram.jpg for details of CPG numbering

k0n_CPG = 1/18.*vec_CPG;        %% parameters k for the different gating variable steady state equations w_infinity
k0m_CPG = 0.1.*vec_CPG;         %%
k0c_CPG = 0.8.*vec_CPG;         %%
vTHn_CPG = -1.2.*vec_CPG;       %% parameters for reversal potentials of time constants tau_i for different variables
vTHm_CPG = 2.*vec_CPG;          %%
vTHc_CPG = -27.*vec_CPG;        %%

gCa_CPG = 4.4.*vec_CPG;         %% Calcium conductance for CPG
gK_CPG  = 9.*vec_CPG;           %% Potassium conductance for CPG
gL_CPG  = 2.*vec_CPG;           %% Leakage conductance for CPG
gKS_CPG = [0.126; 0.126.*ones(5,1)];%0.19.*vec_CPG; %% Slow Potassium conductance for CPG

C_CPG = 1.2.*vec_CPG;           %% Capacitance of CPG         

ECa_CPG = 120.*vec_CPG;         %% Calcium Reversal Potential
EK_CPG  = -80.*vec_CPG;         %% Potassium Reversal Potential 
EL_CPG  = -60.*vec_CPG;         %% Leakage Current Reversal Potential
EKS_CPG = EK_CPG; % where is this %% Slow Potassium Reversal Potential

epsil_CPG = 4.9.*vec_CPG;       %% Epsilon for the dynamics of d/dt (m) of CPG
delt_CPG  = [0.0522; 0.0522.*ones(5,1)];%0.052.*vec_CPG;    %%Delta for the dynamics of d/dt (c)

alph_CPG = 5000.*vec_CPG;       %% alpha found within synaptic dynamics of CPG
bet_CPG  = 0.18.*vec_CPG;       %% beta found within synaptic dynamics of CPG
Tmax_CPG = 2e-3.*vec_CPG;       %% Parameter T_max found within transmitter concentration G Equation
kpre_CPG = 0.22.*vec_CPG;       %% Parameter k_pre found within transmitter concentration G Equation
Esynpre_CPG  = 2.*vec_CPG;      %% Parameter E_synpre within transmitter concentration G Equation
Esynpost_CPG = -70.*vec_CPG;    %% Reversal Potential for the synaptic input current in the postsynaptic neuron
gsyn_CPG = 0.01.*vec_CPG; % Technically, there are 14 synapses. %% Conductance of synaptic current in post neuron


% MNs %%%%%

mns = 24;                           %% Number of Motoneurons
vec_MN = ones(mns,1);               %% Vector of Motoneurons refer to the diagram to see numbering formula
% inhibvec_MN = (-1).^(1:1:mns).';
inhibvec_MN = (-1).*ones(mns,1);

k0n_MN = 0.056.*vec_MN;         %% parameters k for the different gating variable steady state equations w_infinity
k0m_MN = 0.1.*vec_MN;           %%
k0c_MN = 0.8.*vec_MN;           %%
vTHn_MN = -1.2.*vec_MN;         %% parameters for reversal potentials of time constants tau_i for different variables
vTHm_MN = 2.*vec_MN;            %%
vTHc_MN = -27.*vec_MN;          %%

gCa_MN = 4.4.*vec_MN;           %% Calcium conductance for Motoneuron
gK_MN  = 9.*vec_MN;             %% Potassium conductance for for Motoneuron
gL_MN  = 2.*vec_MN;             %% Leakage conductance for Motoneuron
gKS_MN = [0.15;  0.25; 0.4; 0.19;   %% Slow Potassium conductance for Motoneurons
        0.15; 0.19; 0.77; 0.7;
        0.33; 0.25; 0.25; 0.73;
        0.15; 0.25; 0.4; 0.19;
        0.15; 0.19; 0.77; 0.7;
        0.33; 0.25; 0.25; 0.73];%0.25.*vec_MN;

C_MN = [1.35; 1.39; 1.35; 1.39;     %% Capacitances of CPGs
        1.35; 1.39; 1.35; 1.28;
        1.35; 1.39; 1.35; 1.28;
        1.35; 1.39; 1.35; 1.39;
        1.35; 1.39; 1.35; 1.28;
        1.35; 1.39; 1.35; 1.28];%1.39.*vec_MN;

ECa_MN = 120.*vec_MN;               %% Calcium Reversal Potential
EK_MN  = -80.*vec_MN;               %% Potassium Reversal Potential 
EL_MN  = -60.*vec_MN;               %% Leakage Current Reversal Potential
EKS_MN = EK_MN;                     %% Slow Potassium Reversal Potential

epsil_MN = [3.59; 1.3; 3.84; 3.78;      %% Epsilon for the dynamics of d/dt (m) of Motoneurons
        3.577; 3.9; 3.9; 3.85;
        3.83; 1.3; 3.84; 3.85;
        3.59; 1.3; 3.84; 3.78;
        3.577; 3.9; 3.9; 3.85;
        3.83; 1.3; 3.84; 3.85];%4.18.*vec_MN;
delt_MN = [0.011; 0.044; 0.044; 0.034;  %%Delta for the dynamics of d/dt (c) of Motoneurons
        0.011; 0.032; 0.032; 0.016;
        0.025; 0.044; 0.03; 0.016;
        0.011; 0.044; 0.044; 0.034;
        0.011; 0.032; 0.032; 0.016;
        0.025; 0.044; 0.03; 0.016];%0.044.*vec_MN;

alph_MN = [5000; 5000; 5000.*ones(22,1)];%5000.*vec_MN;
bet_MN  = 0.18.*vec_MN;
Tmax_MN = 2e-3.*vec_MN;
kpre_MN = 0.22.*vec_MN;
Esynpre_MN  = 2.*vec_MN;

Esynpost_MN = zeros(mns,1) + -70.*vec_MN.*( mod(vec2ind(vec_MN).',2)==1 );
gsyn_MN = 0.2.*vec_MN; % Technically there are 24 synapses.
% Are these incoming or outgoing synaptic strengths. Ans: For incoming
% right now. gsyn_MN is right now used for the synapses from CPG to MN.


% Separate equations for CPG to MN synaptic dynamics %%%%%

alph_CPG2MN = [1000; 1000; 1000; 1000;      %% alpha found within synaptic dynamics of Motoneurons
            1000; 1000; 1000; 1000;
            1000; 1000; 1000; 1000;
            1000; 1000; 1000; 1000;
            1000; 1000; 1000; 1000;
            1000; 1000; 1000; 1000] + ...
    4000.*vec_MN.*( mod(vec2ind(vec_MN).',2)==1 );
bet_CPG2MN  = [0.18; 0.18; 0.18; 0.18;      %% beta found within synaptic dynamics of Motoneurons
            0.18; 0.18; 0.18; 0.18;
            0.18; 0.18; 0.18; 0.18;
            0.18; 0.18; 0.18; 0.18;
            0.18; 0.18; 0.18; 0.18;
            0.18; 0.18; 0.18; 0.18] + ...
    0.01.*vec_MN.*( mod(vec2ind(vec_MN).',2)==0 );
Tmax_CPG2MN = 2e-3.*vec_MN;                 %% Parameter T_max found within transmitter concentration G Equation
kpre_CPG2MN = 0.22.*vec_MN;                 %% Parameter k_pre found within transmitter concentration G Equation
Esynpre_CPG2MN  = 2.*vec_MN;                %% Parameter E_synpre within transmitter concentration G Equation



% Parameters for the excitation dynamics equations of the muscles

muscles = 24;                               %% Number of Muscles
vec_musc = ones(muscles,1);                 %% Vector of muscles

a1=6*2.*vec_musc;                           %% Parameters for the muscle model
a2=0.1804*4.*vec_musc;                      %% a1-a3 are for Beta
a3=0.8923*3.2.*vec_musc; 

a4=19*2.*vec_musc;                          %%Parameters for the muscle model
a5=0.6023*4.*vec_musc;                      %% a4-a6 are for Gammeta
a6=5.2469*3.2.*vec_musc; 

rho=1.272.*vec_musc;                        %% rho and a0 are parameters for the activation equation
a0=0.*vec_musc;


% Muscle scaling for the various muscles.

muscS = [2.1959; 1.1349; 2.1397; 0.5109;            
        5.7958; 0.2369; 1.9435; 0.5565;
        1.8442; 2.0354; 2.2027; 1.0849;
        2.1959; 1.1349; 2.1397; 0.5109;
        5.7958; 0.2369; 1.9435; 0.5565;
        1.8442; 2.0354; 2.2027; 1.0849];
    
    
% Sensory Neuron Parameters    

SensoryNeurons=12;                          %%Number of Sensory Neurons
vec_SN=ones(SensoryNeurons,1);
vec2_SN=ones(SensoryNeurons*2,1);

ENa_Sens=55*vec_SN;                         %%Sodium Reversal Potential
EK_Sens=-72*vec_SN;                         %%Potassium Reversal Potential
EL_Sens=-17*vec_SN;                         %%Leakage Reversal Potential
gNa_Sens=120*vec_SN;                        %%Sodium Conductance
gK_Sens=20*vec_SN;                          %%Potassium Conductance 
gL_Sens=.3*vec_SN;                          %%Leakage Condcutance
gA_Sens=47.7*vec_SN;                        %%Another conductance 
C_Sens=1*vec_SN;                            %%Capacitance of Sensory Neurons
Iib_Sens=5*vec_SN;                          %%External current
gammab_Sens=.069*vec_SN;                    %%Neuron parameter
Tb_Sens=1*vec_SN;                           %%Another Neuron parameter
Tn_Sens=.52*vec_SN;                         %%Neuron parameter
B_Sens=.21*gA_Sens/gK_Sens*vec_SN;          %%Neuron parameter
Iiext_Sens=2.6*vec_SN;                      %%External current

Epres_Sens=2*vec2_SN;                       %%Parameter for presynaptic current equations
Eposts_Sensi=-70;                           %%Inhibitory reversal potential in postsynaptic current
Eposts_Sense=0;                             %%Excitatory reversal potential in postsynaptic current
ksyn_Sens=.22*vec2_SN;                      %%Parameter for presynaptic current equations
gsyn_Sens=.9*vec2_SN;                       %%Conductance of postsynaptic current
Tmax_Sens=2e-3*vec2_SN;                     %%Parameter for presynaptic current equations

QuickRise=5000;     %Typically Inhibitory
SlowRise=1000;      %Typically Excitatory

%Forming a vector for the different alpha, beta values for excitatory or
%inhibitory synapses, using either quickrise or slowrise parameters

for j=1:SensoryNeurons*2    

alph_Sens2MN(j)=QuickRise*(mod(j,2)==1)+QuickRise*(mod(j,2)==0);  %mod(j,2)==1 excitatory, mod(j,2)==0 inhibitory
bet_Sens2MN(j)=.19*(mod(j,2)==1)+.18*(mod(j,2)==0);

end

%These values govern the conductance of the postsynaptic current from the
%sensory neurons to the motoneurons and can be changed to allow more
%current or less or no feedback at all.

gFe=.15; 
gFi=0.15;
gMe=.15; 
gMi=0.15;
gHe=.15; 
gHi=0.15;
% gFe=0; 
% gFi=0;
% gMe=0; 
% gMi=0;
% gHe=0; 
% gHi=0;


%The following are a set of Matrices that are formed to add the correct
%sensory current inputs to the correct Motoneurons as described in the
%paper.

Ze=zeros(2,24);
Nza=[gFe,0,0,gFi,zeros(1,20);0,gFi,gFe,0,zeros(1,20)];
Nzb=[gMe,0,0,gMi,zeros(1,20);0,gMi,gMe,0,zeros(1,20)];
Nzc=[gHe,0,0,gHi,zeros(1,20);0,gHi,gHe,0,zeros(1,20)];

Nz2a=[gFe*Eposts_Sense,0,0,gFi*Eposts_Sensi,zeros(1,20);0,gFi*Eposts_Sensi,gFe*Eposts_Sense,0,zeros(1,20)];
Nz2b=[gMe*Eposts_Sense,0,0,gMi*Eposts_Sensi,zeros(1,20);0,gMi*Eposts_Sensi,gMe*Eposts_Sense,0,zeros(1,20)];
Nz2c=[gHe*Eposts_Sense,0,0,gHi*Eposts_Sensi,zeros(1,20);0,gHi*Eposts_Sensi,gHe*Eposts_Sense,0,zeros(1,20)];


Ma=[Ze;Nza;Ze;circshift(Nzb,[0 4]);Ze;circshift(Nzc,[0 8]);Ze;circshift(Nza,[0 12]);Ze;circshift(Nzb,[0 16]);Ze;circshift(Nzc,[0 20])];
Ma2=[Ze;Nz2a;Ze;circshift(Nz2b,[0 4]);Ze;circshift(Nz2c,[0 8]);Ze;circshift(Nz2a,[0 12]);Ze;circshift(Nz2b,[0 16]);Ze;circshift(Nz2c,[0 20])];

%This is the yint and slope that regulates the external current applied to
%the sensory neurons given a particular torque.

yint=[1 10 1 10 1 10]';
Slope=[35/(10e-5) 70/(10e-5) 35/(1e-4) 70/(10e-5) 35/(10e-5) 70/(1e-4)]';

for j=1:SensoryNeurons
th(j)=(-1)^(j+1);
end

%SensInit=[-30*ones(SensoryNeurons,1);.1*ones(SensoryNeurons,1);.1*ones(SensoryNeurons*2,1)];




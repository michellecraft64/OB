function []=obscOt(TrLabSt,nTrials,flName,lamO,crG,crM,crMP,timeVars,tevok,tspon)
%function to run OBSC (OB Single Compartment) model, 7 cells, implementing
%Ortho stimulus
% INPUTS: TrLabSt=label trial to start, nTrials=# trials per call, 
% flName=file name, with realiz # appended (x1.mat)
% Script to run # of realizations varying lambda_AMPA and lambda_GABA for
% Poisson kick mediation and time-varying crG,crM,crMP
% OUTPUTS: saves spike times in Realz/flName[TrLabSt+jInd-1].mat (spT_[G/M/P])
%set input params
tauO=lamO.tauO;
tauO2=lamO.tauO2;
ssV=lamO.ssV;
lmOevk=lamO.lmOevk;
lmOsp=lamO.lmOsp;
tsft=lamO.tsft;
dt=timeVars.dt;
tSpon=timeVars.tSpon; %ms
tEvok=timeVars.tEvok; %ms
tEnd=tSpon+tEvok;
tmShift=50; %in ms
tevokS=[tevok; tevok(end)+(dt:dt:tmShift)']; %take nuevok out further by tmShift
nmShift=length((dt:dt:tmShift)); %# elem must shift by
nuevokO=lmOevk*(-(tevokS+tsft).*exp(-(tevokS+tsft)/tauO) + (tevokS+tsft).*exp(-(tevokS+tsft)/tauO2))+ssV;
nuevokO=nuevokO(nmShift+1:end); %shift over by tmShift
nu_A=[lmOsp*ones(length(tspon)-1,1);nuevokO];
nu_A=repmat(nu_A,1,7);
nu_G=0.75*nu_A;
spT_G=cell(1,3);
spT_M=cell(1,2);
spT_P=cell(1,2);
%Calculate trial data
tic
for jInd=1:nTrials
    flgGd=0;
    while(flgGd==0) %run until volt isn't complex
        rng('shuffle');
    [flgGd,spTimes_G1,spTimes_G2,spTimes_G3,spTimes_M1,spTimes_M2,spTimes_P1,...
        spTimes_P2]=crr_bnAG_7cells(nu_A,nu_G,crG,crM,crMP,tSpon,tEvok,dt);
    end
    %save results in a cell
    spT_G{1,1}=spTimes_G1;
    spT_G{1,2}=spTimes_G2;
    spT_G{1,3}=spTimes_G3;
    spT_M{1,1}=spTimes_M1;
    spT_M{1,2}=spTimes_M2;
    spT_P{1,1}=spTimes_P1;
    spT_P{1,2}=spTimes_P2;
    
    save([flName,num2str(TrLabSt+jInd-1)],'spT_G','spT_M','spT_P');
end
toc
function [flgGd,spTimes_G1,spTimes_G2,spTimes_G3,spTimes_M1,spTimes_M2,...
        spTimes_P1,spTimes_P2]=crr_bnAG_7cells(nu_A,nu_G,crG,crM,crMP,tSpon,tEvok,dt)
%Function to compute GC, MC, and PGC voltages as well as spike counts and 
%times. Cells are assumed to be coupled by synaptic currents and include
%ORN input. GC, MC, and PGC voltages are approximated using fourth-order 
%Runge Kutta, synapses are approximated used Forward-Euler, and ORN input
%modeled by random Poisson kicks mediated by excitatory and inhibitory
%synapses. Includes transient time iterations, begins counting spikes after.
%Inputs: nu_A, nu_G (Poiss rate of E/I), crG=correl betweeen GC, 
%crM=correl between 2 MCs, crMP=correl between each MC/PGC pair
% Spontaneous state time in ms (tSpon); Evoked state time in ms
%(tEvok); time step (dt).
%
flgGd=1; %start if off as 'good' until complex voltage
spTimes_G1=[]; spTimes_G2=[]; spTimes_G3=[]; spTimes_M1=[]; spTimes_M2=[]; spTimes_P1=[]; spTimes_P2=[];

trn=1000; %transient time of 1s
tEnd=tSpon+tEvok;
Ltrn=round(trn/dt);
Lt=round((tEnd)/dt)+1;
tnu_A=nu_A(1,:); %transient Poisson kick rate
tnu_G=nu_G(1,:); %transient Poisson kick rate
atrn=ones(Ltrn,7);
a=ones(Lt,7); %jump sizes
%Correlate PP
numinA_G12=min(nu_A(:,1),nu_A(:,2)).*crG;
numinA_G23=min(nu_A(:,2),nu_A(:,3)).*crG;
numinA_G13=min(nu_A(:,1),nu_A(:,3)).*crG;
nuaugA_G12=nu_A(:,1)-numinA_G12;nuaugA_G21=nu_A(:,2)-numinA_G12;
nuaugA_G23=nu_A(:,2)-numinA_G23;nuaugA_G32=nu_A(:,3)-numinA_G23;
nuaugA_G13=nu_A(:,1)-numinA_G13;nuaugA_G31=nu_A(:,3)-numinA_G13;
numinA_M=min(nu_A(:,4),nu_A(:,5)).*crM;
nuaugA_M1=nu_A(:,4)-numinA_M;nuaugA_M2=nu_A(:,5)-numinA_M;
numinA_MP11=min(nu_A(:,4),nu_A(:,6)).*crMP;
numinA_MP12=min(nu_A(:,4),nu_A(:,7)).*crMP;
numinA_MP21=min(nu_A(:,5),nu_A(:,6)).*crMP;
numinA_MP22=min(nu_A(:,5),nu_A(:,7)).*crMP;
nuaugA_MP11=nu_A(:,4)-numinA_MP11;nuaugA_PM11=nu_A(:,6)-numinA_MP11;
nuaugA_MP12=nu_A(:,4)-numinA_MP12;nuaugA_PM12=nu_A(:,7)-numinA_MP12;
nuaugA_MP21=nu_A(:,5)-numinA_MP21;nuaugA_PM21=nu_A(:,6)-numinA_MP21;
nuaugA_MP22=nu_A(:,5)-numinA_MP22;nuaugA_PM22=nu_A(:,7)-numinA_MP22;
numinG_G12=min(nu_G(:,1),nu_G(:,2)).*crG;
numinG_G23=min(nu_G(:,2),nu_G(:,3)).*crG;
numinG_G13=min(nu_G(:,1),nu_G(:,3)).*crG;
nuaugG_G12=nu_G(:,1)-numinG_G12;nuaugG_G21=nu_G(:,2)-numinG_G12;
nuaugG_G23=nu_G(:,2)-numinG_G23;nuaugG_G32=nu_G(:,3)-numinG_G23;
nuaugG_G13=nu_G(:,1)-numinG_G13;nuaugG_G31=nu_G(:,3)-numinG_G13;
numinG_M=min(nu_G(:,4),nu_G(:,5)).*crM;
nuaugG_M1=nu_G(:,4)-numinG_M;nuaugG_M2=nu_G(:,5)-numinG_M;
numinG_MP11=min(nu_G(:,4),nu_G(:,6)).*crMP;
numinG_MP12=min(nu_G(:,4),nu_G(:,7)).*crMP;
numinG_MP21=min(nu_G(:,5),nu_G(:,6)).*crMP;
numinG_MP22=min(nu_G(:,5),nu_G(:,7)).*crMP;
nuaugG_MP11=nu_G(:,4)-numinG_MP11;nuaugG_PM11=nu_G(:,6)-numinG_MP11;
nuaugG_MP12=nu_G(:,4)-numinG_MP12;nuaugG_PM12=nu_G(:,7)-numinG_MP12;
nuaugG_MP21=nu_G(:,5)-numinG_MP21;nuaugG_PM21=nu_G(:,6)-numinG_MP21;
nuaugG_MP22=nu_G(:,5)-numinG_MP22;nuaugG_PM22=nu_G(:,7)-numinG_MP22;
%Transient PP correlation
numinA_G12t=min(tnu_A(:,1),tnu_A(:,2)).*crG(1);
numinA_G23t=min(tnu_A(:,2),tnu_A(:,3)).*crG(1);
numinA_G13t=min(tnu_A(:,1),tnu_A(:,3)).*crG(1);
nuaugA_G12t=tnu_A(:,1)-numinA_G12t;nuaugA_G21t=tnu_A(:,2)-numinA_G12t;
nuaugA_G23t=tnu_A(:,2)-numinA_G23t;nuaugA_G32t=tnu_A(:,3)-numinA_G23t;
nuaugA_G13t=tnu_A(:,1)-numinA_G13t;nuaugA_G31t=tnu_A(:,3)-numinA_G13t;
numinA_Mt=min(tnu_A(:,4),tnu_A(:,5)).*crM(1);
nuaugA_M1t=tnu_A(:,4)-numinA_Mt;nuaugA_M2t=tnu_A(:,5)-numinA_Mt;
numinA_MP11t=min(tnu_A(:,4),tnu_A(:,6)).*crMP(1);
numinA_MP12t=min(tnu_A(:,4),tnu_A(:,7)).*crMP(1);
numinA_MP21t=min(tnu_A(:,5),tnu_A(:,6)).*crMP(1);
numinA_MP22t=min(tnu_A(:,5),tnu_A(:,7)).*crMP(1);
nuaugA_MP11t=tnu_A(:,4)-numinA_MP11t;nuaugA_PM11t=tnu_A(:,6)-numinA_MP11t;
nuaugA_MP12t=tnu_A(:,4)-numinA_MP12t;nuaugA_PM12t=tnu_A(:,7)-numinA_MP12t;
nuaugA_MP21t=tnu_A(:,5)-numinA_MP21t;nuaugA_PM21t=tnu_A(:,6)-numinA_MP21t;
nuaugA_MP22t=tnu_A(:,5)-numinA_MP22t;nuaugA_PM22t=tnu_A(:,7)-numinA_MP22t;
numinG_G12t=min(tnu_G(:,1),tnu_G(:,2)).*crG(1);
numinG_G23t=min(tnu_G(:,2),tnu_G(:,3)).*crG(1);
numinG_G13t=min(tnu_G(:,1),tnu_G(:,3)).*crG(1);
nuaugG_G12t=tnu_G(:,1)-numinG_G12t;nuaugG_G21t=tnu_G(:,2)-numinG_G12t;
nuaugG_G23t=tnu_G(:,2)-numinG_G23t;nuaugG_G32t=tnu_G(:,3)-numinG_G23t;
nuaugG_G13t=tnu_G(:,1)-numinG_G13t;nuaugG_G31t=tnu_G(:,3)-numinG_G13t;
numinG_Mt=min(tnu_G(:,4),tnu_G(:,5)).*crM(1);
nuaugG_M1t=tnu_G(:,4)-numinG_Mt;nuaugG_M2t=tnu_G(:,5)-numinG_Mt;
numinG_MP11t=min(tnu_G(:,4),tnu_G(:,6)).*crMP(1);
numinG_MP12t=min(tnu_G(:,4),tnu_G(:,7)).*crMP(1);
numinG_MP21t=min(tnu_G(:,5),tnu_G(:,6)).*crMP(1);
numinG_MP22t=min(tnu_G(:,5),tnu_G(:,7)).*crMP(1);
nuaugG_MP11t=tnu_G(:,4)-numinG_MP11t;nuaugG_PM11t=tnu_G(:,6)-numinG_MP11t;
nuaugG_MP12t=tnu_G(:,4)-numinG_MP12t;nuaugG_PM12t=tnu_G(:,7)-numinG_MP12t;
nuaugG_MP21t=tnu_G(:,5)-numinG_MP21t;nuaugG_PM21t=tnu_G(:,6)-numinG_MP21t;
nuaugG_MP22t=tnu_G(:,5)-numinG_MP22t;nuaugG_PM22t=tnu_G(:,7)-numinG_MP22t;

%Defined values
F=9.64853329e4; %Coulombs per Mole
T=273.15+35; %Kelvin
RTovF=8314.472*T/F; %R*1000 so in mV
thick_GP=0.2; %thickeness of membrane shell, 0.2micron for PGC and GC
thick_M=1; %thickeness of membrane shell, 1micron for Mitral Cell

% PARAMS
Iapp_G=10; %mV can be adjusted
Iapp_M=130; %mV can be adjusted
Iapp_P=50; %mV can be adjusted
C_G=2.0; %micro-F/cm^2
C_MP=1.2; %micro-F/cm^2
Rm_GM=30; %kilo-Ohm/cm^2
Rm_P=20; %kilo-Ohm/cm^2
w_G=3*ones(2,1);
w_Gc=0.3; %comm GC weight to MC
w_MGA=ones(3,1);
    w_MGA(2)=.5*w_MGA(2); %so comm GC fires less
w_MGN=ones(3,1);
    w_MGN(2)=.5*w_MGN(2);
w_MPA=ones(2,1);
w_MPN=ones(2,1);
w_P=2*ones(2,1); %weight = 8 Li & Cleland (2013)
tau_sGP=5.5; %ms for PGC and GC
tau_sM=10; %ms for MC

%maximal conduct
gNa_GP=70;
gNa_M=120; %mS/cm^2
gNaP=0.42;
gDR_GP=25;
gDR_M=70;
gM_G=0.5; %only in soma
gM_P=1.0; %only in soma
gA_G=80;
gA_M=10; %only in soma
gA_P=40;
gKS=84;
gH=0.2; %only in dendrite and spine
gCaL=0.85;
gCaPN_G=0.2; %only in dendrite and spine
gCaPN_P=1.0; %only in dendrite and spine
gCaT_G=0.1; %only in dendrite and spine
gCaT_P=3.0; %only in dendrite and spine
gCAN=1.0; %only in dendrite and spine
gKCa_G=0.5; %only in dendrite and spine
gKCa_M=5; %only in soma
gKCa_P=2.0; %only in dendrite and spine
g_GABA_G=1.5; %ns
g_AMPA=2; %ns
g_NMDA=1;
g_GABA_P=2; %ns
%reversal Poten
El_GM=-60; %mV
El_P=-65; %mV
Ek=-80;
Ena=45;
Ecat=10;
Eh=0;
Esyn_GP=-80;
Esyn_M=0;
% keep track of tLstSp
tLstSp_G=zeros(1,3);
tLstSp_M=zeros(1,2);
tLstSp_P=zeros(1,2);
idGspk=zeros(1,3);
idMspk=zeros(1,2);
idPspk=zeros(1,2);
mnTim = round(10/dt)+1; %min # time steps for another spike
%Initial variables
vPt_G=-70*ones(1,3);
vPt_M=-70*ones(1,2);
vPt_P=-70*ones(1,2);
vltThres=10; %(mV) voltage threshold for spike
caCon_G=0.05*ones(1,3); %only in soma, micrMol/l
caCon_M=0.05*ones(1,2); %only in soma, micrMol/l
caCon_P=0.05*ones(1,2); %only in soma, micrMol/l
Eca_G=70*ones(1,3); %dyn variable
Eca_M=70*ones(1,2); %dyn variable
Eca_P=70*ones(1,2); %dyn variable
xVr_Gt=zeros(12,3); %gating variables
xVr_Mt=zeros(11,2); %gating variables
xVr_Pt=zeros(12,2); %gating variables 
xiR1_G=zeros(12,3); %aux for RK4
xiR2_G=zeros(12,3); %aux for RK4
xiR3_G=zeros(12,3); %aux for RK4
xiR4_G=zeros(12,3); %aux for RK4
xiR1_M=zeros(11,2); %aux for RK4
xiR2_M=zeros(11,2); %aux for RK4
xiR3_M=zeros(11,2); %aux for RK4
xiR4_M=zeros(11,2); %aux for RK4
xiR1_P=zeros(12,2); %aux for RK4
xiR2_P=zeros(12,2); %aux for RK4
xiR3_P=zeros(12,2); %aux for RK4
xiR4_P=zeros(12,2); %aux for RK4
s_GABA_G=zeros(1,3);
s_AMPA=zeros(1,2);
s_NMDA=zeros(1,2);
s_GABA_P=zeros(1,2);
sA_GC=zeros(1,3);
sA_MC=zeros(1,2);
sA_PGC=zeros(1,2);
sG_GC=zeros(1,3);
sG_MC=zeros(1,2);
sG_PGC=zeros(1,2);
%Generating random distribution array for Poisson kicks
randA_Gt=rand(Ltrn,3);
randA_Gct=rand(Ltrn,3);
randA_Mt=rand(Ltrn,2);
randA_Mct=rand(Ltrn,1);
randA_MPct=rand(Ltrn,2);
randA_PMt=rand(Ltrn,2);
randG_Gt=rand(Ltrn,3);
randG_Gct=rand(Ltrn,3);
randG_Mt=rand(Ltrn,2);
randG_Mct=rand(Ltrn,2);
randG_MPct=rand(Ltrn,2);
randG_PMt=rand(Ltrn,2);
%transient time loop, overwrite voltages & state vars
for j=1:Ltrn    
    %Synaptic currents (GABA/NMDA/AMPA) aprx by Forward-Euler
    s_GABA_G=s_GABA_G+dt*dsGABA(s_GABA_G,vPt_G);
    Isyn_GC=[g_GABA_G*(vPt_M(1)-Esyn_GP)*(w_G(1)*s_GABA_G(1)+w_Gc*s_GABA_G(2)),...
        g_GABA_G*(vPt_M(2)-Esyn_GP)*(w_Gc*s_GABA_G(2)+w_G(2)*s_GABA_G(3))];
    s_AMPA=s_AMPA+dt*dsAMPA(s_AMPA,vPt_M);
    s_NMDA=s_NMDA+dt*dsNMDA(s_NMDA,vPt_M);
    Isyn_MCG=[w_MGA(1)*g_AMPA*s_AMPA(1)*(vPt_G(1)-Esyn_M)+w_MGN(1)*g_NMDA*s_NMDA(1)*B(vPt_G(1))*(vPt_G(1)-Esyn_M),...
        w_MGA(2)*g_AMPA*(vPt_G(2)-Esyn_M)*(s_AMPA(1)+s_AMPA(2))+w_MGN(2)*g_NMDA*B(vPt_G(2))*(vPt_G(2)-Esyn_M)*(s_NMDA(1)+s_NMDA(2)),...
        w_MGA(3)*g_AMPA*s_AMPA(2)*(vPt_G(3)-Esyn_M)+w_MGN(3)*g_NMDA*s_NMDA(2)*B(vPt_G(3))*(vPt_G(3)-Esyn_M)];
    Isyn_MCP=[w_MPA(1)*g_AMPA*s_AMPA(1)*(vPt_P(1)-Esyn_M)+w_MPN(1)*g_NMDA*s_NMDA(1)*B(vPt_P(1))*(vPt_P(1)-Esyn_M),...
        w_MPA(2)*g_AMPA*s_AMPA(2)*(vPt_P(2)-Esyn_M)+w_MPN(2)*g_NMDA*s_NMDA(2)*B(vPt_P(2))*(vPt_P(2)-Esyn_M)];
    s_GABA_P=s_GABA_P+dt*dsGABA(s_GABA_P,vPt_P);
    Isyn_PGC=[w_P(1)*g_GABA_P*s_GABA_P(1).*(vPt_M(1)-Esyn_GP),...
        w_P(2)*g_GABA_P*s_GABA_P(2).*(vPt_M(2)-Esyn_GP)];
    
    %Background noise Poisson kicks with AMPA synapse
    sA_GC=sA_GC+(dt/tau_sGP)*-sA_GC;
    sA_GC=[sA_GC(1)+atrn(1)*(randA_Gt(j,1)<nuaugA_G12t*dt)+atrn(1)*(randA_Gt(j,1)<nuaugA_G13t*dt)+atrn(1)*(randA_Gct(j,1)<numinA_G12t*dt)+atrn(1)*(randA_Gct(j,3)<numinA_G13t*dt),...
        sA_GC(2)+atrn(2)*(randA_Gt(j,2)<nuaugA_G21t*dt)+atrn(2)*(randA_Gt(j,2)<nuaugA_G23t*dt)+atrn(2)*(randA_Gct(j,1)<numinA_G12t*dt)+1*(randA_Gct(j,2)<numinA_G23t*dt),...
        sA_GC(3)+atrn(3)*(randA_Gt(j,3)<nuaugA_G31t*dt)+atrn(3)*(randA_Gt(j,3)<nuaugA_G32t*dt)+atrn(3)*(randA_Gct(j,2)<numinA_G23t*dt)+1*(randA_Gct(j,3)<numinA_G13t*dt)];
    xiA_GC=sA_GC.*(vPt_G-Esyn_M);
    sA_MC=sA_MC+(dt/tau_sM)*-sA_MC;
    sA_MC=[sA_MC(1)+atrn(4)*(randA_Mt(j,1)<nuaugA_M1t*dt)+atrn(4)*(randA_Mct(j)<numinA_Mt*dt)+atrn(4)*(randA_Mt(j,1)<nuaugA_MP11t*dt)+atrn(4)*(randA_Mt(j,1)<nuaugA_MP12t*dt)+atrn(4)*(randA_MPct(j,1)<numinA_MP11t*dt)+atrn(4)*(randA_MPct(j,1)<numinA_MP12t*dt),...
        sA_MC(2)+atrn(5)*(randA_Mt(j,2)<nuaugA_M2t*dt)+atrn(5)*(randA_Mct(j)<numinA_Mt*dt)+atrn(5)*(randA_Mt(j,2)<nuaugA_MP21t*dt)+atrn(5)*(randA_Mt(j,2)<nuaugA_MP22t*dt)+atrn(5)*(randA_MPct(j,2)<numinA_MP21t*dt)+atrn(5)*(randA_MPct(j,2)<numinA_MP22t*dt)];
    xiA_MC=sA_MC.*(vPt_M-Esyn_M);
    sA_PGC=sA_PGC+(dt/tau_sGP)*-sA_PGC;
    sA_PGC=[sA_PGC(1)+atrn(6)*(randA_PMt(j,1)<nuaugA_PM11t*dt)+atrn(6)*(randA_PMt(j,1)<nuaugA_PM21t*dt)+atrn(6)*(randA_MPct(j,1)<numinA_MP11t*dt)+atrn(6)*(randA_MPct(j,2)<numinA_MP21t*dt),...
        sA_PGC(2)+atrn(7)*(randA_PMt(j,2)<nuaugA_PM12t*dt)+atrn(7)*(randA_PMt(j,2)<nuaugA_PM22t*dt)+atrn(7)*(randA_MPct(j,1)<numinA_MP12t*dt)+atrn(7)*(randA_MPct(j,2)<numinA_MP22t*dt)];
    xiA_PGC=sA_PGC.*(vPt_P-Esyn_M);
    %Background noise Poisson kicks with GABA synapse
    sG_GC=sG_GC+(dt/tau_sGP)*-sG_GC;
    sG_GC=[sG_GC(1)+atrn(1)*(randG_Gt(j,1)<nuaugG_G12t*dt)+atrn(1)*(randG_Gt(j,1)<nuaugG_G13t*dt)+atrn(1)*(randG_Gct(j,1)<numinG_G12t*dt)+atrn(1)*(randG_Gct(j,3)<numinG_G13t*dt),...
        sG_GC(2)+atrn(2)*(randG_Gt(j,2)<nuaugG_G21t*dt)+atrn(2)*(randG_Gt(j,2)<nuaugG_G23t*dt)+atrn(2)*(randG_Gct(j,1)<numinG_G12t*dt)+1*(randG_Gct(j,2)<numinG_G23t*dt),...
        sG_GC(3)+atrn(3)*(randG_Gt(j,3)<nuaugG_G31t*dt)+atrn(3)*(randG_Gt(j,3)<nuaugG_G32t*dt)+atrn(3)*(randG_Gct(j,2)<numinG_G23t*dt)+1*(randG_Gct(j,3)<numinG_G13t*dt)];
    xiG_GC=sG_GC.*(vPt_G-Esyn_GP);
    sG_MC=sG_MC+(dt/tau_sM)*-sG_MC;
    sG_MC=[sG_MC(1)+atrn(4)*(randG_Mt(j,1)<nuaugG_M1t*dt)+atrn(4)*(randG_Mct(j)<numinG_Mt*dt)+atrn(4)*(randG_Mt(j,1)<nuaugG_MP11t*dt)+atrn(4)*(randG_Mt(j,1)<nuaugG_MP12t*dt)+atrn(4)*(randG_MPct(j,1)<numinG_MP11t*dt)+atrn(4)*(randG_MPct(j,1)<numinG_MP12t*dt),...
        sG_MC(2)+atrn(5)*(randG_Mt(j,2)<nuaugG_M2t*dt)+atrn(5)*(randG_Mct(j)<numinG_Mt*dt)+atrn(5)*(randG_Mt(j,2)<nuaugG_MP21t*dt)+atrn(5)*(randG_Mt(j,2)<nuaugG_MP22t*dt)+atrn(5)*(randG_MPct(j,2)<numinG_MP21t*dt)+atrn(5)*(randG_MPct(j,2)<numinG_MP22t*dt)];
    xiG_MC=sG_MC.*(vPt_M-Esyn_GP);
    sG_PGC=sG_PGC+(dt/tau_sGP)*-sG_PGC;
    sG_PGC=[sG_PGC(1)+atrn(6)*(randG_PMt(j,1)<nuaugG_PM11t*dt)+atrn(6)*(randG_PMt(j,1)<nuaugG_PM21t*dt)+atrn(6)*(randG_MPct(j,1)<numinG_MP11t*dt)+atrn(6)*(randG_MPct(j,2)<numinG_MP21t*dt),...
        sG_PGC(2)+atrn(7)*(randG_PMt(j,2)<nuaugG_PM12t*dt)+atrn(7)*(randG_PMt(j,2)<nuaugG_PM22t*dt)+atrn(7)*(randG_MPct(j,1)<numinG_MP12t*dt)+atrn(7)*(randG_MPct(j,2)<numinG_MP22t*dt)];
    xiG_PGC=sG_PGC.*(vPt_P-Esyn_GP);
    %Sum background noise
    xi_GC=xiA_GC+xiG_GC;
    xi_MC=xiA_MC+xiG_MC;
    xi_PGC=xiA_PGC+xiG_PGC;
    
    % --- step 1 ---
    rhsV1_G=1/Rm_GM*(vPt_G-El_GM)+gNa_GP*xVr_Gt(1,:).^3.*xVr_Gt(2,:).*(vPt_G-Ena)+...
        gA_G*xVr_Gt(3,:).*xVr_Gt(4,:).*(vPt_G-Ek)+gM_G.*xVr_Gt(5,:).*(vPt_G-Ek)+...
        gCaPN_G.*xVr_Gt(6,:).^2.*xVr_Gt(7,:).*(vPt_G-Eca_G)+gCaT_G.*xVr_Gt(8,:).^2.*xVr_Gt(9,:).*(vPt_G-Eca_G)...
        +gCAN*(caCon_G/(200+caCon_G)).*xVr_Gt(10,:).*(vPt_G-Ecat)+gKCa_G*xVr_Gt(11,:).*(vPt_G-Ek)...
        +gDR_GP*xVr_Gt(12,:).*(vPt_G-Ek);
    rhsV1_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV1_G-xi_GC);
    rhsV1_M=1/Rm_GM*(vPt_M-El_GM)+gNa_M*xVr_Mt(1,:).^3.*xVr_Mt(2,:).*(vPt_M-Ena)+gNaP*minf(vPt_M).*(vPt_M-Ena)+...
        gA_M*xVr_Mt(3,:).*xVr_Mt(4,:).*(vPt_M-Ek)+gKS*xVr_Mt(5,:).*xVr_Mt(6,:).*(vPt_M-Ek)+...
        gCaL*xVr_Mt(7,:).*xVr_Mt(8,:).*(vPt_M-Eca_M)+gKCa_M*xVr_Mt(9,:).*(vPt_M-Ek)+gDR_M*xVr_Mt(10,:).^2.*xVr_Mt(11,:).*(vPt_M-Ek);
    rhsV1_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV1_M-xi_MC);
    rhsV1_P=1/Rm_P*(vPt_P-El_P)+gNa_GP*xVr_Pt(1,:).^3.*xVr_Pt(2,:).*(vPt_P-Ena)+...
        gA_P*xVr_Pt(3,:).*xVr_Pt(4,:).*(vPt_P-Ek)+gM_P.*xVr_Pt(5,:).*(vPt_P-Ek)+gH.*xVr_Pt(6,:).*(vPt_P-Eh)+...
        gCaPN_P.*xVr_Pt(7,:).^2.*xVr_Pt(8,:).*(vPt_P-Eca_P)+gCaT_P.*xVr_Pt(9,:).^2.*xVr_Pt(10,:).*(vPt_P-Eca_P)...
        +gKCa_P*xVr_Pt(11,:).*(vPt_P-Ek)+gDR_GP*xVr_Pt(12,:).*(vPt_P-Ek);
    rhsV1_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV1_P-xi_PGC);
    
    xiR1_G(1,:)=dt*2.1./tauM(vPt_G).*(mInf(vPt_G)-xVr_Gt(1,:)); %I_Na
    xiR1_G(2,:)=dt*2.1./tauH(vPt_G).*(hInf(vPt_G)-xVr_Gt(2,:)); %I_Na
    xiR1_G(3,:)=dt*3.3./atau(vPt_G).*(aminf_g(vPt_G)-xVr_Gt(3,:)); %I_A
    xiR1_G(4,:)=dt*3.3./hatau_g(vPt_G).*(ahinf_g(vPt_G)-xVr_Gt(4,:)); %I_A
    xiR1_G(5,:)=dt./tauMusc(vPt_G).*(minfMusc(vPt_G)-xVr_Gt(5,:)); %I_M
    xiR1_G(6,:)=dt./capmtau(vPt_G).*(capminf(vPt_G)-xVr_Gt(6,:)); %I_CaPN
    xiR1_G(7,:)=dt./caphtau(vPt_G).*(caphinf(vPt_G)-xVr_Gt(7,:)); %I_CaPN
    xiR1_G(8,:)=dt./cattaumG(vPt_G).*(catminfG(vPt_G)-xVr_Gt(8,:)); %I_CaT
    xiR1_G(9,:)=dt./cattauh(vPt_G).*(cathinf(vPt_G)-xVr_Gt(9,:)); %I_CaT
    xiR1_G(10,:)=dt./cantau(vPt_G).*(canminf(vPt_G)-xVr_Gt(10,:)); %I_CAN
    xiR1_G(11,:)=dt*(kcaa(vPt_G,caCon_G).*(1-xVr_Gt(11,:)) - 0.05*xVr_Gt(11,:)); %I_KCa
    xiR1_G(12,:)=dt*3.3*(mss_dr(vPt_G)-xVr_Gt(12,:))./taum_dr(vPt_G); %I_DR
    
    xiR1_M(1,:)=dt*(ana(vPt_M).*(1-xVr_Mt(1,:)) - bna(vPt_M).*xVr_Mt(1,:)); %I_Na
    xiR1_M(2,:)=dt*(ahna(vPt_M).*(1-xVr_Mt(2,:)) - bhna(vPt_M).*xVr_Mt(2,:)); %I_Na
    xiR1_M(3,:)=dt*3.3./atau(vPt_M).*(aminf(vPt_M)-xVr_Mt(3,:)); %I_A
    xiR1_M(4,:)=dt*3.3./hatau(vPt_M).*(ahinf(vPt_M)-xVr_Mt(4,:)); %I_A
    xiR1_M(5,:)=dt/10*(ksminf(vPt_M)-xVr_Mt(5,:)); %I_KS
    xiR1_M(6,:)=dt./kshtau(vPt_M).*(kshinf(vPt_M)-xVr_Mt(6,:));
    xiR1_M(7,:)=dt*(cala(vPt_M).*(1-xVr_Mt(7,:)) - calb(vPt_M).*xVr_Mt(7,:)); %I_CaL
    xiR1_M(8,:)=dt*(calha(vPt_M).*(1-xVr_Mt(8,:)) - calhb(vPt_M).*xVr_Mt(8,:));
    xiR1_M(9,:)=dt*(kcaa(vPt_M,caCon_M).*(1-xVr_Mt(9,:)) - 0.05*xVr_Mt(9,:));
    xiR1_M(10,:)=dt*(nss_dr(vPt_M)-xVr_Mt(10,:))./taun_dr(vPt_M);    %I_DR
    xiR1_M(11,:)=dt*(kss_dr(vPt_M)-xVr_Mt(11,:))./50;
    
    xiR1_P(1,:)=dt*2.1./tauM(vPt_P).*(mInf(vPt_P)-xVr_Pt(1,:)); %I_Na
    xiR1_P(2,:)=dt*2.1./tauH(vPt_P).*(hInf(vPt_P)-xVr_Pt(2,:)); %I_Na
    xiR1_P(3,:)=dt*3.3./atau(vPt_P).*(aminf_g(vPt_P)-xVr_Pt(3,:)); %I_A
    xiR1_P(4,:)=dt*3.3./hatau_g(vPt_P).*(ahinf_g(vPt_P)-xVr_Pt(4,:)); %I_A
    xiR1_P(5,:)=dt./tauMusc(vPt_P).*(minfMusc(vPt_P)-xVr_Pt(5,:)); %I_M
    xiR1_P(6,:)=dt*2.1./hcurmtau(vPt_P).*(hcurminf(vPt_P)-xVr_Pt(6,:)); %I_H
    xiR1_P(7,:)=dt./capmtau(vPt_P).*(capminf(vPt_P)-xVr_Pt(7,:)); %I_CaPN
    xiR1_P(8,:)=dt./caphtau(vPt_P).*(caphinf(vPt_P)-xVr_Pt(8,:)); %I_CaPN
    xiR1_P(9,:)=dt./cattaumP(vPt_P).*(catminfP(vPt_P)-xVr_Pt(9,:)); %I_CaT
    xiR1_P(10,:)=dt./cattauh(vPt_P).*(cathinf(vPt_P)-xVr_Pt(10,:)); %I_CaT
    xiR1_P(11,:)=dt*(kcaa(vPt_P,caCon_P).*(1-xVr_Pt(11,:)) - 0.05*xVr_Pt(11,:)); %I_KCa
    xiR1_P(12,:)=dt*3.3*(mss_dr(vPt_P)-xVr_Pt(12,:))./taum_dr(vPt_P); %I_DR

    %update calcium concentration
    Ica_G=gCaPN_G*xVr_Gt(6,:).^2.*xVr_Gt(7,:).*(vPt_G-Eca_G)+gCaT_G*xVr_Gt(8,:).^2.*xVr_Gt(9,:).*(vPt_G-Eca_G);
    caCr1_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-caCon_G)/800); %tauCa=800ms for PGC & GC
    
    Ica_M=gCaL*xVr_Mt(7,:).*xVr_Mt(8,:).*(vPt_M-Eca_M);
    caCr1_M=dt*(-Ica_M/(2*F*thick_M) + (.05-caCon_M)/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*xVr_Pt(7,:).^2.*xVr_Pt(8,:).*(vPt_P-Eca_P)+gCaT_P*xVr_Pt(9,:).^2.*xVr_Pt(10,:).*(vPt_P-Eca_P);
    caCr1_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-caCon_P)/800); %tauCa=800ms for PGC & GC
    
    % --- step 2 ---
    rhsV2_G=1/Rm_GM*(vPt_G+.5*rhsV1_G-El_GM)+gNa_GP*(xVr_Gt(1,:)+.5*xiR1_G(1,:)).^3.*(xVr_Gt(2,:)+.5*xiR1_G(2,:)).*(vPt_G+.5*rhsV1_G-Ena)+...
        gA_G*(xVr_Gt(3,:)+.5*xiR1_G(3,:)).*(xVr_Gt(4,:)+.5*xiR1_G(4,:)).*(vPt_G+.5*rhsV1_G-Ek)+gM_G.*(xVr_Gt(5,:)+.5*xiR1_G(5,:)).*(vPt_G+.5*rhsV1_G-Ek)+...
        gCaPN_G.*(xVr_Gt(6,:)+.5*xiR1_G(6,:)).^2.*(xVr_Gt(7,:)+.5*xiR1_G(7,:)).*(vPt_G+.5*rhsV1_G-Eca_G)+gCaT_G.*(xVr_Gt(8,:)+.5*xiR1_G(8,:)).^2.*(xVr_Gt(9,:)+.5*xiR1_G(9,:)).*(vPt_G+.5*rhsV1_G-Eca_G)...
        +gCAN.*((caCon_G+.5*caCr1_G)/(200+(caCon_G+.5*caCr1_G))).*(xVr_Gt(10,:)+.5*xiR1_G(10,:)).*(vPt_G+.5*rhsV1_G-Ecat)+gKCa_G*(xVr_Gt(11,:)+.5*xiR1_G(11,:)).*(vPt_G+.5*rhsV1_G-Ek)...
        +gDR_GP*(xVr_Gt(12,:)+.5*xiR1_G(12,:)).*(vPt_G+.5*rhsV1_G-Ek);
    rhsV2_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV2_G-xi_GC);
    rhsV2_M=1/Rm_GM*(vPt_M+.5*rhsV1_M-El_GM)+gNa_M*(xVr_Mt(1,:)+.5*xiR1_M(1,:)).^3.*(xVr_Mt(2,:)+.5*xiR1_M(2,:)).*(vPt_M+.5*rhsV1_M-Ena)+...
        gNaP*minf(vPt_M+.5*rhsV1_M).*(vPt_M+.5*rhsV1_M-Ena)+gA_M*(xVr_Mt(3,:)+.5*xiR1_M(3,:)).*(xVr_Mt(4,:)+...
        .5*xiR1_M(4,:)).*(vPt_M+.5*rhsV1_M-Ek)+gKS*(xVr_Mt(5,:)+.5*xiR1_M(5,:)).*(xVr_Mt(6,:)+.5*xiR1_M(6,:)).*(vPt_M+.5*rhsV1_M-Ek)+...
        gCaL*(xVr_Mt(7,:)+.5*xiR1_M(7,:)).*(xVr_Mt(8,:)+.5*xiR1_M(8,:)).*(vPt_M+.5*rhsV1_M-Eca_M)+gKCa_M*(xVr_Mt(9,:)+...
        .5*xiR1_M(9,:)).*(vPt_M+.5*rhsV1_M-Ek)+gDR_M*(xVr_Mt(10,:)+.5*xiR1_M(10,:)).^2.*(xVr_Mt(11,:)+.5*xiR1_M(11,:)).*(vPt_M+.5*rhsV1_M-Ek);
    rhsV2_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV2_M-xi_MC);
    rhsV2_P=1/Rm_P*(vPt_P+.5*rhsV1_P-El_P)+gNa_GP*(xVr_Pt(1,:)+.5*xiR1_P(1,:)).^3.*(xVr_Pt(2,:)+.5*xiR1_P(2,:)).*(vPt_P+.5*rhsV1_P-Ena)+...
        gA_P*(xVr_Pt(3,:)+.5*xiR1_P(3,:)).*(xVr_Pt(4,:)+.5*xiR1_P(4,:)).*(vPt_P+.5*rhsV1_P-Ek)+gM_P.*(xVr_Pt(5,:)+.5*xiR1_P(5,:)).*(vPt_P+.5*rhsV1_P-Ek)...
        +gH.*(xVr_Pt(6,:)+.5*xiR1_P(6,:)).*(vPt_P+.5*rhsV1_P-Eh)+gCaPN_P.*(xVr_Pt(7,:)+.5*xiR1_P(7,:)).^2.*(xVr_Pt(8,:)+.5*xiR1_P(8,:)).*(vPt_P+.5*rhsV1_P-Eca_P)...
        +gCaT_P.*(xVr_Pt(9,:)+.5*xiR1_P(9,:)).^2.*(xVr_Pt(10,:)+.5*xiR1_P(10,:)).*(vPt_P+.5*rhsV1_P-Eca_P)...
        +gKCa_P*(xVr_Pt(11,:)+.5*xiR1_P(11,:)).*(vPt_P+.5*rhsV1_P-Ek)+gDR_GP*(xVr_Pt(12,:)+.5*xiR1_P(12,:)).*(vPt_P+.5*rhsV1_P-Ek);
    rhsV2_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV2_P-xi_PGC);
    
    xiR2_G(1,:)=dt*2.1./tauM(vPt_G+.5*rhsV1_G).*(mInf(vPt_G+.5*rhsV1_G)-(xVr_Gt(1,:)+.5*xiR1_G(1,:))); %I_Na
    xiR2_G(2,:)=dt*2.1./tauH(vPt_G+.5*rhsV1_G).*(hInf(vPt_G+.5*rhsV1_G)-(xVr_Gt(2,:)+.5*xiR1_G(2,:))); %I_Na
    xiR2_G(3,:)=dt*3.3./atau(vPt_G+.5*rhsV1_G).*(aminf_g(vPt_G+.5*rhsV1_G)-(xVr_Gt(3,:)+.5*xiR1_G(3,:))); %I_A
    xiR2_G(4,:)=dt*3.3./hatau_g(vPt_G+.5*rhsV1_G).*(ahinf_g(vPt_G+.5*rhsV1_G)-(xVr_Gt(4,:)+.5*xiR1_G(4,:))); %I_A
    xiR2_G(5,:)=dt./tauMusc(vPt_G+.5*rhsV1_G).*(minfMusc(vPt_G+.5*rhsV1_G)-(xVr_Gt(5,:)+.5*xiR1_G(5,:))); %I_M
    xiR2_G(6,:)=dt./capmtau(vPt_G+.5*rhsV1_G).*(capminf(vPt_G+.5*rhsV1_G)-(xVr_Gt(6,:)+.5*xiR1_G(6,:))); %I_CaPN
    xiR2_G(7,:)=dt./caphtau(vPt_G+.5*rhsV1_G).*(caphinf(vPt_G+.5*rhsV1_G)-(xVr_Gt(7,:)+.5*xiR1_G(7,:))); %I_CaPN
    xiR2_G(8,:)=dt./cattaumG(vPt_G+.5*rhsV1_G).*(catminfG(vPt_G+.5*rhsV1_G)-(xVr_Gt(8,:)+.5*xiR1_G(8,:))); %I_CaT
    xiR2_G(9,:)=dt./cattauh(vPt_G+.5*rhsV1_G).*(cathinf(vPt_G+.5*rhsV1_G)-(xVr_Gt(9,:)+.5*xiR1_G(9,:))); %I_CaT
    xiR2_G(10,:)=dt./cantau(vPt_G+.5*rhsV1_G).*(canminf(vPt_G+.5*rhsV1_G)-(xVr_Gt(10,:)+.5*xiR1_G(10,:))); %I_CAN
    xiR2_G(11,:)=dt*(kcaa(vPt_G+.5*rhsV1_G,caCon_G+.5*caCr1_G).*(1-(xVr_Gt(11,:)+.5*xiR1_G(11,:))) - 0.05*(xVr_Gt(11,:)+.5*xiR1_G(11,:))); %I_KCa
    xiR2_G(12,:)=dt*3.3*(mss_dr(vPt_G+.5*rhsV1_G)-(xVr_Gt(12,:)+.5*xiR1_G(12,:)))./taum_dr(vPt_G+.5*rhsV1_G); %I_DR
    
    xiR2_M(1,:)=dt*(ana(vPt_M+.5*rhsV1_M).*(1-(xVr_Mt(1,:)+.5*xiR1_M(1,:))) - bna(vPt_M+.5*rhsV1_M).*(xVr_Mt(1,:)+.5*xiR1_M(1,:)));
    xiR2_M(2,:)=dt*(ahna(vPt_M+.5*rhsV1_M).*(1-(xVr_Mt(2,:)+.5*xiR1_M(2,:))) - bhna(vPt_M+.5*rhsV1_M).*(xVr_Mt(2,:)+.5*xiR1_M(2,:)));
    xiR2_M(3,:)=dt*3.3./atau(vPt_M+.5*rhsV1_M).*(aminf(vPt_M+.5*rhsV1_M)-(xVr_Mt(3,:)+.5*xiR1_M(3,:))); %I_A
    xiR2_M(4,:)=dt*3.3./hatau(vPt_M+.5*rhsV1_M).*(ahinf(vPt_M+.5*rhsV1_M)-(xVr_Mt(4,:)+.5*xiR1_M(4,:))); %I_A
    xiR2_M(5,:)=dt/10*(ksminf(vPt_M+.5*rhsV1_M)-(xVr_Mt(5,:)+.5*xiR1_M(5,:))); %I_KS
    xiR2_M(6,:)=dt./kshtau(vPt_M+.5*rhsV1_M).*(kshinf(vPt_M+.5*rhsV1_M)-(xVr_Mt(6,:)+.5*xiR1_M(6,:)));
    xiR2_M(7,:)=dt*(cala(vPt_M+.5*rhsV1_M).*(1-(xVr_Mt(7,:)+.5*xiR1_M(7,:))) - calb(vPt_M+.5*rhsV1_M).*(xVr_Mt(7,:)+.5*xiR1_M(7,:))); %I_CaL
    xiR2_M(8,:)=dt*(calha(vPt_M+.5*rhsV1_M).*(1-(xVr_Mt(8,:)+.5*xiR1_M(8,:))) - calhb(vPt_M+.5*rhsV1_M).*(xVr_Mt(8,:)+.5*xiR1_M(8,:)));
    xiR2_M(9,:)=dt*(kcaa(vPt_M+.5*rhsV1_M,caCon_M+.5*caCr1_M).*(1-(xVr_Mt(9,:)+.5*xiR1_M(9,:))) - 0.05*(xVr_Mt(9,:)+.5*xiR1_M(9,:)));
    xiR2_M(10,:)=dt*(nss_dr(vPt_M+.5*rhsV1_M)-(xVr_Mt(10,:)+.5*xiR1_M(10,:)))./taun_dr(vPt_M+.5*rhsV1_M);    %I_DR
    xiR2_M(11,:)=dt*(kss_dr(vPt_M+.5*rhsV1_M)-(xVr_Mt(11,:)+.5*xiR1_M(11,:)))./50;
    
    xiR2_P(1,:)=dt*2.1./tauM(vPt_P+.5*rhsV1_P).*(mInf(vPt_P+.5*rhsV1_P)-(xVr_Pt(1,:)+.5*xiR1_P(1,:))); %I_Na
    xiR2_P(2,:)=dt*2.1./tauH(vPt_P+.5*rhsV1_P).*(hInf(vPt_P+.5*rhsV1_P)-(xVr_Pt(2,:)+.5*xiR1_P(2,:))); %I_Na
    xiR2_P(3,:)=dt*3.3./atau(vPt_P+.5*rhsV1_P).*(aminf_g(vPt_P+.5*rhsV1_P)-(xVr_Pt(3,:)+.5*xiR1_P(3,:))); %I_A
    xiR2_P(4,:)=dt*3.3./hatau_g(vPt_P+.5*rhsV1_P).*(ahinf_g(vPt_P+.5*rhsV1_P)-(xVr_Pt(4,:)+.5*xiR1_P(4,:))); %I_A
    xiR2_P(5,:)=dt./tauMusc(vPt_P+.5*rhsV1_P).*(minfMusc(vPt_P+.5*rhsV1_P)-(xVr_Pt(5,:)+.5*xiR1_P(5,:))); %I_M
    xiR2_P(6,:)=dt*2.1./hcurmtau(vPt_P+.5*rhsV1_P).*(hcurminf(vPt_P+.5*rhsV1_P)-(xVr_Pt(6,:))+.5*xiR1_P(6,:)); %I_H
    xiR2_P(7,:)=dt./capmtau(vPt_P+.5*rhsV1_P).*(capminf(vPt_P+.5*rhsV1_P)-(xVr_Pt(7,:)+.5*xiR1_P(7,:))); %I_CaPN
    xiR2_P(8,:)=dt./caphtau(vPt_P+.5*rhsV1_P).*(caphinf(vPt_P+.5*rhsV1_P)-(xVr_Pt(8,:)+.5*xiR1_P(8,:))); %I_CaPN
    xiR2_P(9,:)=dt./cattaumP(vPt_P+.5*rhsV1_P).*(catminfP(vPt_P+.5*rhsV1_P)-(xVr_Pt(9,:)+.5*xiR1_P(9,:))); %I_CaT
    xiR2_P(10,:)=dt./cattauh(vPt_P+.5*rhsV1_P).*(cathinf(vPt_P+.5*rhsV1_P)-(xVr_Pt(10,:)+.5*xiR1_P(10,:))); %I_CaT
    xiR2_P(11,:)=dt*(kcaa(vPt_P+.5*rhsV1_P,caCon_P+.5*caCr1_P).*(1-(xVr_Pt(11,:)+.5*xiR1_P(11,:))) - 0.05*(xVr_Pt(11,:)+.5*xiR1_P(11,:))); %I_KCa
    xiR2_P(12,:)=dt*3.3*(mss_dr(vPt_P+.5*rhsV1_P)-(xVr_Pt(12,:)+.5*xiR1_P(12,:)))./taum_dr(vPt_P+.5*rhsV1_P); %I_DR

    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_Gt(6,:)+.5*xiR1_G(6,:)).^2.*(xVr_Gt(7,:)+.5*xiR1_G(7,:)).*(vPt_G+.5*rhsV1_G-Eca_G)+gCaT_G*(xVr_Gt(8,:)+.5*xiR1_G(8,:)).^2.*(xVr_Gt(9,:)+.5*xiR1_G(9,:)).*(vPt_G+.5*rhsV1_G-Eca_G);
    caCr2_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+.5*caCr1_G))/800); %tauCa=800ms for PGC & GC

    Ica_M=gCaL*(xVr_Mt(7,:)+.5*xiR1_M(7,:)).*(xVr_Mt(8,:)+.5*xiR1_M(8,:)).*(vPt_M+.5*rhsV1_M-Eca_M);
    caCr2_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr1_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_Pt(7,:)+.5*xiR1_P(7,:)).^2.*(xVr_Pt(8,:)+.5*xiR1_P(8,:)).*(vPt_P+.5*rhsV1_P-Eca_P)+gCaT_P*(xVr_Pt(9,:)+.5*xiR1_P(9,:)).^2.*(xVr_Pt(10,:)+.5*xiR1_P(10,:)).*(vPt_P+.5*rhsV1_P-Eca_P);
    caCr2_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+.5*caCr1_P))/800); %tauCa=800ms for PGC & GC
    
    % --- step 3 ---
    rhsV3_G=1/Rm_GM*(vPt_G+.5*rhsV2_G-El_GM)+gNa_GP*(xVr_Gt(1,:)+.5*xiR2_G(1,:)).^3.*(xVr_Gt(2,:)+.5*xiR2_G(2,:)).*(vPt_G+.5*rhsV2_G-Ena)+...
        gA_G*(xVr_Gt(3,:)+.5*xiR2_G(3,:)).*(xVr_Gt(4,:)+.5*xiR2_G(4,:)).*(vPt_G+.5*rhsV2_G-Ek)+gM_G.*(xVr_Gt(5,:)+.5*xiR2_G(5,:)).*(vPt_G+.5*rhsV2_G-Ek)+...
        gCaPN_G.*(xVr_Gt(6,:)+.5*xiR2_G(6,:)).^2.*(xVr_Gt(7,:)+.5*xiR2_G(7,:)).*(vPt_G+.5*rhsV2_G-Eca_G)+gCaT_G.*(xVr_Gt(8,:)+.5*xiR2_G(8,:)).^2.*(xVr_Gt(9,:)+.5*xiR2_G(9,:)).*(vPt_G+.5*rhsV2_G-Eca_G)...
        +gCAN.*((caCon_G+.5*caCr2_G)/(200+(caCon_G+.5*caCr2_G))).*(xVr_Gt(10,:)+.5*xiR2_G(10,:)).*(vPt_G+.5*rhsV2_G-Ecat)+gKCa_G*(xVr_Gt(11,:)+.5*xiR2_G(11,:)).*(vPt_G+.5*rhsV2_G-Ek)...
        +gDR_GP*(xVr_Gt(12,:)+.5*xiR2_G(12,:)).*(vPt_G+.5*rhsV2_G-Ek);
    rhsV3_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV3_G-xi_GC);
    rhsV3_M=1/Rm_GM*(vPt_M+.5*rhsV2_M-El_GM)+gNa_M*(xVr_Mt(1,:)+.5*xiR2_M(1,:)).^3.*(xVr_Mt(2,:)+.5*xiR2_M(2,:)).*(vPt_M+.5*rhsV2_M-Ena)+...
        gNaP*minf(vPt_M+.5*rhsV2_M).*(vPt_M+.5*rhsV2_M-Ena)+gA_M*(xVr_Mt(3,:)+.5*xiR2_M(3,:)).*(xVr_Mt(4,:)+...
        .5*xiR2_M(4,:)).*(vPt_M+.5*rhsV2_M-Ek)+gKS*(xVr_Mt(5,:)+.5*xiR2_M(5,:)).*(xVr_Mt(6,:)+.5*xiR2_M(6,:)).*(vPt_M+.5*rhsV2_M-Ek)+...
        gCaL*(xVr_Mt(7,:)+.5*xiR2_M(7,:)).*(xVr_Mt(8,:)+.5*xiR2_M(8,:)).*(vPt_M+.5*rhsV2_M-Eca_M)+gKCa_M*(xVr_Mt(9,:)+...
        .5*xiR2_M(9,:)).*(vPt_M+.5*rhsV2_M-Ek)+gDR_M*(xVr_Mt(10,:)+.5*xiR2_M(10,:)).^2.*(xVr_Mt(11,:)+.5*xiR2_M(11,:)).*(vPt_M+.5*rhsV2_M-Ek);
    rhsV3_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV3_M-xi_MC);
    rhsV3_P=1/Rm_P*(vPt_P+.5*rhsV2_P-El_P)+gNa_GP*(xVr_Pt(1,:)+.5*xiR2_P(1,:)).^3.*(xVr_Pt(2,:)+.5*xiR2_P(2,:)).*(vPt_P+.5*rhsV2_P-Ena)+...
        gA_P*(xVr_Pt(3,:)+.5*xiR2_P(3,:)).*(xVr_Pt(4,:)+.5*xiR2_P(4,:)).*(vPt_P+.5*rhsV2_P-Ek)+gM_P.*(xVr_Pt(5,:)+.5*xiR2_P(5,:)).*(vPt_P+.5*rhsV2_P-Ek)...
        +gH.*(xVr_Pt(6,:)+.5*xiR2_P(6,:)).*(vPt_P+.5*rhsV2_P-Eh)+gCaPN_P.*(xVr_Pt(7,:)+.5*xiR2_P(7,:)).^2.*(xVr_Pt(8,:)+.5*xiR2_P(8,:)).*(vPt_P+.5*rhsV2_P-Eca_P)...
        +gCaT_P.*(xVr_Pt(9,:)+.5*xiR2_P(9,:)).^2.*(xVr_Pt(10,:)+.5*xiR2_P(10,:)).*(vPt_P+.5*rhsV2_P-Eca_P)...
        +gKCa_P*(xVr_Pt(11,:)+.5*xiR2_P(11,:)).*(vPt_P+.5*rhsV2_P-Ek)+gDR_GP*(xVr_Pt(12,:)+.5*xiR2_P(12,:)).*(vPt_P+.5*rhsV2_P-Ek);
    rhsV3_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV3_P-xi_PGC);
    
    xiR3_G(1,:)=dt*2.1./tauM(vPt_G+.5*rhsV2_G).*(mInf(vPt_G+.5*rhsV2_G)-(xVr_Gt(1,:)+.5*xiR2_G(1,:))); %I_Na
    xiR3_G(2,:)=dt*2.1./tauH(vPt_G+.5*rhsV2_G).*(hInf(vPt_G+.5*rhsV2_G)-(xVr_Gt(2,:)+.5*xiR2_G(2,:))); %I_Na
    xiR3_G(3,:)=dt*3.3./atau(vPt_G+.5*rhsV2_G).*(aminf_g(vPt_G+.5*rhsV2_G)-(xVr_Gt(3,:)+.5*xiR2_G(3,:))); %I_A
    xiR3_G(4,:)=dt*3.3./hatau_g(vPt_G+.5*rhsV2_G).*(ahinf_g(vPt_G+.5*rhsV2_G)-(xVr_Gt(4,:)+.5*xiR2_G(4,:))); %I_A
    xiR3_G(5,:)=dt./tauMusc(vPt_G+.5*rhsV2_G).*(minfMusc(vPt_G+.5*rhsV2_G)-(xVr_Gt(5,:)+.5*xiR2_G(5,:))); %I_M
    xiR3_G(6,:)=dt./capmtau(vPt_G+.5*rhsV2_G).*(capminf(vPt_G+.5*rhsV2_G)-(xVr_Gt(6,:)+.5*xiR2_G(6,:))); %I_CaPN
    xiR3_G(7,:)=dt./caphtau(vPt_G+.5*rhsV2_G).*(caphinf(vPt_G+.5*rhsV2_G)-(xVr_Gt(7,:)+.5*xiR2_G(7,:))); %I_CaPN
    xiR3_G(8,:)=dt./cattaumG(vPt_G+.5*rhsV2_G).*(catminfG(vPt_G+.5*rhsV2_G)-(xVr_Gt(8,:)+.5*xiR2_G(8,:))); %I_CaT
    xiR3_G(9,:)=dt./cattauh(vPt_G+.5*rhsV2_G).*(cathinf(vPt_G+.5*rhsV2_G)-(xVr_Gt(9,:)+.5*xiR2_G(9,:))); %I_CaT
    xiR3_G(10,:)=dt./cantau(vPt_G+.5*rhsV2_G).*(canminf(vPt_G+.5*rhsV2_G)-(xVr_Gt(10,:)+.5*xiR2_G(10,:))); %I_CAN
    xiR3_G(11,:)=dt*(kcaa(vPt_G+.5*rhsV2_G,caCon_G+.5*caCr2_G).*(1-(xVr_Gt(11,:)+.5*xiR2_G(11,:))) - 0.05*(xVr_Gt(11,:)+.5*xiR2_G(11,:))); %I_KCa
    xiR3_G(12,:)=dt*3.3*(mss_dr(vPt_G+.5*rhsV2_G)-(xVr_Gt(12,:)+.5*xiR2_G(12,:)))./taum_dr(vPt_G+.5*rhsV2_G); %I_DR
    
    xiR3_M(1,:)=dt*(ana(vPt_M+.5*rhsV2_M).*(1-(xVr_Mt(1,:)+.5*xiR2_M(1,:))) - bna(vPt_M+.5*rhsV2_M).*(xVr_Mt(1,:)+.5*xiR2_M(1,:)));
    xiR3_M(2,:)=dt*(ahna(vPt_M+.5*rhsV2_M).*(1-(xVr_Mt(2,:)+.5*xiR2_M(2,:))) - bhna(vPt_M+.5*rhsV2_M).*(xVr_Mt(2,:)+.5*xiR2_M(2,:)));
    xiR3_M(3,:)=dt*3.3./atau(vPt_M+.5*rhsV2_M).*(aminf(vPt_M+.5*rhsV2_M)-(xVr_Mt(3,:)+.5*xiR2_M(3,:))); %I_A
    xiR3_M(4,:)=dt*3.3./hatau(vPt_M+.5*rhsV2_M).*(ahinf(vPt_M+.5*rhsV2_M)-(xVr_Mt(4,:)+.5*xiR2_M(4,:))); %I_A
    xiR3_M(5,:)=dt/10*(ksminf(vPt_M+.5*rhsV2_M)-(xVr_Mt(5,:)+.5*xiR2_M(5,:))); %I_KS
    xiR3_M(6,:)=dt./kshtau(vPt_M+.5*rhsV2_M).*(kshinf(vPt_M+.5*rhsV2_M)-(xVr_Mt(6,:)+.5*xiR2_M(6,:)));
    xiR3_M(7,:)=dt*(cala(vPt_M+.5*rhsV2_M).*(1-(xVr_Mt(7,:)+.5*xiR2_M(7,:))) - calb(vPt_M+.5*rhsV2_M).*(xVr_Mt(7,:)+.5*xiR2_M(7,:))); %I_CaL
    xiR3_M(8,:)=dt*(calha(vPt_M+.5*rhsV2_M).*(1-(xVr_Mt(8,:)+.5*xiR2_M(8,:))) - calhb(vPt_M+.5*rhsV2_M).*(xVr_Mt(8,:)+.5*xiR2_M(8,:)));
    xiR3_M(9,:)=dt*(kcaa(vPt_M+.5*rhsV2_M,caCon_M+.5*caCr2_M).*(1-(xVr_Mt(9,:)+.5*xiR2_M(9,:))) - 0.05*(xVr_Mt(9,:)+.5*xiR2_M(9,:)));
    xiR3_M(10,:)=dt*(nss_dr(vPt_M+.5*rhsV2_M)-(xVr_Mt(10,:)+.5*xiR2_M(10,:)))./taun_dr(vPt_M+.5*rhsV2_M);    %I_DR
    xiR3_M(11,:)=dt*(kss_dr(vPt_M+.5*rhsV2_M)-(xVr_Mt(11,:)+.5*xiR2_M(11,:)))./50;
    
    xiR3_P(1,:)=dt*2.1./tauM(vPt_P+.5*rhsV2_P).*(mInf(vPt_P+.5*rhsV2_P)-(xVr_Pt(1,:)+.5*xiR2_P(1,:))); %I_Na
    xiR3_P(2,:)=dt*2.1./tauH(vPt_P+.5*rhsV2_P).*(hInf(vPt_P+.5*rhsV2_P)-(xVr_Pt(2,:)+.5*xiR2_P(2,:))); %I_Na
    xiR3_P(3,:)=dt*3.3./atau(vPt_P+.5*rhsV2_P).*(aminf_g(vPt_P+.5*rhsV2_P)-(xVr_Pt(3,:)+.5*xiR2_P(3,:))); %I_A
    xiR3_P(4,:)=dt*3.3./hatau_g(vPt_P+.5*rhsV2_P).*(ahinf_g(vPt_P+.5*rhsV2_P)-(xVr_Pt(4,:)+.5*xiR2_P(4,:))); %I_A
    xiR3_P(5,:)=dt./tauMusc(vPt_P+.5*rhsV2_P).*(minfMusc(vPt_P+.5*rhsV2_P)-(xVr_Pt(5,:)+.5*xiR2_P(5,:))); %I_M
    xiR3_P(6,:)=dt*2.1./hcurmtau(vPt_P+.5*rhsV2_P).*(hcurminf(vPt_P+.5*rhsV2_P)-(xVr_Pt(6,:))+.5*xiR2_P(6,:)); %I_H
    xiR3_P(7,:)=dt./capmtau(vPt_P+.5*rhsV2_P).*(capminf(vPt_P+.5*rhsV2_P)-(xVr_Pt(7,:)+.5*xiR2_P(7,:))); %I_CaPN
    xiR3_P(8,:)=dt./caphtau(vPt_P+.5*rhsV2_P).*(caphinf(vPt_P+.5*rhsV2_P)-(xVr_Pt(8,:)+.5*xiR2_P(8,:))); %I_CaPN
    xiR3_P(9,:)=dt./cattaumP(vPt_P+.5*rhsV2_P).*(catminfP(vPt_P+.5*rhsV2_P)-(xVr_Pt(9,:)+.5*xiR2_P(9,:))); %I_CaT
    xiR3_P(10,:)=dt./cattauh(vPt_P+.5*rhsV2_P).*(cathinf(vPt_P+.5*rhsV2_P)-(xVr_Pt(10,:)+.5*xiR2_P(10,:))); %I_CaT
    xiR3_P(11,:)=dt*(kcaa(vPt_P+.5*rhsV2_P,caCon_P+.5*caCr2_P).*(1-(xVr_Pt(11,:)+.5*xiR2_P(11,:))) - 0.05*(xVr_Pt(11,:)+.5*xiR2_P(11,:))); %I_KCa
    xiR3_P(12,:)=dt*3.3*(mss_dr(vPt_P+.5*rhsV2_P)-(xVr_Pt(12,:)+.5*xiR2_P(12,:)))./taum_dr(vPt_P+.5*rhsV2_P); %I_DR
    
    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_Gt(6,:)+.5*xiR2_G(6,:)).^2.*(xVr_Gt(7,:)+.5*xiR2_G(7,:)).*(vPt_G+.5*rhsV2_G-Eca_G)+gCaT_G*(xVr_Gt(8,:)+.5*xiR2_G(8,:)).^2.*(xVr_Gt(9,:)+.5*xiR2_G(9,:)).*(vPt_G+.5*rhsV2_G-Eca_G);
    caCr3_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+.5*caCr2_G))/800); %tauCa=800ms for PGC & GC
    
    Ica_M=gCaL*(xVr_Mt(7,:)+.5*xiR2_M(7,:)).*(xVr_Mt(8,:)+.5*xiR2_M(8,:)).*(vPt_M+.5*rhsV2_M-Eca_M);
    caCr3_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr2_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_Pt(7,:)+.5*xiR2_P(7,:)).^2.*(xVr_Pt(8,:)+.5*xiR2_P(8,:)).*(vPt_P+.5*rhsV2_P-Eca_P)+gCaT_P*(xVr_Pt(9,:)+.5*xiR2_P(9,:)).^2.*(xVr_Pt(10,:)+.5*xiR2_P(10,:)).*(vPt_P+.5*rhsV2_P-Eca_P);
    caCr3_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+.5*caCr2_P))/800); %tauCa=800ms for PGC & GC
    
    % --- step 4 ---
    rhsV4_G=1/Rm_GM*(vPt_G+rhsV3_G-El_GM)+gNa_GP*(xVr_Gt(1,:)+xiR3_G(1,:)).^3.*(xVr_Gt(2,:)+xiR3_G(2,:)).*(vPt_G+rhsV3_G-Ena)+...
        gA_G*(xVr_Gt(3,:)+xiR3_G(3,:)).*(xVr_Gt(4,:)+xiR3_G(4,:)).*(vPt_G+rhsV3_G-Ek)+gM_G.*(xVr_Gt(5,:)+xiR3_G(5,:)).*(vPt_G+rhsV3_G-Ek)+...
        gCaPN_G.*(xVr_Gt(6,:)+xiR3_G(6,:)).^2.*(xVr_Gt(7,:)+xiR3_G(7,:)).*(vPt_G+rhsV3_G-Eca_G)+gCaT_G.*(xVr_Gt(8,:)+xiR3_G(8,:)).^2.*(xVr_Gt(9,:)+xiR3_G(9,:)).*(vPt_G+rhsV3_G-Eca_G)...
        +gCAN.*((caCon_G+caCr3_G)/(200+(caCon_G+caCr3_G))).*(xVr_Gt(10,:)+xiR3_G(10,:)).*(vPt_G+rhsV3_G-Ecat)+gKCa_G*(xVr_Gt(11,:)+xiR3_G(11,:)).*(vPt_G+rhsV3_G-Ek)...
        +gDR_GP*(xVr_Gt(12,:)+xiR3_G(12,:)).*(vPt_G+rhsV3_G-Ek);
    rhsV4_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV4_G-xi_GC);
    rhsV4_M=1/Rm_GM*(vPt_M+rhsV3_M-El_GM)+gNa_M*(xVr_Mt(1,:)+xiR3_M(1,:)).^3.*(xVr_Mt(2,:)+xiR3_M(2,:)).*(vPt_M+rhsV3_M-Ena)+...
        gNaP*minf(vPt_M+rhsV3_M).*(vPt_M+rhsV3_M-Ena)+gA_M*(xVr_Mt(3,:)+xiR3_M(3,:)).*(xVr_Mt(4,:)+...
        xiR3_M(4,:)).*(vPt_M+rhsV3_M-Ek)+gKS*(xVr_Mt(5,:)+xiR3_M(5,:)).*(xVr_Mt(6,:)+xiR3_M(6,:)).*(vPt_M+rhsV3_M-Ek)+...
        gCaL*(xVr_Mt(7,:)+xiR3_M(7,:)).*(xVr_Mt(8,:)+xiR3_M(8,:)).*(vPt_M+rhsV3_M-Eca_M)+gKCa_M*(xVr_Mt(9,:)+...
        xiR3_M(9,:)).*(vPt_M+rhsV3_M-Ek)+gDR_M*(xVr_Mt(10,:)+xiR3_M(10,:)).^2.*(xVr_Mt(11,:)+xiR3_M(11,:)).*(vPt_M+rhsV3_M-Ek);
    rhsV4_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV4_M-xi_MC);
    rhsV4_P=1/Rm_P*(vPt_P+rhsV3_P-El_P)+gNa_GP*(xVr_Pt(1,:)+xiR3_P(1,:)).^3.*(xVr_Pt(2,:)+xiR3_P(2,:)).*(vPt_P+rhsV3_P-Ena)+...
        gA_P*(xVr_Pt(3,:)+xiR3_P(3,:)).*(xVr_Pt(4,:)+xiR3_P(4,:)).*(vPt_P+rhsV3_P-Ek)+gM_P.*(xVr_Pt(5,:)+xiR3_P(5,:)).*(vPt_P+rhsV3_P-Ek)...
        +gH.*(xVr_Pt(6,:)+xiR3_P(6,:)).*(vPt_P+rhsV3_P-Eh)+gCaPN_P.*(xVr_Pt(7,:)+xiR3_P(7,:)).^2.*(xVr_Pt(8,:)+xiR3_P(8,:)).*(vPt_P+rhsV3_P-Eca_P)...
        +gCaT_P.*(xVr_Pt(9,:)+xiR3_P(9,:)).^2.*(xVr_Pt(10,:)+xiR3_P(10,:)).*(vPt_P+rhsV3_P-Eca_P)...
        +gKCa_P*(xVr_Pt(11,:)+xiR3_P(11,:)).*(vPt_P+rhsV3_P-Ek)+gDR_GP*(xVr_Pt(12,:)+xiR3_P(12,:)).*(vPt_P+rhsV3_P-Ek);
    rhsV4_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV4_P-xi_PGC);
    
    xiR4_G(1,:)=dt*2.1./tauM(vPt_G+rhsV3_G).*(mInf(vPt_G+rhsV3_G)-(xVr_Gt(1,:)+xiR3_G(1,:))); %I_Na
    xiR4_G(2,:)=dt*2.1./tauH(vPt_G+rhsV3_G).*(hInf(vPt_G+rhsV3_G)-(xVr_Gt(4,:)+xiR3_G(4,:))); %I_Na
    xiR4_G(3,:)=dt*3.3./atau(vPt_G+rhsV3_G).*(aminf_g(vPt_G+rhsV3_G)-(xVr_Gt(3,:)+xiR3_G(3,:))); %I_A
    xiR4_G(4,:)=dt*3.3./hatau_g(vPt_G+rhsV3_G).*(ahinf_g(vPt_G+rhsV3_G)-(xVr_Gt(4,:)+xiR3_G(4,:))); %I_A
    xiR4_G(5,:)=dt./tauMusc(vPt_G+rhsV3_G).*(minfMusc(vPt_G+rhsV3_G)-(xVr_Gt(5,:)+xiR3_G(5,:))); %I_M
    xiR4_G(6,:)=dt./capmtau(vPt_G+rhsV3_G).*(capminf(vPt_G+rhsV3_G)-(xVr_Gt(6,:)+xiR3_G(6,:))); %I_CaPN
    xiR4_G(7,:)=dt./caphtau(vPt_G+rhsV3_G).*(caphinf(vPt_G+rhsV3_G)-(xVr_Gt(7,:)+xiR3_G(7,:))); %I_CaPN
    xiR4_G(8,:)=dt./cattaumG(vPt_G+rhsV3_G).*(catminfG(vPt_G+rhsV3_G)-(xVr_Gt(8,:)+xiR3_G(8,:))); %I_CaT
    xiR4_G(9,:)=dt./cattauh(vPt_G+rhsV3_G).*(cathinf(vPt_G+rhsV3_G)-(xVr_Gt(9,:)+xiR3_G(9,:))); %I_CaT
    xiR4_G(10,:)=dt./cantau(vPt_G+rhsV3_G).*(canminf(vPt_G+rhsV3_G)-(xVr_Gt(10,:)+xiR3_G(10,:))); %I_CAN
    xiR4_G(11,:)=dt*(kcaa(vPt_G+rhsV3_G,caCon_G+caCr3_G).*(1-(xVr_Gt(11,:)+xiR3_G(11,:))) - 0.05*(xVr_Gt(11,:)+xiR3_G(11,:))); %I_KCa
    xiR4_G(12,:)=dt*3.3*(mss_dr(vPt_G+rhsV3_G)-(xVr_Gt(12,:)+xiR3_G(12,:)))./taum_dr(vPt_G+rhsV3_G); %I_DR
    
    xiR4_M(1,:)=dt*(ana(vPt_M+rhsV3_M).*(1-(xVr_Mt(1,:)+xiR3_M(1,:))) - bna(vPt_M+rhsV3_M).*(xVr_Mt(1,:)+xiR3_M(1,:)));
    xiR4_M(2,:)=dt*(ahna(vPt_M+rhsV3_M).*(1-(xVr_Mt(2,:)+xiR3_M(2,:))) - bhna(vPt_M+rhsV3_M).*(xVr_Mt(2,:)+xiR3_M(2,:)));
    xiR4_M(3,:)=dt*3.3./atau(vPt_M+rhsV3_M).*(aminf(vPt_M+rhsV3_M)-(xVr_Mt(3,:)+xiR3_M(3,:))); %I_A
    xiR4_M(4,:)=dt*3.3./hatau(vPt_M+rhsV3_M).*(ahinf(vPt_M+rhsV3_M)-(xVr_Mt(4,:)+xiR3_M(4,:))); %I_A
    xiR4_M(5,:)=dt/10*(ksminf(vPt_M+rhsV3_M)-(xVr_Mt(5,:)+xiR3_M(5,:))); %I_KS
    xiR4_M(6,:)=dt./kshtau(vPt_M+rhsV3_M).*(kshinf(vPt_M+rhsV3_M)-(xVr_Mt(6,:)+xiR3_M(6,:)));
    xiR4_M(7,:)=dt*(cala(vPt_M+rhsV3_M).*(1-(xVr_Mt(7,:)+xiR3_M(7,:))) - calb(vPt_M+rhsV3_M).*(xVr_Mt(7,:)+xiR3_M(7,:))); %I_CaL
    xiR4_M(8,:)=dt*(calha(vPt_M+rhsV3_M).*(1-(xVr_Mt(8,:)+xiR3_M(8,:))) - calhb(vPt_M+rhsV3_M).*(xVr_Mt(8,:)+xiR3_M(8,:)));
    xiR4_M(9,:)=dt*(kcaa(vPt_M+rhsV3_M,caCon_M+caCr3_M).*(1-(xVr_Mt(9,:)+xiR3_M(9,:))) - 0.05*(xVr_Mt(9,:)+xiR3_M(9,:)));
    xiR4_M(10,:)=dt*(nss_dr(vPt_M+rhsV3_M)-(xVr_Mt(10,:)+xiR3_M(10,:)))./taun_dr(vPt_M+rhsV3_M);    %I_DR
    xiR4_M(11,:)=dt*(kss_dr(vPt_M+rhsV3_M)-(xVr_Mt(11,:)+xiR3_M(11,:)))./50;
    
    xiR4_P(1,:)=dt*2.1./tauM(vPt_P+rhsV3_P).*(mInf(vPt_P+rhsV3_P)-(xVr_Pt(1,:)+xiR3_P(1,:))); %I_Na
    xiR4_P(2,:)=dt*2.1./tauH(vPt_P+rhsV3_P).*(hInf(vPt_P+rhsV3_P)-(xVr_Pt(2,:)+xiR3_P(2,:))); %I_Na
    xiR4_P(3,:)=dt*3.3./atau(vPt_P+rhsV3_P).*(aminf_g(vPt_P+rhsV3_P)-(xVr_Pt(3,:)+xiR3_P(3,:))); %I_A
    xiR4_P(4,:)=dt*3.3./hatau_g(vPt_P+rhsV3_P).*(ahinf_g(vPt_P+rhsV3_P)-(xVr_Pt(4,:)+xiR3_P(4,:))); %I_A
    xiR4_P(5,:)=dt./tauMusc(vPt_P+rhsV3_P).*(minfMusc(vPt_P+rhsV3_P)-(xVr_Pt(5,:)+xiR3_P(5,:))); %I_M
    xiR4_P(6,:)=dt*2.1./hcurmtau(vPt_P+rhsV3_P).*(hcurminf(vPt_P+rhsV3_P)-(xVr_Pt(6,:))+xiR3_P(6,:)); %I_H
    xiR4_P(7,:)=dt./capmtau(vPt_P+rhsV3_P).*(capminf(vPt_P+rhsV3_P)-(xVr_Pt(7,:)+xiR3_P(7,:))); %I_CaPN
    xiR4_P(8,:)=dt./caphtau(vPt_P+rhsV3_P).*(caphinf(vPt_P+rhsV3_P)-(xVr_Pt(8,:)+xiR3_P(8,:))); %I_CaPN
    xiR4_P(9,:)=dt./cattaumP(vPt_P+rhsV3_P).*(catminfP(vPt_P+rhsV3_P)-(xVr_Pt(9,:)+xiR3_P(9,:))); %I_CaT
    xiR4_P(10,:)=dt./cattauh(vPt_P+rhsV3_P).*(cathinf(vPt_P+rhsV3_P)-(xVr_Pt(10,:)+xiR3_P(10,:))); %I_CaT
    xiR4_P(11,:)=dt*(kcaa(vPt_P+rhsV3_P,caCon_P+caCr3_P).*(1-(xVr_Pt(11,:)+xiR3_P(11,:))) - 0.05*(xVr_Pt(11,:)+xiR3_P(11,:))); %I_KCa
    xiR4_P(12,:)=dt*3.3*(mss_dr(vPt_P+rhsV3_P)-(xVr_Pt(12,:)+xiR3_P(12,:)))./taum_dr(vPt_P+rhsV3_P); %I_DR
    
    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_Gt(6,:)+xiR3_G(6,:)).^2.*(xVr_Gt(7,:)+xiR3_G(7,:)).*(vPt_G+rhsV3_G-Eca_G)+gCaT_G*(xVr_Gt(8,:)+xiR3_G(8,:)).^2.*(xVr_Gt(9,:)+xiR3_G(9,:)).*(vPt_G+rhsV3_G-Eca_G);
    caCr4_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+caCr3_G))/800); %tauCa=800ms for PGC & GC

    Ica_M=gCaL*(xVr_Mt(7,:)+xiR3_M(7,:)).*(xVr_Mt(8,:)+xiR3_M(8,:)).*(vPt_M+rhsV3_M-Eca_M);
    caCr4_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+caCr3_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_Pt(7,:)+xiR3_P(7,:)).^2.*(xVr_Pt(8,:)+xiR3_P(8,:)).*(vPt_P+rhsV3_P-Eca_P)+gCaT_P*(xVr_Pt(9,:)+xiR3_P(9,:)).^2.*(xVr_Pt(10,:)+xiR3_P(10,:)).*(vPt_P+rhsV3_P-Eca_P);
    caCr4_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+caCr3_P))/800); %tauCa=800ms for PGC & GC
    
    % FINAL Step... update voltage
    vPt_G=vPt_G+1/6*(rhsV1_G+2*rhsV2_G+2*rhsV3_G+rhsV4_G);
    vPt_M=vPt_M+1/6*(rhsV1_M+2*rhsV2_M+2*rhsV3_M+rhsV4_M);
    vPt_P=vPt_P+1/6*(rhsV1_P+2*rhsV2_P+2*rhsV3_P+rhsV4_P);
    
    %update gating variables
    xVr_Gt=xVr_Gt+1/6*(xiR1_G+2*xiR2_G+2*xiR3_G+xiR4_G);
    xVr_Mt=xVr_Mt+1/6*(xiR1_M+2*xiR2_M+2*xiR3_M+xiR4_M);
    xVr_Pt=xVr_Pt+1/6*(xiR1_P+2*xiR2_P+2*xiR3_P+xiR4_P);
    
    %update calcium concentration
    caCon_G=caCon_G+1/6*(caCr1_G+2*caCr2_G+2*caCr3_G+caCr4_G);
    caCon_M=caCon_M+1/6*(caCr1_M+2*caCr2_M+2*caCr3_M+caCr4_M);
    caCon_P=caCon_P+1/6*(caCr1_P+2*caCr2_P+2*caCr3_P+caCr4_P);
    
    %update calcium reversal potential, Nernst eqn
    Eca_G=RTovF./2*log(10./caCon_G); %get 70mV when caCon=0.05 
    Eca_M=RTovF./2*log(10./caCon_M); %get 70mV when caCon=0.05
    Eca_P=RTovF./2*log(10./caCon_P); %get 70mV when caCon=0.05 
    
    %spiking
    idGspk=(vPt_G>vltThres).*(j-tLstSp_G>=mnTim);
    idMspk=(vPt_M>vltThres).*(j-tLstSp_M>=mnTim);
    idPspk=(vPt_P>vltThres).*(j-tLstSp_P>=mnTim);
    %update spiking vars
    if(sum(idGspk)>0) % at least 1 G spikes
        tLstSp_G(idGspk==1)=j;
    end
    if(sum(idMspk)>0)
        tLstSp_M(idMspk==1)=j;
    end
    if(sum(idPspk)>0)
        tLstSp_P(idPspk==1)=j;
    end
end
%shift tLstSp_ by Ltrn
tLstSp_G=tLstSp_G-Ltrn; tLstSp_M=tLstSp_M-Ltrn; tLstSp_P=tLstSp_P-Ltrn;
spTimes_G=[]; spTimes_M=[]; spTimes_P=[]; %make empty; dyn append

% Use SAME volt name, overwrite at each step (vPt_[G/M/P])
vP_G=vPt_G;  vP_M=vPt_M;    vP_P=vPt_P;
xVr_G=xVr_Gt;
xVr_M=xVr_Mt;
xVr_P=xVr_Pt;
%Generating random distribution array for Poisson kicks
randA_G=rand(Lt,3);
randA_Gc=rand(Lt,3);
randA_M=rand(Lt,2);
randA_Mc=rand(Lt,1);
randA_MPc=rand(Lt,2);
randA_PM=rand(Lt,2);
randG_G=rand(Lt,3);
randG_Gc=rand(Lt,3);
randG_M=rand(Lt,2);
randG_Mc=rand(Lt,2);
randG_MPc=rand(Lt,2);
randG_PM=rand(Lt,2);

%Post-transient iterations
for j=2:Lt
    
    %Synaptic currents (GABA/NMDA/AMPA) aprx by Forward-Euler
    s_GABA_G=s_GABA_G+dt*dsGABA(s_GABA_G,vP_G);
    Isyn_GC=[g_GABA_G*(vP_M(1)-Esyn_GP)*(w_G(1)*s_GABA_G(1)+w_Gc*s_GABA_G(2)),...
        g_GABA_G*(vP_M(2)-Esyn_GP)*(w_Gc*s_GABA_G(2)+w_G(2)*s_GABA_G(3))];
    s_AMPA=s_AMPA+dt*dsAMPA(s_AMPA,vP_M);
    s_NMDA=s_NMDA+dt*dsNMDA(s_NMDA,vP_M);
    Isyn_MCG=[w_MGA(1)*g_AMPA*s_AMPA(1)*(vP_G(1)-Esyn_M)+w_MGN(1)*g_NMDA*s_NMDA(1)*B(vP_G(1))*(vP_G(1)-Esyn_M),...
        w_MGA(2)*g_AMPA*(vP_G(2)-Esyn_M)*(s_AMPA(1)+s_AMPA(2))+w_MGN(2)*g_NMDA*B(vP_G(2))*(vP_G(2)-Esyn_M)*(s_NMDA(1)+s_NMDA(2)),...
        w_MGA(3)*g_AMPA*s_AMPA(2)*(vP_G(3)-Esyn_M)+w_MGN(3)*g_NMDA*s_NMDA(2)*B(vP_G(3))*(vP_G(3)-Esyn_M)];
    Isyn_MCP=[w_MPA(1)*g_AMPA*s_AMPA(1)*(vP_P(1)-Esyn_M)+w_MPN(1)*g_NMDA*s_NMDA(1)*B(vP_P(1))*(vP_P(1)-Esyn_M),...
        w_MPA(2)*g_AMPA*s_AMPA(2)*(vP_P(2)-Esyn_M)+w_MPN(2)*g_NMDA*s_NMDA(2)*B(vP_P(2))*(vP_P(2)-Esyn_M)];
    s_GABA_P=s_GABA_P+dt*dsGABA(s_GABA_P,vP_P);
    Isyn_PGC=[w_P(1)*g_GABA_P*s_GABA_P(1).*(vP_M(1)-Esyn_GP),...
        w_P(2)*g_GABA_P*s_GABA_P(2).*(vP_M(2)-Esyn_GP)];
    
     %Background noise Poisson kicks with AMPA synapse
    sA_GC=sA_GC+(dt/tau_sGP)*-sA_GC;
    sA_GC=[sA_GC(1)+a(1)*(randA_G(j,1)<nuaugA_G12(j)*dt)+a(1)*(randA_G(j,1)<nuaugA_G13(j)*dt)+a(1)*(randA_Gc(j,1)<numinA_G12(j)*dt)+a(1)*(randA_Gc(j,3)<numinA_G13(j)*dt),...
        sA_GC(2)+a(2)*(randA_G(j,2)<nuaugA_G21(j)*dt)+a(2)*(randA_G(j,2)<nuaugA_G23(j)*dt)+a(2)*(randA_Gc(j,1)<numinA_G12(j)*dt)+1*(randA_Gc(j,2)<numinA_G23(j)*dt),...
        sA_GC(3)+a(3)*(randA_G(j,3)<nuaugA_G31(j)*dt)+a(3)*(randA_G(j,3)<nuaugA_G32(j)*dt)+a(3)*(randA_Gc(j,2)<numinA_G23(j)*dt)+1*(randA_Gc(j,3)<numinA_G13(j)*dt)];
    xiA_GC=sA_GC.*(vP_G-Esyn_M);
    sA_MC=sA_MC+(dt/tau_sM)*-sA_MC;
    sA_MC=[sA_MC(1)+a(4)*(randA_M(j,1)<nuaugA_M1(j)*dt)+a(4)*(randA_Mc(j)<numinA_M(j)*dt)+a(4)*(randA_M(j,1)<nuaugA_MP11(j)*dt)+a(4)*(randA_M(j,1)<nuaugA_MP12(j)*dt)+a(4)*(randA_MPc(j,1)<numinA_MP11(j)*dt)+a(4)*(randA_MPc(j,1)<numinA_MP12(j)*dt),...
        sA_MC(2)+a(5)*(randA_M(j,2)<nuaugA_M2(j)*dt)+a(5)*(randA_Mc(j)<numinA_M(j)*dt)+a(5)*(randA_M(j,2)<nuaugA_MP21(j)*dt)+a(5)*(randA_M(j,2)<nuaugA_MP22(j)*dt)+a(5)*(randA_MPc(j,2)<numinA_MP21(j)*dt)+a(5)*(randA_MPc(j,2)<numinA_MP22(j)*dt)];
    xiA_MC=sA_MC.*(vP_M-Esyn_M);
    sA_PGC=sA_PGC+(dt/tau_sGP)*-sA_PGC;
    sA_PGC=[sA_PGC(1)+a(6)*(randA_PM(j,1)<nuaugA_PM11(j)*dt)+a(6)*(randA_PM(j,1)<nuaugA_PM21(j)*dt)+a(6)*(randA_MPc(j,1)<numinA_MP11(j)*dt)+a(6)*(randA_MPc(j,2)<numinA_MP21(j)*dt),...
        sA_PGC(2)+a(7)*(randA_PM(j,2)<nuaugA_PM12(j)*dt)+a(7)*(randA_PM(j,2)<nuaugA_PM22(j)*dt)+a(7)*(randA_MPc(j,1)<numinA_MP12(j)*dt)+a(7)*(randA_MPc(j,2)<numinA_MP22(j)*dt)];
    xiA_PGC=sA_PGC.*(vP_P-Esyn_M);
    %Background noise Poisson kicks with GABA synapse
    sG_GC=sG_GC+(dt/tau_sGP)*-sG_GC;
    sG_GC=[sG_GC(1)+a(1)*(randG_G(j,1)<nuaugG_G12(j)*dt)+a(1)*(randG_G(j,1)<nuaugG_G13(j)*dt)+a(1)*(randG_Gc(j,1)<numinG_G12(j)*dt)+a(1)*(randG_Gc(j,3)<numinG_G13(j)*dt),...
        sG_GC(2)+a(2)*(randG_G(j,2)<nuaugG_G21(j)*dt)+a(2)*(randG_G(j,2)<nuaugG_G23(j)*dt)+a(2)*(randG_Gc(j,1)<numinG_G12(j)*dt)+1*(randG_Gc(j,2)<numinG_G23(j)*dt),...
        sG_GC(3)+a(3)*(randG_G(j,3)<nuaugG_G31(j)*dt)+a(3)*(randG_G(j,3)<nuaugG_G32(j)*dt)+a(3)*(randG_Gc(j,2)<numinG_G23(j)*dt)+1*(randG_Gc(j,3)<numinG_G13(j)*dt)];
    xiG_GC=sG_GC.*(vP_G-Esyn_GP);
    sG_MC=sG_MC+(dt/tau_sM)*-sG_MC;
    sG_MC=[sG_MC(1)+a(4)*(randG_M(j,1)<nuaugG_M1(j)*dt)+a(4)*(randG_Mc(j)<numinG_M(j)*dt)+a(4)*(randG_M(j,1)<nuaugG_MP11(j)*dt)+a(4)*(randG_M(j,1)<nuaugG_MP12(j)*dt)+a(4)*(randG_MPc(j,1)<numinG_MP11(j)*dt)+a(4)*(randG_MPc(j,1)<numinG_MP12(j)*dt),...
        sG_MC(2)+a(5)*(randG_M(j,2)<nuaugG_M2(j)*dt)+a(5)*(randG_Mc(j)<numinG_M(j)*dt)+a(5)*(randG_M(j,2)<nuaugG_MP21(j)*dt)+a(5)*(randG_M(j,2)<nuaugG_MP22(j)*dt)+a(5)*(randG_MPc(j,2)<numinG_MP21(j)*dt)+a(5)*(randG_MPc(j,2)<numinG_MP22(j)*dt)];
    xiG_MC=sG_MC.*(vP_M-Esyn_GP);
    sG_PGC=sG_PGC+(dt/tau_sGP)*-sG_PGC;
    sG_PGC=[sG_PGC(1)+a(6)*(randG_PM(j,1)<nuaugG_PM11(j)*dt)+a(6)*(randG_PM(j,1)<nuaugG_PM21(j)*dt)+a(6)*(randG_MPc(j,1)<numinG_MP11(j)*dt)+a(6)*(randG_MPc(j,2)<numinG_MP21(j)*dt),...
        sG_PGC(2)+a(7)*(randG_PM(j,2)<nuaugG_PM12(j)*dt)+a(7)*(randG_PM(j,2)<nuaugG_PM22(j)*dt)+a(7)*(randG_MPc(j,1)<numinG_MP12(j)*dt)+a(7)*(randG_MPc(j,2)<numinG_MP22(j)*dt)];
    xiG_PGC=sG_PGC.*(vP_P-Esyn_GP);
    %Sum background noise
    xi_GC=xiA_GC+xiG_GC;
    xi_MC=xiA_MC+xiG_MC;
    xi_PGC=xiA_PGC+xiG_PGC;

    % --- step 1 ---
    rhsV1_G=1/Rm_GM*(vP_G-El_GM)+gNa_GP*xVr_G(1,:).^3.*xVr_G(2,:).*(vP_G-Ena)+...
        gA_G*xVr_G(3,:).*xVr_G(4,:).*(vP_G-Ek)+gM_G.*xVr_G(5,:).*(vP_G-Ek)+...
        gCaPN_G.*xVr_G(6,:).^2.*xVr_G(7,:).*(vP_G-Eca_G)+gCaT_G.*xVr_G(8,:).^2.*xVr_G(9,:).*(vP_G-Eca_G)...
        +gCAN*(caCon_G/(200+caCon_G)).*xVr_G(10,:).*(vP_G-Ecat)+gKCa_G*xVr_G(11,:).*(vP_G-Ek)...
        +gDR_GP*xVr_G(12,:).*(vP_G-Ek);
    rhsV1_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV1_G-xi_GC);
    rhsV1_M=1/Rm_GM*(vP_M-El_GM)+gNa_M*xVr_M(1,:).^3.*xVr_M(2,:).*(vP_M-Ena)+gNaP*minf(vP_M).*(vP_M-Ena)+...
        gA_M*xVr_M(3,:).*xVr_M(4,:).*(vP_M-Ek)+gKS*xVr_M(5,:).*xVr_M(6,:).*(vP_M-Ek)+...
        gCaL*xVr_M(7,:).*xVr_M(8,:).*(vP_M-Eca_M)+gKCa_M*xVr_M(9,:).*(vP_M-Ek)+gDR_M*xVr_M(10,:).^2.*xVr_M(11,:).*(vP_M-Ek);
    rhsV1_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV1_M-xi_MC);
    rhsV1_P=1/Rm_P*(vP_P-El_P)+gNa_GP*xVr_P(1,:).^3.*xVr_P(2,:).*(vP_P-Ena)+...
        gA_P*xVr_P(3,:).*xVr_P(4,:).*(vP_P-Ek)+gM_P.*xVr_P(5,:).*(vP_P-Ek)+gH.*xVr_P(6,:).*(vP_P-Eh)+...
        gCaPN_P.*xVr_P(7,:).^2.*xVr_P(8,:).*(vP_P-Eca_P)+gCaT_P.*xVr_P(9,:).^2.*xVr_P(10,:).*(vP_P-Eca_P)...
        +gKCa_P*xVr_P(11,:).*(vP_P-Ek)+gDR_GP*xVr_P(12,:).*(vP_P-Ek);
    rhsV1_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV1_P-xi_PGC);
    
    xiR1_G(1,:)=dt*2.1./tauM(vP_G).*(mInf(vP_G)-xVr_G(1,:)); %I_Na
    xiR1_G(2,:)=dt*2.1./tauH(vP_G).*(hInf(vP_G)-xVr_G(2,:)); %I_Na
    xiR1_G(3,:)=dt*3.3./atau(vP_G).*(aminf_g(vP_G)-xVr_G(3,:)); %I_A
    xiR1_G(4,:)=dt*3.3./hatau_g(vP_G).*(ahinf_g(vP_G)-xVr_G(4,:)); %I_A
    xiR1_G(5,:)=dt./tauMusc(vP_G).*(minfMusc(vP_G)-xVr_G(5,:)); %I_M
    xiR1_G(6,:)=dt./capmtau(vP_G).*(capminf(vP_G)-xVr_G(6,:)); %I_CaPN
    xiR1_G(7,:)=dt./caphtau(vP_G).*(caphinf(vP_G)-xVr_G(7,:)); %I_CaPN
    xiR1_G(8,:)=dt./cattaumG(vP_G).*(catminfG(vP_G)-xVr_G(8,:)); %I_CaT
    xiR1_G(9,:)=dt./cattauh(vP_G).*(cathinf(vP_G)-xVr_G(9,:)); %I_CaT
    xiR1_G(10,:)=dt./cantau(vP_G).*(canminf(vP_G)-xVr_G(10,:)); %I_CAN
    xiR1_G(11,:)=dt*(kcaa(vP_G,caCon_G).*(1-xVr_G(11,:)) - 0.05*xVr_G(11,:)); %I_KCa
    xiR1_G(12,:)=dt*3.3*(mss_dr(vP_G)-xVr_G(12,:))./taum_dr(vP_G); %I_DR
    
    xiR1_M(1,:)=dt*(ana(vP_M).*(1-xVr_M(1,:)) - bna(vP_M).*xVr_M(1,:)); %I_Na
    xiR1_M(2,:)=dt*(ahna(vP_M).*(1-xVr_M(2,:)) - bhna(vP_M).*xVr_M(2,:)); %I_Na
    xiR1_M(3,:)=dt*3.3./atau(vP_M).*(aminf(vP_M)-xVr_M(3,:)); %I_A
    xiR1_M(4,:)=dt*3.3./hatau(vP_M).*(ahinf(vP_M)-xVr_M(4,:)); %I_A
    xiR1_M(5,:)=dt/10*(ksminf(vP_M)-xVr_M(5,:)); %I_KS
    xiR1_M(6,:)=dt./kshtau(vP_M).*(kshinf(vP_M)-xVr_M(6,:));
    xiR1_M(7,:)=dt*(cala(vP_M).*(1-xVr_M(7,:)) - calb(vP_M).*xVr_M(7,:)); %I_CaL
    xiR1_M(8,:)=dt*(calha(vP_M).*(1-xVr_M(8,:)) - calhb(vP_M).*xVr_M(8,:));
    xiR1_M(9,:)=dt*(kcaa(vP_M,caCon_M).*(1-xVr_M(9,:)) - 0.05*xVr_M(9,:));
    xiR1_M(10,:)=dt*(nss_dr(vP_M)-xVr_M(10,:))./taun_dr(vP_M);    %I_DR
    xiR1_M(11,:)=dt*(kss_dr(vP_M)-xVr_M(11,:))./50;
    
    xiR1_P(1,:)=dt*2.1./tauM(vP_P).*(mInf(vP_P)-xVr_P(1,:)); %I_Na
    xiR1_P(2,:)=dt*2.1./tauH(vP_P).*(hInf(vP_P)-xVr_P(2,:)); %I_Na
    xiR1_P(3,:)=dt*3.3./atau(vP_P).*(aminf_g(vP_P)-xVr_P(3,:)); %I_A
    xiR1_P(4,:)=dt*3.3./hatau_g(vP_P).*(ahinf_g(vP_P)-xVr_P(4,:)); %I_A
    xiR1_P(5,:)=dt./tauMusc(vP_P).*(minfMusc(vP_P)-xVr_P(5,:)); %I_M
    xiR1_P(6,:)=dt*2.1./hcurmtau(vP_P).*(hcurminf(vP_P)-xVr_P(6,:)); %I_H
    xiR1_P(7,:)=dt./capmtau(vP_P).*(capminf(vP_P)-xVr_P(7,:)); %I_CaPN
    xiR1_P(8,:)=dt./caphtau(vP_P).*(caphinf(vP_P)-xVr_P(8,:)); %I_CaPN
    xiR1_P(9,:)=dt./cattaumP(vP_P).*(catminfP(vP_P)-xVr_P(9,:)); %I_CaT
    xiR1_P(10,:)=dt./cattauh(vP_P).*(cathinf(vP_P)-xVr_P(10,:)); %I_CaT
    xiR1_P(11,:)=dt*(kcaa(vP_P,caCon_P).*(1-xVr_P(11,:)) - 0.05*xVr_P(11,:)); %I_KCa
    xiR1_P(12,:)=dt*3.3*(mss_dr(vP_P)-xVr_P(12,:))./taum_dr(vP_P); %I_DR

    %update calcium concentration
    Ica_G=gCaPN_G*xVr_G(6,:).^2.*xVr_G(7,:).*(vP_G-Eca_G)+gCaT_G*xVr_G(8,:).^2.*xVr_G(9,:).*(vP_G-Eca_G);
    caCr1_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-caCon_G)/800); %tauCa=800ms for PGC & GC
    
    Ica_M=gCaL*xVr_M(7,:).*xVr_M(8,:).*(vP_M-Eca_M);
    caCr1_M=dt*(-Ica_M/(2*F*thick_M) + (.05-caCon_M)/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*xVr_P(7,:).^2.*xVr_P(8,:).*(vP_P-Eca_P)+gCaT_P*xVr_P(9,:).^2.*xVr_P(10,:).*(vP_P-Eca_P);
    caCr1_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-caCon_P)/800); %tauCa=800ms for PGC & GC
    
    % --- step 2 ---
    rhsV2_G=1/Rm_GM*(vP_G+.5*rhsV1_G-El_GM)+gNa_GP*(xVr_G(1,:)+.5*xiR1_G(1,:)).^3.*(xVr_G(2,:)+.5*xiR1_G(2,:)).*(vP_G+.5*rhsV1_G-Ena)+...
        gA_G*(xVr_G(3,:)+.5*xiR1_G(3,:)).*(xVr_G(4,:)+.5*xiR1_G(4,:)).*(vP_G+.5*rhsV1_G-Ek)+gM_G.*(xVr_G(5,:)+.5*xiR1_G(5,:)).*(vP_G+.5*rhsV1_G-Ek)+...
        gCaPN_G.*(xVr_G(6,:)+.5*xiR1_G(6,:)).^2.*(xVr_G(7,:)+.5*xiR1_G(7,:)).*(vP_G+.5*rhsV1_G-Eca_G)+gCaT_G.*(xVr_G(8,:)+.5*xiR1_G(8,:)).^2.*(xVr_G(9,:)+.5*xiR1_G(9,:)).*(vP_G+.5*rhsV1_G-Eca_G)...
        +gCAN.*((caCon_G+.5*caCr1_G)/(200+(caCon_G+.5*caCr1_G))).*(xVr_G(10,:)+.5*xiR1_G(10,:)).*(vP_G+.5*rhsV1_G-Ecat)+gKCa_G*(xVr_G(11,:)+.5*xiR1_G(11,:)).*(vP_G+.5*rhsV1_G-Ek)...
        +gDR_GP*(xVr_G(12,:)+.5*xiR1_G(12,:)).*(vP_G+.5*rhsV1_G-Ek);
    rhsV2_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV2_G-xi_GC);
    rhsV2_M=1/Rm_GM*(vP_M+.5*rhsV1_M-El_GM)+gNa_M*(xVr_M(1,:)+.5*xiR1_M(1,:)).^3.*(xVr_M(2,:)+.5*xiR1_M(2,:)).*(vP_M+.5*rhsV1_M-Ena)+...
        gNaP*minf(vP_M+.5*rhsV1_M).*(vP_M+.5*rhsV1_M-Ena)+gA_M*(xVr_M(3,:)+.5*xiR1_M(3,:)).*(xVr_M(4,:)+...
        .5*xiR1_M(4,:)).*(vP_M+.5*rhsV1_M-Ek)+gKS*(xVr_M(5,:)+.5*xiR1_M(5,:)).*(xVr_M(6,:)+.5*xiR1_M(6,:)).*(vP_M+.5*rhsV1_M-Ek)+...
        gCaL*(xVr_M(7,:)+.5*xiR1_M(7,:)).*(xVr_M(8,:)+.5*xiR1_M(8,:)).*(vP_M+.5*rhsV1_M-Eca_M)+gKCa_M*(xVr_M(9,:)+...
        .5*xiR1_M(9,:)).*(vP_M+.5*rhsV1_M-Ek)+gDR_M*(xVr_M(10,:)+.5*xiR1_M(10,:)).^2.*(xVr_M(11,:)+.5*xiR1_M(11,:)).*(vP_M+.5*rhsV1_M-Ek);
    rhsV2_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV2_M-xi_MC);
    rhsV2_P=1/Rm_P*(vP_P+.5*rhsV1_P-El_P)+gNa_GP*(xVr_P(1,:)+.5*xiR1_P(1,:)).^3.*(xVr_P(2,:)+.5*xiR1_P(2,:)).*(vP_P+.5*rhsV1_P-Ena)+...
        gA_P*(xVr_P(3,:)+.5*xiR1_P(3,:)).*(xVr_P(4,:)+.5*xiR1_P(4,:)).*(vP_P+.5*rhsV1_P-Ek)+gM_P.*(xVr_P(5,:)+.5*xiR1_P(5,:)).*(vP_P+.5*rhsV1_P-Ek)...
        +gH.*(xVr_P(6,:)+.5*xiR1_P(6,:)).*(vP_P+.5*rhsV1_P-Eh)+gCaPN_P.*(xVr_P(7,:)+.5*xiR1_P(7,:)).^2.*(xVr_P(8,:)+.5*xiR1_P(8,:)).*(vP_P+.5*rhsV1_P-Eca_P)...
        +gCaT_P.*(xVr_P(9,:)+.5*xiR1_P(9,:)).^2.*(xVr_P(10,:)+.5*xiR1_P(10,:)).*(vP_P+.5*rhsV1_P-Eca_P)...
        +gKCa_P*(xVr_P(11,:)+.5*xiR1_P(11,:)).*(vP_P+.5*rhsV1_P-Ek)+gDR_GP*(xVr_P(12,:)+.5*xiR1_P(12,:)).*(vP_P+.5*rhsV1_P-Ek);
    rhsV2_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV2_P-xi_PGC);
    
    xiR2_G(1,:)=dt*2.1./tauM(vP_G+.5*rhsV1_G).*(mInf(vP_G+.5*rhsV1_G)-(xVr_G(1,:)+.5*xiR1_G(1,:))); %I_Na
    xiR2_G(2,:)=dt*2.1./tauH(vP_G+.5*rhsV1_G).*(hInf(vP_G+.5*rhsV1_G)-(xVr_G(2,:)+.5*xiR1_G(2,:))); %I_Na
    xiR2_G(3,:)=dt*3.3./atau(vP_G+.5*rhsV1_G).*(aminf_g(vP_G+.5*rhsV1_G)-(xVr_G(3,:)+.5*xiR1_G(3,:))); %I_A
    xiR2_G(4,:)=dt*3.3./hatau_g(vP_G+.5*rhsV1_G).*(ahinf_g(vP_G+.5*rhsV1_G)-(xVr_G(4,:)+.5*xiR1_G(4,:))); %I_A
    xiR2_G(5,:)=dt./tauMusc(vP_G+.5*rhsV1_G).*(minfMusc(vP_G+.5*rhsV1_G)-(xVr_G(5,:)+.5*xiR1_G(5,:))); %I_M
    xiR2_G(6,:)=dt./capmtau(vP_G+.5*rhsV1_G).*(capminf(vP_G+.5*rhsV1_G)-(xVr_G(6,:)+.5*xiR1_G(6,:))); %I_CaPN
    xiR2_G(7,:)=dt./caphtau(vP_G+.5*rhsV1_G).*(caphinf(vP_G+.5*rhsV1_G)-(xVr_G(7,:)+.5*xiR1_G(7,:))); %I_CaPN
    xiR2_G(8,:)=dt./cattaumG(vP_G+.5*rhsV1_G).*(catminfG(vP_G+.5*rhsV1_G)-(xVr_G(8,:)+.5*xiR1_G(8,:))); %I_CaT
    xiR2_G(9,:)=dt./cattauh(vP_G+.5*rhsV1_G).*(cathinf(vP_G+.5*rhsV1_G)-(xVr_G(9,:)+.5*xiR1_G(9,:))); %I_CaT
    xiR2_G(10,:)=dt./cantau(vP_G+.5*rhsV1_G).*(canminf(vP_G+.5*rhsV1_G)-(xVr_G(10,:)+.5*xiR1_G(10,:))); %I_CAN
    xiR2_G(11,:)=dt*(kcaa(vP_G+.5*rhsV1_G,caCon_G+.5*caCr1_G).*(1-(xVr_G(11,:)+.5*xiR1_G(11,:))) - 0.05*(xVr_G(11,:)+.5*xiR1_G(11,:))); %I_KCa
    xiR2_G(12,:)=dt*3.3*(mss_dr(vP_G+.5*rhsV1_G)-(xVr_G(12,:)+.5*xiR1_G(12,:)))./taum_dr(vP_G+.5*rhsV1_G); %I_DR
    
    xiR2_M(1,:)=dt*(ana(vP_M+.5*rhsV1_M).*(1-(xVr_M(1,:)+.5*xiR1_M(1,:))) - bna(vP_M+.5*rhsV1_M).*(xVr_M(1,:)+.5*xiR1_M(1,:)));
    xiR2_M(2,:)=dt*(ahna(vP_M+.5*rhsV1_M).*(1-(xVr_M(2,:)+.5*xiR1_M(2,:))) - bhna(vP_M+.5*rhsV1_M).*(xVr_M(2,:)+.5*xiR1_M(2,:)));
    xiR2_M(3,:)=dt*3.3./atau(vP_M+.5*rhsV1_M).*(aminf(vP_M+.5*rhsV1_M)-(xVr_M(3,:)+.5*xiR1_M(3,:))); %I_A
    xiR2_M(4,:)=dt*3.3./hatau(vP_M+.5*rhsV1_M).*(ahinf(vP_M+.5*rhsV1_M)-(xVr_M(4,:)+.5*xiR1_M(4,:))); %I_A
    xiR2_M(5,:)=dt/10*(ksminf(vP_M+.5*rhsV1_M)-(xVr_M(5,:)+.5*xiR1_M(5,:))); %I_KS
    xiR2_M(6,:)=dt./kshtau(vP_M+.5*rhsV1_M).*(kshinf(vP_M+.5*rhsV1_M)-(xVr_M(6,:)+.5*xiR1_M(6,:)));
    xiR2_M(7,:)=dt*(cala(vP_M+.5*rhsV1_M).*(1-(xVr_M(7,:)+.5*xiR1_M(7,:))) - calb(vP_M+.5*rhsV1_M).*(xVr_M(7,:)+.5*xiR1_M(7,:))); %I_CaL
    xiR2_M(8,:)=dt*(calha(vP_M+.5*rhsV1_M).*(1-(xVr_M(8,:)+.5*xiR1_M(8,:))) - calhb(vP_M+.5*rhsV1_M).*(xVr_M(8,:)+.5*xiR1_M(8,:)));
    xiR2_M(9,:)=dt*(kcaa(vP_M+.5*rhsV1_M,caCon_M+.5*caCr1_M).*(1-(xVr_M(9,:)+.5*xiR1_M(9,:))) - 0.05*(xVr_M(9,:)+.5*xiR1_M(9,:)));
    xiR2_M(10,:)=dt*(nss_dr(vP_M+.5*rhsV1_M)-(xVr_M(10,:)+.5*xiR1_M(10,:)))./taun_dr(vP_M+.5*rhsV1_M);    %I_DR
    xiR2_M(11,:)=dt*(kss_dr(vP_M+.5*rhsV1_M)-(xVr_M(11,:)+.5*xiR1_M(11,:)))./50;
    
    xiR2_P(1,:)=dt*2.1./tauM(vP_P+.5*rhsV1_P).*(mInf(vP_P+.5*rhsV1_P)-(xVr_P(1,:)+.5*xiR1_P(1,:))); %I_Na
    xiR2_P(2,:)=dt*2.1./tauH(vP_P+.5*rhsV1_P).*(hInf(vP_P+.5*rhsV1_P)-(xVr_P(2,:)+.5*xiR1_P(2,:))); %I_Na
    xiR2_P(3,:)=dt*3.3./atau(vP_P+.5*rhsV1_P).*(aminf_g(vP_P+.5*rhsV1_P)-(xVr_P(3,:)+.5*xiR1_P(3,:))); %I_A
    xiR2_P(4,:)=dt*3.3./hatau_g(vP_P+.5*rhsV1_P).*(ahinf_g(vP_P+.5*rhsV1_P)-(xVr_P(4,:)+.5*xiR1_P(4,:))); %I_A
    xiR2_P(5,:)=dt./tauMusc(vP_P+.5*rhsV1_P).*(minfMusc(vP_P+.5*rhsV1_P)-(xVr_P(5,:)+.5*xiR1_P(5,:))); %I_M
    xiR2_P(6,:)=dt*2.1./hcurmtau(vP_P+.5*rhsV1_P).*(hcurminf(vP_P+.5*rhsV1_P)-(xVr_P(6,:))+.5*xiR1_P(6,:)); %I_H
    xiR2_P(7,:)=dt./capmtau(vP_P+.5*rhsV1_P).*(capminf(vP_P+.5*rhsV1_P)-(xVr_P(7,:)+.5*xiR1_P(7,:))); %I_CaPN
    xiR2_P(8,:)=dt./caphtau(vP_P+.5*rhsV1_P).*(caphinf(vP_P+.5*rhsV1_P)-(xVr_P(8,:)+.5*xiR1_P(8,:))); %I_CaPN
    xiR2_P(9,:)=dt./cattaumP(vP_P+.5*rhsV1_P).*(catminfP(vP_P+.5*rhsV1_P)-(xVr_P(9,:)+.5*xiR1_P(9,:))); %I_CaT
    xiR2_P(10,:)=dt./cattauh(vP_P+.5*rhsV1_P).*(cathinf(vP_P+.5*rhsV1_P)-(xVr_P(10,:)+.5*xiR1_P(10,:))); %I_CaT
    xiR2_P(11,:)=dt*(kcaa(vP_P+.5*rhsV1_P,caCon_P+.5*caCr1_P).*(1-(xVr_P(11,:)+.5*xiR1_P(11,:))) - 0.05*(xVr_P(11,:)+.5*xiR1_P(11,:))); %I_KCa
    xiR2_P(12,:)=dt*3.3*(mss_dr(vP_P+.5*rhsV1_P)-(xVr_P(12,:)+.5*xiR1_P(12,:)))./taum_dr(vP_P+.5*rhsV1_P); %I_DR

    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_G(6,:)+.5*xiR1_G(6,:)).^2.*(xVr_G(7,:)+.5*xiR1_G(7,:)).*(vP_G+.5*rhsV1_G-Eca_G)+gCaT_G*(xVr_G(8,:)+.5*xiR1_G(8,:)).^2.*(xVr_G(9,:)+.5*xiR1_G(9,:)).*(vP_G+.5*rhsV1_G-Eca_G);
    caCr2_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+.5*caCr1_G))/800); %tauCa=800ms for PGC & GC

    Ica_M=gCaL*(xVr_M(7,:)+.5*xiR1_M(7,:)).*(xVr_M(8,:)+.5*xiR1_M(8,:)).*(vP_M+.5*rhsV1_M-Eca_M);
    caCr2_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr1_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_P(7,:)+.5*xiR1_P(7,:)).^2.*(xVr_P(8,:)+.5*xiR1_P(8,:)).*(vP_P+.5*rhsV1_P-Eca_P)+gCaT_P*(xVr_P(9,:)+.5*xiR1_P(9,:)).^2.*(xVr_P(10,:)+.5*xiR1_P(10,:)).*(vP_P+.5*rhsV1_P-Eca_P);
    caCr2_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+.5*caCr1_P))/800); %tauCa=800ms for PGC & GC
    
    % --- step 3 ---
    rhsV3_G=1/Rm_GM*(vP_G+.5*rhsV2_G-El_GM)+gNa_GP*(xVr_G(1,:)+.5*xiR2_G(1,:)).^3.*(xVr_G(2,:)+.5*xiR2_G(2,:)).*(vP_G+.5*rhsV2_G-Ena)+...
        gA_G*(xVr_G(3,:)+.5*xiR2_G(3,:)).*(xVr_G(4,:)+.5*xiR2_G(4,:)).*(vP_G+.5*rhsV2_G-Ek)+gM_G.*(xVr_G(5,:)+.5*xiR2_G(5,:)).*(vP_G+.5*rhsV2_G-Ek)+...
        gCaPN_G.*(xVr_G(6,:)+.5*xiR2_G(6,:)).^2.*(xVr_G(7,:)+.5*xiR2_G(7,:)).*(vP_G+.5*rhsV2_G-Eca_G)+gCaT_G.*(xVr_G(8,:)+.5*xiR2_G(8,:)).^2.*(xVr_G(9,:)+.5*xiR2_G(9,:)).*(vP_G+.5*rhsV2_G-Eca_G)...
        +gCAN.*((caCon_G+.5*caCr2_G)/(200+(caCon_G+.5*caCr2_G))).*(xVr_G(10,:)+.5*xiR2_G(10,:)).*(vP_G+.5*rhsV2_G-Ecat)+gKCa_G*(xVr_G(11,:)+.5*xiR2_G(11,:)).*(vP_G+.5*rhsV2_G-Ek)...
        +gDR_GP*(xVr_G(12,:)+.5*xiR2_G(12,:)).*(vP_G+.5*rhsV2_G-Ek);
    rhsV3_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV3_G-xi_GC);
    rhsV3_M=1/Rm_GM*(vP_M+.5*rhsV2_M-El_GM)+gNa_M*(xVr_M(1,:)+.5*xiR2_M(1,:)).^3.*(xVr_M(2,:)+.5*xiR2_M(2,:)).*(vP_M+.5*rhsV2_M-Ena)+...
        gNaP*minf(vP_M+.5*rhsV2_M).*(vP_M+.5*rhsV2_M-Ena)+gA_M*(xVr_M(3,:)+.5*xiR2_M(3,:)).*(xVr_M(4,:)+...
        .5*xiR2_M(4,:)).*(vP_M+.5*rhsV2_M-Ek)+gKS*(xVr_M(5,:)+.5*xiR2_M(5,:)).*(xVr_M(6,:)+.5*xiR2_M(6,:)).*(vP_M+.5*rhsV2_M-Ek)+...
        gCaL*(xVr_M(7,:)+.5*xiR2_M(7,:)).*(xVr_M(8,:)+.5*xiR2_M(8,:)).*(vP_M+.5*rhsV2_M-Eca_M)+gKCa_M*(xVr_M(9,:)+...
        .5*xiR2_M(9,:)).*(vP_M+.5*rhsV2_M-Ek)+gDR_M*(xVr_M(10,:)+.5*xiR2_M(10,:)).^2.*(xVr_M(11,:)+.5*xiR2_M(11,:)).*(vP_M+.5*rhsV2_M-Ek);
    rhsV3_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV3_M-xi_MC);
    rhsV3_P=1/Rm_P*(vP_P+.5*rhsV2_P-El_P)+gNa_GP*(xVr_P(1,:)+.5*xiR2_P(1,:)).^3.*(xVr_P(2,:)+.5*xiR2_P(2,:)).*(vP_P+.5*rhsV2_P-Ena)+...
        gA_P*(xVr_P(3,:)+.5*xiR2_P(3,:)).*(xVr_P(4,:)+.5*xiR2_P(4,:)).*(vP_P+.5*rhsV2_P-Ek)+gM_P.*(xVr_P(5,:)+.5*xiR2_P(5,:)).*(vP_P+.5*rhsV2_P-Ek)...
        +gH.*(xVr_P(6,:)+.5*xiR2_P(6,:)).*(vP_P+.5*rhsV2_P-Eh)+gCaPN_P.*(xVr_P(7,:)+.5*xiR2_P(7,:)).^2.*(xVr_P(8,:)+.5*xiR2_P(8,:)).*(vP_P+.5*rhsV2_P-Eca_P)...
        +gCaT_P.*(xVr_P(9,:)+.5*xiR2_P(9,:)).^2.*(xVr_P(10,:)+.5*xiR2_P(10,:)).*(vP_P+.5*rhsV2_P-Eca_P)...
        +gKCa_P*(xVr_P(11,:)+.5*xiR2_P(11,:)).*(vP_P+.5*rhsV2_P-Ek)+gDR_GP*(xVr_P(12,:)+.5*xiR2_P(12,:)).*(vP_P+.5*rhsV2_P-Ek);
    rhsV3_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV3_P-xi_PGC);
    
    xiR3_G(1,:)=dt*2.1./tauM(vP_G+.5*rhsV2_G).*(mInf(vP_G+.5*rhsV2_G)-(xVr_G(1,:)+.5*xiR2_G(1,:))); %I_Na
    xiR3_G(2,:)=dt*2.1./tauH(vP_G+.5*rhsV2_G).*(hInf(vP_G+.5*rhsV2_G)-(xVr_G(2,:)+.5*xiR2_G(2,:))); %I_Na
    xiR3_G(3,:)=dt*3.3./atau(vP_G+.5*rhsV2_G).*(aminf_g(vP_G+.5*rhsV2_G)-(xVr_G(3,:)+.5*xiR2_G(3,:))); %I_A
    xiR3_G(4,:)=dt*3.3./hatau_g(vP_G+.5*rhsV2_G).*(ahinf_g(vP_G+.5*rhsV2_G)-(xVr_G(4,:)+.5*xiR2_G(4,:))); %I_A
    xiR3_G(5,:)=dt./tauMusc(vP_G+.5*rhsV2_G).*(minfMusc(vP_G+.5*rhsV2_G)-(xVr_G(5,:)+.5*xiR2_G(5,:))); %I_M
    xiR3_G(6,:)=dt./capmtau(vP_G+.5*rhsV2_G).*(capminf(vP_G+.5*rhsV2_G)-(xVr_G(6,:)+.5*xiR2_G(6,:))); %I_CaPN
    xiR3_G(7,:)=dt./caphtau(vP_G+.5*rhsV2_G).*(caphinf(vP_G+.5*rhsV2_G)-(xVr_G(7,:)+.5*xiR2_G(7,:))); %I_CaPN
    xiR3_G(8,:)=dt./cattaumG(vP_G+.5*rhsV2_G).*(catminfG(vP_G+.5*rhsV2_G)-(xVr_G(8,:)+.5*xiR2_G(8,:))); %I_CaT
    xiR3_G(9,:)=dt./cattauh(vP_G+.5*rhsV2_G).*(cathinf(vP_G+.5*rhsV2_G)-(xVr_G(9,:)+.5*xiR2_G(9,:))); %I_CaT
    xiR3_G(10,:)=dt./cantau(vP_G+.5*rhsV2_G).*(canminf(vP_G+.5*rhsV2_G)-(xVr_G(10,:)+.5*xiR2_G(10,:))); %I_CAN
    xiR3_G(11,:)=dt*(kcaa(vP_G+.5*rhsV2_G,caCon_G+.5*caCr2_G).*(1-(xVr_G(11,:)+.5*xiR2_G(11,:))) - 0.05*(xVr_G(11,:)+.5*xiR2_G(11,:))); %I_KCa
    xiR3_G(12,:)=dt*3.3*(mss_dr(vP_G+.5*rhsV2_G)-(xVr_G(12,:)+.5*xiR2_G(12,:)))./taum_dr(vP_G+.5*rhsV2_G); %I_DR
    
    xiR3_M(1,:)=dt*(ana(vP_M+.5*rhsV2_M).*(1-(xVr_M(1,:)+.5*xiR2_M(1,:))) - bna(vP_M+.5*rhsV2_M).*(xVr_M(1,:)+.5*xiR2_M(1,:)));
    xiR3_M(2,:)=dt*(ahna(vP_M+.5*rhsV2_M).*(1-(xVr_M(2,:)+.5*xiR2_M(2,:))) - bhna(vP_M+.5*rhsV2_M).*(xVr_M(2,:)+.5*xiR2_M(2,:)));
    xiR3_M(3,:)=dt*3.3./atau(vP_M+.5*rhsV2_M).*(aminf(vP_M+.5*rhsV2_M)-(xVr_M(3,:)+.5*xiR2_M(3,:))); %I_A
    xiR3_M(4,:)=dt*3.3./hatau(vP_M+.5*rhsV2_M).*(ahinf(vP_M+.5*rhsV2_M)-(xVr_M(4,:)+.5*xiR2_M(4,:))); %I_A
    xiR3_M(5,:)=dt/10*(ksminf(vP_M+.5*rhsV2_M)-(xVr_M(5,:)+.5*xiR2_M(5,:))); %I_KS
    xiR3_M(6,:)=dt./kshtau(vP_M+.5*rhsV2_M).*(kshinf(vP_M+.5*rhsV2_M)-(xVr_M(6,:)+.5*xiR2_M(6,:)));
    xiR3_M(7,:)=dt*(cala(vP_M+.5*rhsV2_M).*(1-(xVr_M(7,:)+.5*xiR2_M(7,:))) - calb(vP_M+.5*rhsV2_M).*(xVr_M(7,:)+.5*xiR2_M(7,:))); %I_CaL
    xiR3_M(8,:)=dt*(calha(vP_M+.5*rhsV2_M).*(1-(xVr_M(8,:)+.5*xiR2_M(8,:))) - calhb(vP_M+.5*rhsV2_M).*(xVr_M(8,:)+.5*xiR2_M(8,:)));
    xiR3_M(9,:)=dt*(kcaa(vP_M+.5*rhsV2_M,caCon_M+.5*caCr2_M).*(1-(xVr_M(9,:)+.5*xiR2_M(9,:))) - 0.05*(xVr_M(9,:)+.5*xiR2_M(9,:)));
    xiR3_M(10,:)=dt*(nss_dr(vP_M+.5*rhsV2_M)-(xVr_M(10,:)+.5*xiR2_M(10,:)))./taun_dr(vP_M+.5*rhsV2_M);    %I_DR
    xiR3_M(11,:)=dt*(kss_dr(vP_M+.5*rhsV2_M)-(xVr_M(11,:)+.5*xiR2_M(11,:)))./50;
    
    xiR3_P(1,:)=dt*2.1./tauM(vP_P+.5*rhsV2_P).*(mInf(vP_P+.5*rhsV2_P)-(xVr_P(1,:)+.5*xiR2_P(1,:))); %I_Na
    xiR3_P(2,:)=dt*2.1./tauH(vP_P+.5*rhsV2_P).*(hInf(vP_P+.5*rhsV2_P)-(xVr_P(2,:)+.5*xiR2_P(2,:))); %I_Na
    xiR3_P(3,:)=dt*3.3./atau(vP_P+.5*rhsV2_P).*(aminf_g(vP_P+.5*rhsV2_P)-(xVr_P(3,:)+.5*xiR2_P(3,:))); %I_A
    xiR3_P(4,:)=dt*3.3./hatau_g(vP_P+.5*rhsV2_P).*(ahinf_g(vP_P+.5*rhsV2_P)-(xVr_P(4,:)+.5*xiR2_P(4,:))); %I_A
    xiR3_P(5,:)=dt./tauMusc(vP_P+.5*rhsV2_P).*(minfMusc(vP_P+.5*rhsV2_P)-(xVr_P(5,:)+.5*xiR2_P(5,:))); %I_M
    xiR3_P(6,:)=dt*2.1./hcurmtau(vP_P+.5*rhsV2_P).*(hcurminf(vP_P+.5*rhsV2_P)-(xVr_P(6,:))+.5*xiR2_P(6,:)); %I_H
    xiR3_P(7,:)=dt./capmtau(vP_P+.5*rhsV2_P).*(capminf(vP_P+.5*rhsV2_P)-(xVr_P(7,:)+.5*xiR2_P(7,:))); %I_CaPN
    xiR3_P(8,:)=dt./caphtau(vP_P+.5*rhsV2_P).*(caphinf(vP_P+.5*rhsV2_P)-(xVr_P(8,:)+.5*xiR2_P(8,:))); %I_CaPN
    xiR3_P(9,:)=dt./cattaumP(vP_P+.5*rhsV2_P).*(catminfP(vP_P+.5*rhsV2_P)-(xVr_P(9,:)+.5*xiR2_P(9,:))); %I_CaT
    xiR3_P(10,:)=dt./cattauh(vP_P+.5*rhsV2_P).*(cathinf(vP_P+.5*rhsV2_P)-(xVr_P(10,:)+.5*xiR2_P(10,:))); %I_CaT
    xiR3_P(11,:)=dt*(kcaa(vP_P+.5*rhsV2_P,caCon_P+.5*caCr2_P).*(1-(xVr_P(11,:)+.5*xiR2_P(11,:))) - 0.05*(xVr_P(11,:)+.5*xiR2_P(11,:))); %I_KCa
    xiR3_P(12,:)=dt*3.3*(mss_dr(vP_P+.5*rhsV2_P)-(xVr_P(12,:)+.5*xiR2_P(12,:)))./taum_dr(vP_P+.5*rhsV2_P); %I_DR
    
    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_G(6,:)+.5*xiR2_G(6,:)).^2.*(xVr_G(7,:)+.5*xiR2_G(7,:)).*(vP_G+.5*rhsV2_G-Eca_G)+gCaT_G*(xVr_G(8,:)+.5*xiR2_G(8,:)).^2.*(xVr_G(9,:)+.5*xiR2_G(9,:)).*(vP_G+.5*rhsV2_G-Eca_G);
    caCr3_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+.5*caCr2_G))/800); %tauCa=800ms for PGC & GC
    
    Ica_M=gCaL*(xVr_M(7,:)+.5*xiR2_M(7,:)).*(xVr_M(8,:)+.5*xiR2_M(8,:)).*(vP_M+.5*rhsV2_M-Eca_M);
    caCr3_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+.5*caCr2_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_P(7,:)+.5*xiR2_P(7,:)).^2.*(xVr_P(8,:)+.5*xiR2_P(8,:)).*(vP_P+.5*rhsV2_P-Eca_P)+gCaT_P*(xVr_P(9,:)+.5*xiR2_P(9,:)).^2.*(xVr_P(10,:)+.5*xiR2_P(10,:)).*(vP_P+.5*rhsV2_P-Eca_P);
    caCr3_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+.5*caCr2_P))/800); %tauCa=800ms for PGC & GC
    
    % --- step 4 ---
    rhsV4_G=1/Rm_GM*(vP_G+rhsV3_G-El_GM)+gNa_GP*(xVr_G(1,:)+xiR3_G(1,:)).^3.*(xVr_G(2,:)+xiR3_G(2,:)).*(vP_G+rhsV3_G-Ena)+...
        gA_G*(xVr_G(3,:)+xiR3_G(3,:)).*(xVr_G(4,:)+xiR3_G(4,:)).*(vP_G+rhsV3_G-Ek)+gM_G.*(xVr_G(5,:)+xiR3_G(5,:)).*(vP_G+rhsV3_G-Ek)+...
        gCaPN_G.*(xVr_G(6,:)+xiR3_G(6,:)).^2.*(xVr_G(7,:)+xiR3_G(7,:)).*(vP_G+rhsV3_G-Eca_G)+gCaT_G.*(xVr_G(8,:)+xiR3_G(8,:)).^2.*(xVr_G(9,:)+xiR3_G(9,:)).*(vP_G+rhsV3_G-Eca_G)...
        +gCAN.*((caCon_G+caCr3_G)/(200+(caCon_G+caCr3_G))).*(xVr_G(10,:)+xiR3_G(10,:)).*(vP_G+rhsV3_G-Ecat)+gKCa_G*(xVr_G(11,:)+xiR3_G(11,:)).*(vP_G+rhsV3_G-Ek)...
        +gDR_GP*(xVr_G(12,:)+xiR3_G(12,:)).*(vP_G+rhsV3_G-Ek);
    rhsV4_G=dt/C_G*(Iapp_G-Isyn_MCG-rhsV4_G-xi_GC);
    rhsV4_M=1/Rm_GM*(vP_M+rhsV3_M-El_GM)+gNa_M*(xVr_M(1,:)+xiR3_M(1,:)).^3.*(xVr_M(2,:)+xiR3_M(2,:)).*(vP_M+rhsV3_M-Ena)+...
        gNaP*minf(vP_M+rhsV3_M).*(vP_M+rhsV3_M-Ena)+gA_M*(xVr_M(3,:)+xiR3_M(3,:)).*(xVr_M(4,:)+...
        xiR3_M(4,:)).*(vP_M+rhsV3_M-Ek)+gKS*(xVr_M(5,:)+xiR3_M(5,:)).*(xVr_M(6,:)+xiR3_M(6,:)).*(vP_M+rhsV3_M-Ek)+...
        gCaL*(xVr_M(7,:)+xiR3_M(7,:)).*(xVr_M(8,:)+xiR3_M(8,:)).*(vP_M+rhsV3_M-Eca_M)+gKCa_M*(xVr_M(9,:)+...
        xiR3_M(9,:)).*(vP_M+rhsV3_M-Ek)+gDR_M*(xVr_M(10,:)+xiR3_M(10,:)).^2.*(xVr_M(11,:)+xiR3_M(11,:)).*(vP_M+rhsV3_M-Ek);
    rhsV4_M=dt/C_MP*(Iapp_M-Isyn_GC-Isyn_PGC-rhsV4_M-xi_MC);
    rhsV4_P=1/Rm_P*(vP_P+rhsV3_P-El_P)+gNa_GP*(xVr_P(1,:)+xiR3_P(1,:)).^3.*(xVr_P(2,:)+xiR3_P(2,:)).*(vP_P+rhsV3_P-Ena)+...
        gA_P*(xVr_P(3,:)+xiR3_P(3,:)).*(xVr_P(4,:)+xiR3_P(4,:)).*(vP_P+rhsV3_P-Ek)+gM_P.*(xVr_P(5,:)+xiR3_P(5,:)).*(vP_P+rhsV3_P-Ek)...
        +gH.*(xVr_P(6,:)+xiR3_P(6,:)).*(vP_P+rhsV3_P-Eh)+gCaPN_P.*(xVr_P(7,:)+xiR3_P(7,:)).^2.*(xVr_P(8,:)+xiR3_P(8,:)).*(vP_P+rhsV3_P-Eca_P)...
        +gCaT_P.*(xVr_P(9,:)+xiR3_P(9,:)).^2.*(xVr_P(10,:)+xiR3_P(10,:)).*(vP_P+rhsV3_P-Eca_P)...
        +gKCa_P*(xVr_P(11,:)+xiR3_P(11,:)).*(vP_P+rhsV3_P-Ek)+gDR_GP*(xVr_P(12,:)+xiR3_P(12,:)).*(vP_P+rhsV3_P-Ek);
    rhsV4_P=dt/C_MP*(Iapp_P-Isyn_MCP-rhsV4_P-xi_PGC);
    
    xiR4_G(1,:)=dt*2.1./tauM(vP_G+rhsV3_G).*(mInf(vP_G+rhsV3_G)-(xVr_G(1,:)+xiR3_G(1,:))); %I_Na
    xiR4_G(2,:)=dt*2.1./tauH(vP_G+rhsV3_G).*(hInf(vP_G+rhsV3_G)-(xVr_G(4,:)+xiR3_G(4,:))); %I_Na
    xiR4_G(3,:)=dt*3.3./atau(vP_G+rhsV3_G).*(aminf_g(vP_G+rhsV3_G)-(xVr_G(3,:)+xiR3_G(3,:))); %I_A
    xiR4_G(4,:)=dt*3.3./hatau_g(vP_G+rhsV3_G).*(ahinf_g(vP_G+rhsV3_G)-(xVr_G(4,:)+xiR3_G(4,:))); %I_A
    xiR4_G(5,:)=dt./tauMusc(vP_G+rhsV3_G).*(minfMusc(vP_G+rhsV3_G)-(xVr_G(5,:)+xiR3_G(5,:))); %I_M
    xiR4_G(6,:)=dt./capmtau(vP_G+rhsV3_G).*(capminf(vP_G+rhsV3_G)-(xVr_G(6,:)+xiR3_G(6,:))); %I_CaPN
    xiR4_G(7,:)=dt./caphtau(vP_G+rhsV3_G).*(caphinf(vP_G+rhsV3_G)-(xVr_G(7,:)+xiR3_G(7,:))); %I_CaPN
    xiR4_G(8,:)=dt./cattaumG(vP_G+rhsV3_G).*(catminfG(vP_G+rhsV3_G)-(xVr_G(8,:)+xiR3_G(8,:))); %I_CaT
    xiR4_G(9,:)=dt./cattauh(vP_G+rhsV3_G).*(cathinf(vP_G+rhsV3_G)-(xVr_G(9,:)+xiR3_G(9,:))); %I_CaT
    xiR4_G(10,:)=dt./cantau(vP_G+rhsV3_G).*(canminf(vP_G+rhsV3_G)-(xVr_G(10,:)+xiR3_G(10,:))); %I_CAN
    xiR4_G(11,:)=dt*(kcaa(vP_G+rhsV3_G,caCon_G+caCr3_G).*(1-(xVr_G(11,:)+xiR3_G(11,:))) - 0.05*(xVr_G(11,:)+xiR3_G(11,:))); %I_KCa
    xiR4_G(12,:)=dt*3.3*(mss_dr(vP_G+rhsV3_G)-(xVr_G(12,:)+xiR3_G(12,:)))./taum_dr(vP_G+rhsV3_G); %I_DR
    
    xiR4_M(1,:)=dt*(ana(vP_M+rhsV3_M).*(1-(xVr_M(1,:)+xiR3_M(1,:))) - bna(vP_M+rhsV3_M).*(xVr_M(1,:)+xiR3_M(1,:)));
    xiR4_M(2,:)=dt*(ahna(vP_M+rhsV3_M).*(1-(xVr_M(2,:)+xiR3_M(2,:))) - bhna(vP_M+rhsV3_M).*(xVr_M(2,:)+xiR3_M(2,:)));
    xiR4_M(3,:)=dt*3.3./atau(vP_M+rhsV3_M).*(aminf(vP_M+rhsV3_M)-(xVr_M(3,:)+xiR3_M(3,:))); %I_A
    xiR4_M(4,:)=dt*3.3./hatau(vP_M+rhsV3_M).*(ahinf(vP_M+rhsV3_M)-(xVr_M(4,:)+xiR3_M(4,:))); %I_A
    xiR4_M(5,:)=dt/10*(ksminf(vP_M+rhsV3_M)-(xVr_M(5,:)+xiR3_M(5,:))); %I_KS
    xiR4_M(6,:)=dt./kshtau(vP_M+rhsV3_M).*(kshinf(vP_M+rhsV3_M)-(xVr_M(6,:)+xiR3_M(6,:)));
    xiR4_M(7,:)=dt*(cala(vP_M+rhsV3_M).*(1-(xVr_M(7,:)+xiR3_M(7,:))) - calb(vP_M+rhsV3_M).*(xVr_M(7,:)+xiR3_M(7,:))); %I_CaL
    xiR4_M(8,:)=dt*(calha(vP_M+rhsV3_M).*(1-(xVr_M(8,:)+xiR3_M(8,:))) - calhb(vP_M+rhsV3_M).*(xVr_M(8,:)+xiR3_M(8,:)));
    xiR4_M(9,:)=dt*(kcaa(vP_M+rhsV3_M,caCon_M+caCr3_M).*(1-(xVr_M(9,:)+xiR3_M(9,:))) - 0.05*(xVr_M(9,:)+xiR3_M(9,:)));
    xiR4_M(10,:)=dt*(nss_dr(vP_M+rhsV3_M)-(xVr_M(10,:)+xiR3_M(10,:)))./taun_dr(vP_M+rhsV3_M);    %I_DR
    xiR4_M(11,:)=dt*(kss_dr(vP_M+rhsV3_M)-(xVr_M(11,:)+xiR3_M(11,:)))./50;
    
    xiR4_P(1,:)=dt*2.1./tauM(vP_P+rhsV3_P).*(mInf(vP_P+rhsV3_P)-(xVr_P(1,:)+xiR3_P(1,:))); %I_Na
    xiR4_P(2,:)=dt*2.1./tauH(vP_P+rhsV3_P).*(hInf(vP_P+rhsV3_P)-(xVr_P(2,:)+xiR3_P(2,:))); %I_Na
    xiR4_P(3,:)=dt*3.3./atau(vP_P+rhsV3_P).*(aminf_g(vP_P+rhsV3_P)-(xVr_P(3,:)+xiR3_P(3,:))); %I_A
    xiR4_P(4,:)=dt*3.3./hatau_g(vP_P+rhsV3_P).*(ahinf_g(vP_P+rhsV3_P)-(xVr_P(4,:)+xiR3_P(4,:))); %I_A
    xiR4_P(5,:)=dt./tauMusc(vP_P+rhsV3_P).*(minfMusc(vP_P+rhsV3_P)-(xVr_P(5,:)+xiR3_P(5,:))); %I_M
    xiR4_P(6,:)=dt*2.1./hcurmtau(vP_P+rhsV3_P).*(hcurminf(vP_P+rhsV3_P)-(xVr_P(6,:))+xiR3_P(6,:)); %I_H
    xiR4_P(7,:)=dt./capmtau(vP_P+rhsV3_P).*(capminf(vP_P+rhsV3_P)-(xVr_P(7,:)+xiR3_P(7,:))); %I_CaPN
    xiR4_P(8,:)=dt./caphtau(vP_P+rhsV3_P).*(caphinf(vP_P+rhsV3_P)-(xVr_P(8,:)+xiR3_P(8,:))); %I_CaPN
    xiR4_P(9,:)=dt./cattaumP(vP_P+rhsV3_P).*(catminfP(vP_P+rhsV3_P)-(xVr_P(9,:)+xiR3_P(9,:))); %I_CaT
    xiR4_P(10,:)=dt./cattauh(vP_P+rhsV3_P).*(cathinf(vP_P+rhsV3_P)-(xVr_P(10,:)+xiR3_P(10,:))); %I_CaT
    xiR4_P(11,:)=dt*(kcaa(vP_P+rhsV3_P,caCon_P+caCr3_P).*(1-(xVr_P(11,:)+xiR3_P(11,:))) - 0.05*(xVr_P(11,:)+xiR3_P(11,:))); %I_KCa
    xiR4_P(12,:)=dt*3.3*(mss_dr(vP_P+rhsV3_P)-(xVr_P(12,:)+xiR3_P(12,:)))./taum_dr(vP_P+rhsV3_P); %I_DR
    
    %update calcium concentration
    Ica_G=gCaPN_G*(xVr_G(6,:)+xiR3_G(6,:)).^2.*(xVr_G(7,:)+xiR3_G(7,:)).*(vP_G+rhsV3_G-Eca_G)+gCaT_G*(xVr_G(8,:)+xiR3_G(8,:)).^2.*(xVr_G(9,:)+xiR3_G(9,:)).*(vP_G+rhsV3_G-Eca_G);
    caCr4_G=dt*(-Ica_G/(2*F*thick_GP) + (.05-(caCon_G+caCr3_G))/800); %tauCa=800ms for PGC & GC

    Ica_M=gCaL*(xVr_M(7,:)+xiR3_M(7,:)).*(xVr_M(8,:)+xiR3_M(8,:)).*(vP_M+rhsV3_M-Eca_M);
    caCr4_M=dt*(-Ica_M/(2*F*thick_M) + (.05-(caCon_M+caCr3_M))/10); %tauCa=10ms for MC
    
    Ica_P=gCaPN_P*(xVr_P(7,:)+xiR3_P(7,:)).^2.*(xVr_P(8,:)+xiR3_P(8,:)).*(vP_P+rhsV3_P-Eca_P)+gCaT_P*(xVr_P(9,:)+xiR3_P(9,:)).^2.*(xVr_P(10,:)+xiR3_P(10,:)).*(vP_P+rhsV3_P-Eca_P);
    caCr4_P=dt*(-Ica_P/(2*F*thick_GP) + (.05-(caCon_P+caCr3_P))/800); %tauCa=800ms for PGC & GC
    
    % FINAL Step... update voltage
    vP_G=vP_G+1/6*(rhsV1_G+2*rhsV2_G+2*rhsV3_G+rhsV4_G);
    vP_M=vP_M+1/6*(rhsV1_M+2*rhsV2_M+2*rhsV3_M+rhsV4_M);
    vP_P=vP_P+1/6*(rhsV1_P+2*rhsV2_P+2*rhsV3_P+rhsV4_P);
    
    %update gating variables
    xVr_G=xVr_G+1/6*(xiR1_G+2*xiR2_G+2*xiR3_G+xiR4_G);
    xVr_M=xVr_M+1/6*(xiR1_M+2*xiR2_M+2*xiR3_M+xiR4_M);
    xVr_P=xVr_P+1/6*(xiR1_P+2*xiR2_P+2*xiR3_P+xiR4_P);
    
    %update calcium concentration
    caCon_G=caCon_G+1/6*(caCr1_G+2*caCr2_G+2*caCr3_G+caCr4_G);
    caCon_M=caCon_M+1/6*(caCr1_M+2*caCr2_M+2*caCr3_M+caCr4_M);
    caCon_P=caCon_P+1/6*(caCr1_P+2*caCr2_P+2*caCr3_P+caCr4_P);
    
    %update calcium reversal potential, Nernst eqn
    Eca_G=RTovF./2*log(10./caCon_G); %get 70mV when caCon=0.05
    Eca_M=RTovF./2*log(10./caCon_M); %get 70mV when caCon=0.05
    Eca_P=RTovF./2*log(10./caCon_P); %get 70mV when caCon=0.05
    
    %spiking
    idGspk=(vP_G>vltThres).*(j-tLstSp_G>=mnTim);
    idMspk=(vP_M>vltThres).*(j-tLstSp_M>=mnTim);
    idPspk=(vP_P>vltThres).*(j-tLstSp_P>=mnTim);
    %update spiking vars
    if(sum(idGspk)>0) % at least 1 G spikes
        tLstSp_G(idGspk==1)=j;
        spTimes_G=[spTimes_G; find(idGspk') j*dt*ones(sum(idGspk),1)];
    end
    if(sum(idMspk)>0)
        tLstSp_M(idMspk==1)=j;
        spTimes_M=[spTimes_M; find(idMspk') j*dt*ones(sum(idMspk),1)];
    end
    if(sum(idPspk)>0)
        tLstSp_P(idPspk==1)=j;
        spTimes_P=[spTimes_P; find(idPspk') j*dt*ones(sum(idPspk),1)];
    end
    if(~isreal(vP_G) || ~isreal(vP_M) || ~isreal(vP_P)) %imag voltage
       flgGd=0; 
       return
    end
end
%Caluclating spike count
if(~isempty(spTimes_G))
    spTimes_G1=spTimes_G(find(spTimes_G(:,1)==1),2);
    spTimes_G2=spTimes_G(find(spTimes_G(:,1)==2),2);
    spTimes_G3=spTimes_G(find(spTimes_G(:,1)==3),2);
end
if(~isempty(spTimes_M))
    spTimes_M1=spTimes_M(find(spTimes_M(:,1)==1),2);
    spTimes_M2=spTimes_M(find(spTimes_M(:,1)==2),2);
end
if(~isempty(spTimes_P))
    spTimes_P1=spTimes_P(find(spTimes_P(:,1)==1),2);
    spTimes_P2=spTimes_P(find(spTimes_P(:,1)==2),2);
end
end

    %Ina m alph & bet
    function a_na=ana_g(v)
        a_na=.4*(v+25)./(1-exp(-(v+25)./7.2));
    end
    function b_na=bna_g(v)
        b_na=-.124*(v+25)./(1-exp((v+25)./7.2));
    end
    function x=mInf(v)
       x=ana_g(v)./(ana_g(v)+bna_g(v));
    end
    function x=tauM(v)
        x=max(1./(ana_g(v)+bna_g(v)),0.02);
    end
    function a_na=ana(v)
        a_na=.32*(v+45)./(1-exp(-(v+45)./4));
    end
    function b_na=bna(v)
        b_na=-.28*(v+18)./(1-exp((v+18)./5));
    end
    %Ina h alph & bet
    function a_hna=ahna_g(v)
        a_hna=.03*(v+40)./(1-exp(-(v+40)./1.5));
    end
    function b_hna=bhna_g(v)
        b_hna=-.01*(v+40)./(1-exp((v+40)./1.5));
    end
    function x=hInf(v)
       x=1./(1+exp((v+45)./4));
    end
    function x=tauH(v)
        x=max(1./(ahna_g(v)+bhna_g(v)),0.5);
    end
    function a_hna=ahna(v)
        a_hna=.128./exp((v+41)./18);
    end
    function b_hna=bhna(v)
        b_hna=4./(1+exp(-(v+18)./5));
    end
    %InaP
    function m_inf=minf(v)
        m_inf=1./(1+exp(-(v+50)./5));
    end
    %I_M
    function x=minfMusc(v)
        x=1./(1+exp(-(v+35)./5));
    end
    function x=tauMusc(v)
        x=1000./(3.3*exp((v+35)./40)+exp(-(v+35)./20));
    end
    %I_A
    function a_minf=aminf_g(v)
       a_minf=1./(1+exp(-(v-7.6)./14));
    end
    function a_minf=aminf(v)
       a_minf=1./(1+exp(-(v-17.5)./14));
    end
    function a_tau=atau(v)
        a_tau=25*exp((v+45)./13.3)./(1+exp((v+45)./10));
    end
    function a_hinf=ahinf_g(v)
        a_hinf=1./(1+exp((v+67.4)./6));
    end
    function ha_tau=hatau_g(v)
        ha_tau=138.8*exp((v+70)./5.1)./(1+exp((v+70)./5));
    end
    function a_hinf=ahinf(v)
        a_hinf=1./(1+exp((v+41.7)./6));
    end
    function ha_tau=hatau(v)
        ha_tau=55.5*exp((v+70)./5.1)./(1+exp((v+70)./5));
    end
    %I_KS (tau=10)
    function ks_minf=ksminf(v)
        ks_minf=1./(1+exp(-(v+34)./6.5));
    end
    function ks_hinf=kshinf(v)
        ks_hinf=1./(1+exp((v+68)./6.6));
    end
    function ks_htau=kshtau(v)
        ks_htau=200+330./(1+exp(-(v+71.6)./6.85));
    end
    % H current
    function x=hcurminf(v)
        x=1./(1+exp((v+80)./10));
    end
    function x=hcurmtau(v)
        x=1176.5*exp((v+65)./23.5)./(1+exp((v+65)./11.8));
    end
    %I_CaL
    function cal_a=cala(v)
        cal_a=7.5./(1+exp(-(v-13)./7));
    end
    function cal_b=calb(v)
        cal_b=1.65./(1+exp((v-14)./4));
    end
    function cal_ha=calha(v)
        cal_ha=.0068./(1+exp((v+30)./12));
    end
    function cal_hb=calhb(v)
        cal_hb=.06./(1+exp(-v./11));
    end
%I_CaPN
    function x=capminf(v)
        x=1./(1+exp(-(v+10)./4));
    end
    function x=capmtau(v)
        x=.4+.7./(exp(-(v+5)./15)+exp((v+5)./15));
    end
    function x=caphinf(v)
        x=1./(1+exp((v+25)./2));
    end
    function x=caphtau(v)
        x=300+100./(exp(-(v+40)./9.5)+exp((v+40)./9.5));
    end
    %I_CaT
    function x=catminfG(v) %change for rev PGC!
        x=1./(1+exp(-(v+44)./5.5)); %(v+59) for PGC
    end
    function x=cattaumG(v) %change for rev PGC!
        x=1.5+3.5./(exp(-(v+30)./15)+exp((v+30)./15)); %(v+45) for PGC
    end
    function x=catminfP(v) %change for rev PGC!
        x=1./(1+exp(-(v+59)./5.5)); %(v+59) for PGC
    end
    function x=cattaumP(v) %change for rev PGC!
        x=1.5+3.5./(exp(-(v+45)./15)+exp((v+45)./15)); %(v+45) for PGC
    end
    function x=cathinf(v)
        x=1./(1+exp((v+70)./4)); 
    end
    function x=cattauh(v) 
        x=10+40./(exp(-(v+50)./15)+exp((v+50)./15));
    end
    %I_CAN
    function x=canminf(v)
        x=1./(1+exp(-(v+43)./5.2));
    end
    function x=cantau(v)
        x=1.6+2.7./(exp(-(v+55)./15)+exp((v+55)./15));
    end
    %I_KCA, beta_m=.05
    function kca_a=kcaa(v,CA)
        kca_a=-500*exp((v-65)./27).*(.015-CA)./(1-exp(-(CA-.015)./.0013));
    end
    %I_DR fast delayed rectifier; parms from getKssParm.m, stored in parm_fstDR.mat
    function xx=mss_dr(v)
        xx=1./(1+exp(-(v-21)./10));
    end
    function xx=taum_dr(v)
        xx=285.7*exp((v+50)./36.4)./(1+exp((v+50)./18.2));
    end
    function xx=kss_dr(v)
        xx=.43315*(1+tanh(-(v+13.925)./13.0215))+.1337;
    end
    function xx=nss_dr(v)
        xx=((v+100)./150).^8.5849./(0.5747^8.5849+((v+100)./150).^8.5849);
    end
    function xx=taun_dr(v)
        xx=1./(.27654./exp((v+29.9998)./66.3783)+2.89./(1+exp(-(v-19.0524)./12.8786)));
    end
    function ds_GABA=dsGABA(s,V)
        ds_GABA=(1/1.25)*FGP(V).*(1-s)-(1/18).*s;
    end
    function F_GP=FGP(Vpre)
        F_GP=1./(1+exp(-(Vpre+40)./2)); %sigma=0.2 from MC; sigma=2 to MC 
    end
    function ds_AMPA=dsAMPA(s,V)
        ds_AMPA=FM(V).*(1-s)-(1/5.5)*s;
    end
    function ds_NMDA=dsNMDA(s,V)
        ds_NMDA=(1/52)*FM(V).*(1-s)-(1/343)*s;
    end
    function x=B(V)
        x=(1+exp(-0.062*V)./3.57).^-1;
    end
    function F_M=FM(Vpre)
        F_M=1./(1+exp(-(Vpre)./0.2)); %sigma=0.2 from MC; sigma=2 to MC 
    end

end
% script to check formulas for stats of bg syn
% doing Nt realizations
stim=0; % 0 = ortho; 1 = retro

Nt=1500;

rng('shuffle');

timeVars=struct('dt',0.01,'tSpon',1200,'tEvok',1900);
tevok=(0:timeVars.dt:timeVars.tEvok)';
tspon=(-timeVars.tSpon:timeVars.dt:0)';
tEnd=timeVars.tSpon+timeVars.tEvok;
Lt=round((tEnd)/timeVars.dt)+1;
tm=(-timeVars.tSpon:timeVars.dt:timeVars.tEvok)';

%time varying rate and correlation of PP kicks
if stim==0 %ortho
    lamO=struct('tsft',25.1,'lmOsp',.002,'lmOevk',15/367,'tauO',60,'tauO2',61,'ssV',.007);
    crEv=0.18*(-(tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO) + (tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO2))+.9;
    crM=[0.2*ones(length(tspon)-1,1); crEv];
    tauO=lamO.tauO;
    tauO2=lamO.tauO2;
    ssV=lamO.ssV;
    lmOevk=lamO.lmOevk;
    lmOsp=lamO.lmOsp;
    tsft=lamO.tsft;
    tmShift=50; %in ms
    tevokS=[tevok; tevok(end)+(timeVars.dt:timeVars.dt:tmShift)']; %take nuevok out further by tmShift
    nmShift=length((timeVars.dt:timeVars.dt:tmShift)); %# elem must shift by
    nuevokO=lmOevk*(-(tevokS+tsft).*exp(-(tevokS+tsft)/tauO) + (tevokS+tsft).*exp(-(tevokS+tsft)/tauO2))+ssV;
    nuevokO=nuevokO(nmShift+1:end); %shift over by tmShift
    nu_A=[lmOsp*ones(length(tspon)-1,1);nuevokO];
elseif stim==1 %retro
    lamA=struct('tsft',25.1,'lmAspon',.002,'lmAevk',.095/367,'tau_A',450);
    crEvr=0.2+0.5*(1-exp(-(tevok(10001:end)-tevok(10001))/100)); %Jan 2021
    crM=[0.2*ones(length(tspon)+9999,1);crEvr];
    tau_A=lamA.tau_A;
    lmAevk=lamA.lmAevk;
    lmAspon=lamA.lmAspon;
    tsft=lamA.tsft;
    nuevok_A=(tevok+tsft).*exp(-(tevok+tsft)/tau_A);
    nu_A=[lmAspon*ones(length(tspon)-1,1);lmAevk*nuevok_A];
else
    disp('Variable stim needs to be set to either 0 or 1')
end

tau_sM=[10; 10]; %ms for MC
a=ones(2,1); %size of jumps, ONLY 2 syna
nu_A=[nu_A nu_A]; %making identical nu_A but CAN be different!

numinA=min(nu_A(:,1),nu_A(:,2)).*crM;
nuaugA_1=nu_A(:,1)-numinA;
nuaugA_2=nu_A(:,2)-numinA;

randA_M1=rand(Nt,Lt);
randA_M2=rand(Nt,Lt);
randA_Mc=rand(Nt,Lt);

%state variables
sA_MC1=zeros(Nt,1);
sA_MC2=zeros(Nt,1);
%keep running sum of SIMS
mnM_sA=zeros(2,Lt);
varM_sA=zeros(2,Lt);
covM_sA=zeros(1,Lt);
%analytic form
mnM_Ana=zeros(2,Lt);
vrM_Ana=zeros(2,Lt);
cvM_Ana=zeros(1,Lt);

for j=2:Lt
    
    sA_MC1=sA_MC1+timeVars.dt*(-sA_MC1)./tau_sM(1);
    sA_MC1=sA_MC1+a(1)*(randA_M1(:,j)<nuaugA_1(j)*timeVars.dt)+a(1)*(randA_Mc(:,j)<numinA(j)*timeVars.dt);
    
    sA_MC2=sA_MC2+timeVars.dt*(-sA_MC2)./tau_sM(2);
    sA_MC2=sA_MC2+a(2)*(randA_M2(:,j)<nuaugA_2(j)*timeVars.dt)+a(2)*(randA_Mc(:,j)<numinA(j)*timeVars.dt);
   
    %calc avg, keep running sums for 2nd order stats
    mnM_sA(1,j)=mean(sA_MC1);
    mnM_sA(2,j)=mean(sA_MC1);
    varM_sA(1,j)=sum(sA_MC1.^2)./(Nt-1);
    varM_sA(2,j)=sum(sA_MC2.^2)./(Nt-1);
    covM_sA(j)=sum(sA_MC1.*sA_MC2)./(Nt-1);
    
    %analytic time-step
    mnM_Ana(:,j)=mnM_Ana(:,j-1)+timeVars.dt*(-mnM_Ana(:,j-1)./tau_sM+a.*(nu_A(j,:)'));
    vrM_Ana(:,j)=vrM_Ana(:,j-1)+timeVars.dt*(-vrM_Ana(:,j-1)*2./tau_sM+a.^2.*(nu_A(j,:)'));
    cvM_Ana(:,j)=cvM_Ana(:,j-1)+timeVars.dt*(-cvM_Ana(:,j-1)*2/tau_sM(1)+crM(j)*a(1).^2.*min(nu_A(j,:)));
end
%center stats
varM_sA=varM_sA-Nt/(Nt-1)*mnM_sA.^2;
covM_sA=covM_sA-Nt/(Nt-1)*(mnM_sA(1,:).*mnM_sA(2,:));

figure
hold on
plot(tm,mnM_sA(1,:),'k--','LineWidth',2)
plot(tm,mnM_Ana(1,:),'r-')
set(gca,'FontSize',18)
xlabel('Time (ms)')
ylabel('Mean s(t)')
title('MC<->MC Mean')

figure
hold on
plot(tm,varM_sA(1,:),'k--','LineWidth',2)
plot(tm,vrM_Ana(1,:),'r-')
set(gca,'FontSize',18)
xlabel('Time (ms)')
ylabel('Var s(t)')
title('MC<->MC Var')

figure
hold on
plot(tm,covM_sA(1,:),'k--','LineWidth',2)
plot(tm,cvM_Ana(1,:),'r-')
set(gca,'FontSize',18)
xlabel('Time (ms)')
ylabel('Var s(t)')
title('MC<->MC Cov')
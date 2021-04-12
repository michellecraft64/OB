% script to fit inputs (time-varying syn) to MonteCarlo of ouput stats
% Uses LinFiltShift.m

dt=0.1; %in ms
tvec=(0:dt:3000)' - 1000;
Lt=length(tvec);

%orth part
lamO=struct('tsft',25.1,'lmOsp',.002,'lmOevk',15/367,'tauO',60,'tauO2',61,'ssV',.007);
timeVars=struct('dt',0.01,'tSpon',1200,'tEvok',1900);
tevok=(0:timeVars.dt:timeVars.tEvok)';
tspon=(-timeVars.tSpon:timeVars.dt:0)';
crEv=0.18*(-(tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO) + (tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO2))+.9;
crM_t=[0.2*ones(length(tspon)-1,1); crEv];
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
nu_Ai=[lmOsp*ones(length(tspon)-1,1);nuevokO];
tme=[tspon;tevok(2:end)]; %single time vector
nu_Aort=pchip(tme,nu_Ai,tvec);
crMort=pchip(tme,crM_t,tvec);
%retr part
lamA=struct('tsft',25.1,'lmAspon',.002,'lmAevk',.095/367,'tau_A',450);
timeVars=struct('dt',0.01,'tSpon',1200,'tEvok',1900);
tevok=(0:timeVars.dt:timeVars.tEvok)';
tspon=(-timeVars.tSpon:timeVars.dt:0)';
crEvr=0.2+0.5*(1-exp(-(tevok(10001:end)-tevok(10001))/100)); %Jan 2021
crM_t=[0.2*ones(length(tspon)+9999,1);crEvr];
tau_A=lamA.tau_A;
lmAevk=lamA.lmAevk;
lmAspon=lamA.lmAspon;
tsft=lamA.tsft;
nuevok_A=(tevok+tsft).*exp(-(tevok+tsft)/tau_A);
nu_Ai=[lmAspon*ones(length(tspon)-1,1);lmAevk*nuevok_A];
nu_Aret=pchip(tme,nu_Ai,tvec);
crMret=pchip(tme,crM_t,tvec);
%% calc sE and sI ahead of time because fast to run MC_bg.m
tau_sy=[10 5.5]; %syn time-scales (ms) for MC (E&I); for GC/PGC (E&I) it is 5.5
aJmp=[1 1];  %1st entry is for E, 2nd entry is for I bg synapse
Esyn_GP=-80; %revers potent. for I
Esyn_M=0;    %revers potent. for E
% syn stats for Orth
sMeanE=zeros(Lt,2); %!! second entry is input for GC/PGC!!
    sMeanE(1,:)=nu_Aort(1)*aJmp.*tau_sy; %assume same nu_A
sMeanI=zeros(Lt,2);
    sMeanI(1,:)=0.75*nu_Aort(1)*aJmp.*tau_sy; %assume same nu_A; .* for aJmp,tau_sy dont always MATCH up
sVarE=zeros(Lt,2);
    sVarE(1,:)=nu_Aort(1)*aJmp.^2.*tau_sy./2; %assume same nu_A
sVarI=zeros(Lt,2);
    sVarI(1,:)=0.75*nu_Aort(1)*aJmp.^2.*tau_sy./2; %assume same nu_A
sCovE=zeros(Lt,1);
    sCovE(1)=nu_Aort(1)*aJmp(1)^2*tau_sy(1)/2*crMort(1); %assume same nu_A(min), tau b/c MC-MC
sCovI=zeros(Lt,1);
    sCovI(1)=0.75*nu_Aort(1)*aJmp(2)^2*tau_sy(1)/2*crMort(1); %assume same nu_A(min), tau b/c MC-MC    
%solve tau_sM*s' = -s+a*tau_sM*nu_A with trapezoidal rule
for j=2:Lt
    sMeanE(j,:)=((1-0.5*dt./tau_sy).*sMeanE(j-1,:)+0.5*dt*aJmp*(nu_Aort(j-1)+nu_Aort(j)))./...
        (1+0.5*dt./tau_sy);
    sMeanI(j,:)=((1-0.5*dt./tau_sy).*sMeanI(j-1,:)+0.5*dt*aJmp*0.75*(nu_Aort(j-1)+nu_Aort(j)))./...
        (1+0.5*dt./tau_sy);
    sVarE(j,:)=((1-0.5*dt*2./tau_sy).*sVarE(j-1,:)+0.5*dt*aJmp.^2*(nu_Aort(j-1)+nu_Aort(j)))./...
        (1+0.5*dt*2./tau_sy);
    sVarI(j,:)=((1-0.5*dt*2./tau_sy).*sVarI(j-1,:)+0.5*dt*aJmp.^2*0.75*(nu_Aort(j-1)+nu_Aort(j)))./...
        (1+0.5*dt*2./tau_sy);
    sCovE(j)=((1-0.5*dt*2/tau_sy(1))*sCovE(j-1)+0.5*dt*aJmp(1)^2*(nu_Aort(j-1)*crMort(j-1)+nu_Aort(j)*crMort(j)))/...
        (1+0.5*dt*2/tau_sy(1)); %with symmetry same as sVar but with crM*nu_A
    sCovI(j)=((1-0.5*dt*2/tau_sy(1))*sCovI(j-1)+0.5*dt*aJmp(2)^2*0.75*(nu_Aort(j-1)*crMort(j-1)+nu_Aort(j)*crMort(j)))/...
        (1+0.5*dt*2/tau_sy(1)); %with symmetry same as sVar but with crM*nu_A
end

% syn stats for Retr
sMeanEr=zeros(Lt,2); %!! second entry is input for GC/PGC!!
    sMeanEr(1,:)=nu_Aret(1)*aJmp.*tau_sy; %assume same nu_A
sMeanIr=zeros(Lt,2);
    sMeanIr(1,:)=0.75*nu_Aret(1)*aJmp.*tau_sy; %assume same nu_A; .* for aJmp,tau_sy dont always MATCH up
sVarEr=zeros(Lt,2);
    sVarEr(1,:)=nu_Aret(1)*aJmp.^2.*tau_sy./2; %assume same nu_A
sVarIr=zeros(Lt,2);
    sVarIr(1,:)=0.75*nu_Aret(1)*aJmp.^2.*tau_sy./2; %assume same nu_A
sCovEr=zeros(Lt,1);
    sCovEr(1)=nu_Aret(1)*aJmp(1)^2*tau_sy(1)/2*crMret(1); %assume same nu_A(min), tau b/c MC-MC
sCovIr=zeros(Lt,1);
    sCovIr(1)=0.75*nu_Aret(1)*aJmp(2)^2*tau_sy(1)/2*crMret(1); %assume same nu_A(min), tau b/c MC-MC    
%solve tau_sM*s' = -s+a*tau_sM*nu_A with trapezoidal rule
for j=2:Lt
    sMeanEr(j,:)=((1-0.5*dt./tau_sy).*sMeanEr(j-1,:)+0.5*dt*aJmp*(nu_Aret(j-1)+nu_Aret(j)))./...
        (1+0.5*dt./tau_sy);
    sMeanIr(j,:)=((1-0.5*dt./tau_sy).*sMeanIr(j-1,:)+0.5*dt*aJmp*0.75*(nu_Aret(j-1)+nu_Aret(j)))./...
        (1+0.5*dt./tau_sy);
    sVarEr(j,:)=((1-0.5*dt*2./tau_sy).*sVarEr(j-1,:)+0.5*dt*aJmp.^2*(nu_Aret(j-1)+nu_Aret(j)))./...
        (1+0.5*dt*2./tau_sy);
    sVarIr(j,:)=((1-0.5*dt*2./tau_sy).*sVarIr(j-1,:)+0.5*dt*aJmp.^2*0.75*(nu_Aret(j-1)+nu_Aret(j)))./...
        (1+0.5*dt*2./tau_sy);
    sCovEr(j)=((1-0.5*dt*2/tau_sy(1))*sCovEr(j-1)+0.5*dt*aJmp(1)^2*(nu_Aret(j-1)*crMret(j-1)+nu_Aret(j)*crMret(j)))/...
        (1+0.5*dt*2/tau_sy(1)); %with symmetry same as sVar but with crM*nu_A
    sCovIr(j)=((1-0.5*dt*2/tau_sy(1))*sCovIr(j-1)+0.5*dt*aJmp(2)^2*0.75*(nu_Aret(j-1)*crMret(j-1)+nu_Aret(j)*crMret(j)))/...
        (1+0.5*dt*2/tau_sy(1)); %with symmetry same as sVar but with crM*nu_A
end

% load prev calc Monte Carlo of OB model
flNameOr='../OBsc/dOrth_ct'; %assume these mat files exist
flNameRet='../OBsc/dRetr_ct';
load(flNameOr) 
psthMCor=psthMC;
EBpsthMCor=EB_psthMC;
varMCor=varMC;
EBvarMCor=EB_varMC;
covMCor=covMC;
EBcovMCor=EB_covMC;
load(flNameRet) 
psthMCret=psthMC;
EBpsthMCret=EB_psthMC;
varMCret=varMC;
EBvarMCret=EB_varMC;
covMCret=covMC;
EBcovMCret=EB_covMC;

%interp data so on same scale as prescribed sMean/Var/Cov
mcPsthOr=pchip(tme(start:finish)*1000,psthMCor(1,:)',tvec);
mcVrOr=pchip(tme(start:finish)*1000,varMCor(1,:)',tvec);
mcCvOr=pchip(tme(start:finish)*1000,covMCor',tvec);
mcPsthRet=pchip(tme(start:finish)*1000,psthMCret(1,:)',tvec);
mcVrRet=pchip(tme(start:finish)*1000,varMCret(1,:)',tvec);
mcCvRet=pchip(tme(start:finish)*1000,covMCret',tvec);

sMt_o=[sMeanE(:,1)./norm(sMeanE(:,1)) sVarE(:,1)./norm(sVarE(:,1)) sCovE./norm(sCovE)]; %normalize
sMt_r=[sMeanEr(:,1)./norm(sMeanEr(:,1)) sVarEr(:,1)./norm(sVarEr(:,1)) sCovEr./norm(sCovEr)]; %normalize

%-- need to down sample, else too much noise --
ndt=0.3; %coarser mesh, in ms
tvs=(tvec(1):ndt:tvec(end))';
nPts=round(1000/ndt+1); %go back all the way in spont state
trncPO=0.3; %amount of filter to include (P/V/C=PSTH/Var/Cov, O/R=Ortho/Retro)
trncVO=0.3;
trncCO=0.3;
trncPR=0.3;
trncVR=0.3;
trncCR=0.3;
sMnE_O=pchip(tvec,sMt_o(:,1),tvs);
sVrE_O=pchip(tvec,sMt_o(:,2),tvs);
sCvE_O=pchip(tvec,sMt_o(:,3),tvs);
mcPO=pchip(tvec,mcPsthOr,tvs);
mcVO=pchip(tvec,mcVrOr,tvs);
mcCO=pchip(tvec,mcCvOr,tvs);
[kFilt_mn_orth,Shft_mn_orth,Coeff_mn_orth,lsApprx_mn_orth]=LinFiltShift(nPts,sMnE_O,log(mcPO),trncPO);
[kFilt_vr_orth,Shft_vr_orth,Coeff_vr_orth,lsApprx_vr_orth]=LinFiltShift(nPts,sVrE_O,log(mcVO),trncVO);
[kFilt_cv_orth,Shft_cv_orth,Coeff_cv_orth,lsApprx_cv_orth]=LinFiltShift(nPts,sCvE_O,log(mcCO),trncCO);
%retro
sMnE_R=pchip(tvec,sMt_r(:,1),tvs);
sVrE_R=pchip(tvec,sMt_r(:,2),tvs);
sCvE_R=pchip(tvec,sMt_r(:,3),tvs);
mcPR=pchip(tvec,mcPsthRet,tvs);
mcVR=pchip(tvec,mcVrRet,tvs);
mcCR=pchip(tvec,mcCvRet,tvs);
[kFilt_mn_retr,Shft_mn_retr,Coeff_mn_retr,lsApprx_mn_retr]=LinFiltShift(nPts,sMnE_R,log(mcPR),trncPR);
[kFilt_vr_retr,Shft_vr_retr,Coeff_vr_retr,lsApprx_vr_retr]=LinFiltShift(nPts,sVrE_R,log(mcVR),trncVR);
[kFilt_cv_retr,Shft_cv_retr,Coeff_cv_retr,lsApprx_cv_retr]=LinFiltShift(nPts,sCvE_R,log(mcCR),trncCR);

%plot it all
figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',psthMCor(1,:),EBpsthMCor(1,:),'b')
plot(tvec./1000,mcPsthOr,'b') %interpolated result
plot(tvs./1000,exp(lsApprx_mn_orth),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('PSTH (Hz)')
title('Ortho PSTH')

figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',varMCor(1,:),EBvarMCor(1,:),'b')
plot(tvec./1000,mcVrOr,'b') %interpolated result
plot(tvs./1000,exp(lsApprx_vr_orth),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Var Spiking')
title('Ortho Var')

figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',covMCor,EBcovMCor,'b')
plot(tvec./1000,mcCvOr,'b') %interpolated result
plot(tvs./1000,exp(lsApprx_cv_orth),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cov Spiking')
title('Ortho Cov')

figure('Renderer', 'painters', 'Position', [1250 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',psthMCret(1,:),EBpsthMCret(1,:),'r')
plot(tvec./1000,mcPsthRet,'r') %interpolated result
plot(tvs./1000,exp(lsApprx_mn_retr),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('PSTH (Hz)')
title('Retro PSTH')

figure('Renderer', 'painters', 'Position', [1250 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',varMCret(1,:),EBvarMCret(1,:),'r')
plot(tvec./1000,mcVrRet,'r') %interpolated result
plot(tvs./1000,exp(lsApprx_vr_retr),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Var Spiking')
title('Retro Var')

figure('Renderer', 'painters', 'Position', [1250 1000 1500 1000])
hold on;
errorbar(tme(start:finish)',covMCret,EBcovMCret,'r')
plot(tvec./1000,mcCvRet,'r') %interpolated result
plot(tvs./1000,exp(lsApprx_cv_retr),'k--','LineWidth',3)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Cov Spiking')
title('Retro Cov')

%plot the filters
tmFilt=0:-ndt/1000:-(trncPO)+ndt/1000; %assuming trncPO same & divisible by ndt
figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
plot(tmFilt,kFilt_mn_orth(1:round(nPts*trncPO))./ndt,'b','LineWidth',3)
plot(tmFilt,kFilt_mn_retr(1:round(nPts*trncPR))./ndt,'r','LineWidth',3)
set(gca,'FontSize',18,'XLim',[-0.005 0])
xlabel('Time (s)')
legend('Ortho','Retro')
title('Mean Filters')

figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
plot(tmFilt,kFilt_vr_orth(1:round(nPts*trncVO))./ndt,'b','LineWidth',3)
plot(tmFilt,kFilt_vr_retr(1:round(nPts*trncVR))./ndt,'r','LineWidth',3)
set(gca,'FontSize',18,'XLim',[-0.005 0])
xlabel('Time (s)')
legend('Ortho','Retro')
title('Var Filters')

figure('Renderer', 'painters', 'Position', [20 1000 1500 1000])
hold on;
plot(tmFilt,kFilt_cv_orth(1:round(nPts*trncCO))./ndt,'b','LineWidth',3)
plot(tmFilt,kFilt_cv_retr(1:round(nPts*trncCR))./ndt,'r','LineWidth',3)
set(gca,'FontSize',18,'XLim',[-0.005 0])
xlabel('Time (s)')
legend('Ortho','Retro')
title('Cov Filters')
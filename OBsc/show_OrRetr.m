%Compare OB network model MC stats with Orth/Retr in vivo data

flag_Sep=0; %if=1 sep figs for Ortho & Retro, else combine
flag_Data=1; %if=1 show data, else dont show data

flNameOr='dOrth_ct'; %assume these mat files exist
flNameRet='dRetr_ct';

% load empirical data for plot overlay
load('expDat_orEB_Rat1.mat','mn_allOr','vr_allOr','fano_eachOr',...
    'cov_allOr','tme','Twin','mn_allRet','vr_allRet','fano_eachRet',...
    'cov_allRet','crr_eachOr','crr_eachRet')

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

%Plot firing rate
figure
hold on
errorbar(tme(start:finish)',psthMCor(1,:),EBpsthMCor(1,:),'b')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('PSTH (Hz)')
if(flag_Data==1)
    plot(tme(start:finish)',mean(mn_allOr(:,start:finish)./Twin),'b','LineWidth',2)
end
if(flag_Sep==1)
    figure
    hold on
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('PSTH (Hz)')
end
if(flag_Data==1)
    errorbar(tme(start:finish)',psthMCret(1,:),EBpsthMCret(1,:),'r')
    plot(tme(start:finish)',mean(mn_allRet(:,start:finish)./Twin),'r','LineWidth',2)
else
    errorbar(tme(start:finish)',psthMCret(1,:),EBpsthMCret(1,:),'r')
end

%Plot variances
figure
hold on
errorbar(tme(start:finish)',varMCor(1,:),EBvarMCor(1,:),'b')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Variance of Counts')
if(flag_Data==1)
    plot(tme(start:finish)',mean(vr_allOr(:,start:finish)),'b','LineWidth',2)
end
if(flag_Sep==1)
    figure
    hold on
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('Variance of Counts')
end
if(flag_Data==1)
    errorbar(tme(start:finish)',varMCret(1,:),EBvarMCret(1,:),'r')
    plot(tme(start:finish)',mean(vr_allRet(:,start:finish)),'r','LineWidth',2)
else
    errorbar(tme(start:finish)',varMCret(1,:),EBvarMCret(1,:),'r')
end

%Plot covariance
figure
hold on
errorbar(tme(start:finish)',covMCor,EBcovMCor,'b')
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Covariance of Counts')
if(flag_Data==1)
    plot(tme(start:finish)',mean(cov_allOr(:,start:finish)),'b','LineWidth',2)
end
if(flag_Sep==1)
    figure
    hold on
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('Covariance of Counts')
end
if(flag_Data==1)
    errorbar(tme(start:finish)',covMCret(1,:),EBcovMCret(1,:),'r')
    plot(tme(start:finish)',mean(cov_allRet(:,start:finish)),'r','LineWidth',2)
else
    errorbar(tme(start:finish)',covMCret(1,:),EBcovMCret(1,:),'r')
end

% plot ratio of Ret to Orth to show similar
mnCov_Od=mean(cov_allOr);
mnCov_Rd=mean(cov_allRet);
RetRat_dat=mnCov_Rd(31:finish)./mnCov_Od(31:finish);
RetRat_mod=covMCret(12:end)./covMCor(12:end);
figure
hold on
plot(1,mean(RetRat_dat),'k*')
plot(1,mean(RetRat_dat)+std(RetRat_dat),'ko')
plot(1,mean(RetRat_dat)-std(RetRat_dat),'ko')
plot(2,mean(RetRat_mod),'b*')
plot(2,mean(RetRat_mod)+std(RetRat_mod),'bo')
plot(2,mean(RetRat_mod)-std(RetRat_mod),'bo')
set(gca,'FontSize',18)
set(gca,'XLim',[0 3])

% plot FF
figure
hold on
plot(tme(start:finish)',varMCor(1,:)./(psthMCor(1,:)*Twin),'b','LineWidth',2)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('FF')
if(flag_Data==1)
    plot(tme(start:finish)',mean(fano_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',2)
    plot(tme(start:finish)',mean(fano_eachOr(:,start:finish))+.2*std(fano_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',.5)
    plot(tme(start:finish)',mean(fano_eachOr(:,start:finish))-.2*std(fano_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',.5)
end
if(flag_Sep==1)
    figure
    hold on
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('FF')
end
if(flag_Data==1)
    plot(tme(start:finish)',varMCret(1,:)./(psthMCret(1,:)*Twin),'r','LineWidth',2)
    plot(tme(start:finish)',mean(fano_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',2)
    plot(tme(start:finish)',mean(fano_eachRet(:,start:finish))+.2*std(fano_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',.5)
    plot(tme(start:finish)',mean(fano_eachRet(:,start:finish))-.2*std(fano_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',.5)
else
    plot(tme(start:finish)',varMCret(1,:)./(psthMCret(1,:)*Twin),'r','LineWidth',2)
end

% plot Correl
figure
hold on
plot(tme(start:finish)',covMCor(1,:)./varMCor(1,:),'b','LineWidth',2)
set(gca,'FontSize',18)
xlabel('Time (s)')
ylabel('Corr')
if(flag_Data==1)
    plot(tme(start:finish)',mean(crr_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',2)
    plot(tme(start:finish)',mean(crr_eachOr(:,start:finish))+.2*std(crr_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',.5)
    plot(tme(start:finish)',mean(crr_eachOr(:,start:finish))-.2*std(crr_eachOr(:,start:finish)),'color',[30 0 140]./255,'LineWidth',.5)
end
if(flag_Sep==1)
    figure
    hold on
    set(gca,'FontSize',18)
    xlabel('Time (s)')
    ylabel('Corr')
end
if(flag_Data==1)
    plot(tme(start:finish)',covMCret(1,:)./varMCret(1,:),'r','LineWidth',2)
    plot(tme(start:finish)',mean(crr_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',2)
    plot(tme(start:finish)',mean(crr_eachRet(:,start:finish))+.2*std(crr_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',.5)
    plot(tme(start:finish)',mean(crr_eachRet(:,start:finish))-.2*std(crr_eachRet(:,start:finish)),'color',[.5 0 0],'LineWidth',.5)
else
    plot(tme(start:finish)',covMCret(1,:)./varMCret(1,:),'r','LineWidth',2)
end

% plot ratio of Ret to Orth CORREL to show similar
mnCrr_Od=mean(crr_eachOr);
mnCrr_Rd=mean(crr_eachRet);
RetRat_dat=mnCrr_Rd(31:finish)./mnCrr_Od(31:finish);
RetRat_mod=covMCret(1,12:end)./varMCret(1,12:end)./(covMCor(1,12:end)./varMCor(1,12:end));
figure
hold on
plot(1,mean(RetRat_dat),'k*')
plot(1,mean(RetRat_dat)+std(RetRat_dat),'ko')
plot(1,mean(RetRat_dat)-std(RetRat_dat),'ko')
plot(2,mean(RetRat_mod),'b*')
plot(2,mean(RetRat_mod)+std(RetRat_mod),'bo')
plot(2,mean(RetRat_mod)-std(RetRat_mod),'bo')
set(gca,'FontSize',18)
set(gca,'XLim',[0 3])
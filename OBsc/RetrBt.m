function []=RetrBt()
% calls obscRt.m

nTrials=1;
indR=1;

lamA=struct('tsft',25.1,'lmAspon',.002,'lmAevk',.095/367,'tau_A',450);

timeVars=struct('dt',0.01,'tSpon',1200,'tEvok',1900);
tevok=(0:timeVars.dt:timeVars.tEvok)';
tspon=(-timeVars.tSpon:timeVars.dt:0)';

crG_R=[0.3*ones(length(tspon)-1,1);0.3*ones(length(tevok),1)];
crMP_R=[0.3*ones(length(tspon)-1,1);0.3*ones(length(tevok),1)];

crEvr=0.2+0.5*(1-exp(-(tevok(10001:end)-tevok(10001))/100));
crM_R=[0.2*ones(length(tspon)+9999,1);crEvr];

flName='dRet';

obscRt((indR-1)*nTrials+1,nTrials,flName,lamA,crG_R,crM_R,crMP_R,timeVars,tevok,tspon);

end
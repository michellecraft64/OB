function []=OrthBt()
% calls obscOt.m

nTrials=1;
indR=1;

lamO=struct('tsft',25.1,'lmOsp',.002,'lmOevk',15/367,'tauO',60,'tauO2',61,'ssV',.007);

timeVars=struct('dt',0.01,'tSpon',1200,'tEvok',1900);
tevok=(0:timeVars.dt:timeVars.tEvok)';
tspon=(-timeVars.tSpon:timeVars.dt:0)';

crG_O=[0.3*ones(length(tspon)-1,1);0.3*ones(length(tevok),1)];
crMP_O=[0.3*ones(length(tspon)-1,1);0.3*ones(length(tevok),1)];

crEv=0.18*(-(tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO) + (tevok+lamO.tsft).*exp(-(tevok+lamO.tsft)/lamO.tauO2))+.9;
crM_O=[0.2*ones(length(tspon)-1,1); crEv];
flName='dOrt';

obscOt((indR-1)*nTrials+1,nTrials,flName,lamO,crG_O,crM_O,crMP_O,timeVars,tevok,tspon);

end